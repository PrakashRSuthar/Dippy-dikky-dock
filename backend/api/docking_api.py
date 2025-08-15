# backend/api/docking_api.py

from fastapi import FastAPI, HTTPException, UploadFile, File, Form
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import Optional, List
import json
import uuid
from pathlib import Path
from datetime import datetime

from modules.protein_fetcher import ProteinFetcher
from modules.protein_preprocessor import prepare_protein
from modules.ligand_preparer import process_ligand, _convert_to_pdbqt
from modules.pocket_identifier import identify_binding_site
from modules.docking_engine import DockingEngine
from modules.result_parser import parse_vina_results
from modules.docking_accuracy_evaluator import evaluate_docking_accuracy_from_file

app = FastAPI(title="Dippy Docking API", version="1.0.0")

# CORS for frontend
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # In production, specify your frontend URL
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Data Models
class DockingRequest(BaseModel):
    protein_input: str
    ligand_input: str
    job_name: Optional[str] = None

class JobStatus(BaseModel):
    job_id: str
    status: str  # "queued", "running", "completed", "failed"
    progress: int  # 0-100
    message: str
    created_at: str
    completed_at: Optional[str] = None

class DockingResult(BaseModel):
    job_id: str
    status: str
    best_affinity: Optional[float] = None
    total_poses: Optional[int] = None
    pocket_center: Optional[List[float]] = None
    pocket_size: Optional[List[float]] = None
    method: Optional[str] = None
    confidence: Optional[str] = None
    run_directory: Optional[str] = None
    files: Optional[dict] = None
    error_message: Optional[str] = None

# Global job tracking
jobs = {}

class DockingPipeline:
    def __init__(self):
        self.status_callback = None
    
    def set_status_callback(self, callback):
        self.status_callback = callback
    
    def update_status(self, job_id: str, status: str, progress: int, message: str):
        if self.status_callback:
            self.status_callback(job_id, status, progress, message)
    
    def run_pipeline(self, job_id: str, protein_input: str, ligand_input: str):
        try:
            self.update_status(job_id, "running", 10, "Setting up workspace...")
            
            # Setup
            run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
            run_dir = Path(f"runs/{run_id}")
            raw_dir = run_dir / "raw"
            prepared_dir = run_dir / "prepared"
            
            for d in [raw_dir, prepared_dir]:
                d.mkdir(parents=True, exist_ok=True)

            self.update_status(job_id, "running", 20, "Resolving protein input...")
            
            # Resolve protein
            raw_protein_path = self.resolve_protein_input(protein_input, raw_dir)
            
            self.update_status(job_id, "running", 30, "Resolving ligand input...")
            
            # Resolve ligand
            prepared_ligand_pdbqt = self.resolve_ligand_input(ligand_input, raw_dir, prepared_dir)
            
            self.update_status(job_id, "running", 40, "Preparing protein receptor...")
            
            # Prepare protein
            prepared_protein_pdbqt = prepare_protein(raw_protein_path)
            if not prepared_protein_pdbqt:
                raise Exception("Protein preparation failed")
            
            if not prepared_protein_pdbqt.endswith('.pdbqt'):
                raise Exception("Protein preparation produced invalid format")

            self.update_status(job_id, "running", 60, "Detecting binding pocket...")
            
            # Pocket detection
            pocket_info = identify_binding_site(
                raw_protein_path, prepared_protein_pdbqt, prepared_ligand_pdbqt,
                use_validated=True, return_n=1
            )
            if not pocket_info:
                raise Exception("Pocket detection failed")
            
            center = pocket_info["primary"]

            self.update_status(job_id, "running", 80, "Running molecular docking...")
            
            # Docking
            docking_dir = run_dir / "docking"
            docking_dir.mkdir(exist_ok=True)
            
            engine = DockingEngine(prepared_protein_pdbqt, prepared_ligand_pdbqt)
            docked_file = engine.run_docking(
                center_x=center["center_x"], center_y=center["center_y"], center_z=center["center_z"],
                box_size_x=center["size_x"], box_size_y=center["size_y"], box_size_z=center["size_z"],
                output_dir=str(docking_dir)
            )
            
            if not docked_file:
                raise Exception("Docking simulation failed")

            self.update_status(job_id, "running", 90, "Parsing results...")
            
            # Parse results
            results_csv = docking_dir / "binding_scores.csv"
            df = parse_vina_results(docked_file, output_csv_path=results_csv)
            if df is None or df.empty:
                raise Exception("Results parsing failed")
            
            best_affinity = float(df["Binding Affinity (kcal/mol)"].min())

            self.update_status(job_id, "running", 95, "Finalizing results...")
            
            # Create result summary
            summary = {
                "job_id": job_id,
                "status": "completed",
                "run_dir": str(run_dir.resolve()),
                "best_affinity": best_affinity,
                "total_poses": len(df),
                "pocket_method": center["method"],
                "pocket_center": [center["center_x"], center["center_y"], center["center_z"]],
                "pocket_size": [center["size_x"], center["size_y"], center["size_z"]],
                "confidence": center["confidence"],
                "docked_file": docked_file,
                "results_csv": str(results_csv),
                "files": {
                    "raw_protein": raw_protein_path,
                    "prepared_protein": prepared_protein_pdbqt,
                    "prepared_ligand": prepared_ligand_pdbqt,
                    "docked_poses": docked_file,
                    "results_csv": str(results_csv)
                }
            }

            (run_dir / "results.json").write_text(json.dumps(summary, indent=2, default=str))
            
            self.update_status(job_id, "completed", 100, f"Docking completed! Best affinity: {best_affinity:.2f} kcal/mol")
            
            return summary

        except Exception as e:
            error_msg = f"Pipeline failed: {str(e)}"
            self.update_status(job_id, "failed", 0, error_msg)
            return {"job_id": job_id, "status": "failed", "error_message": error_msg}

    def resolve_protein_input(self, user_input: str, raw_dir: Path) -> str:
        """Resolve protein input: file path or ID"""
        p = Path(user_input)
        if p.exists() and p.is_file():
            dst = raw_dir / p.name
            if dst.resolve() != p.resolve():
                dst.write_bytes(p.read_bytes())
            return str(dst)

        # Fetch from database
        fetcher = ProteinFetcher(download_dir=str(raw_dir))
        if len(user_input) == 4 and user_input.isalnum():
            return fetcher.fetch_from_pdb(user_input)
        else:
            return fetcher.fetch_from_alphafold(user_input)

    def resolve_ligand_input(self, user_input: str, raw_dir: Path, prepared_dir: Path) -> str:
        """Resolve ligand input: file path or PubChem name"""
        p = Path(user_input)
        if p.exists() and p.is_file():
            dst = raw_dir / p.name
            if dst.resolve() != p.resolve():
                dst.write_bytes(p.read_bytes())

            ext = p.suffix.lower()
            if ext == ".pdbqt":
                return str(dst)
            elif ext in {".sdf", ".mol2", ".pdb"}:
                return _convert_to_pdbqt(str(dst), str(prepared_dir))

        # Fetch from PubChem
        return process_ligand(user_input, ligand_dir=str(raw_dir), prepared_dir=str(prepared_dir))

def update_job_status(job_id: str, status: str, progress: int, message: str):
    """Update job status in global storage"""
    if job_id in jobs:
        jobs[job_id].status = status
        jobs[job_id].progress = progress
        jobs[job_id].message = message
        if status in ["completed", "failed"]:
            jobs[job_id].completed_at = datetime.now().isoformat()

# API Endpoints
@app.post("/api/dock", response_model=dict)
async def start_docking(request: DockingRequest):
    """Start a new docking job"""
    job_id = str(uuid.uuid4())
    
    # Create job entry
    jobs[job_id] = JobStatus(
        job_id=job_id,
        status="queued",
        progress=0,
        message="Job queued for processing",
        created_at=datetime.now().isoformat()
    )
    
    # Start pipeline in background (in production, use Celery or similar)
    import threading
    pipeline = DockingPipeline()
    pipeline.set_status_callback(update_job_status)
    
    thread = threading.Thread(
        target=pipeline.run_pipeline,
        args=(job_id, request.protein_input, request.ligand_input)
    )
    thread.daemon = True
    thread.start()
    
    return {"job_id": job_id, "message": "Docking job started", "status": "queued"}

@app.get("/api/jobs/{job_id}/status", response_model=JobStatus)
async def get_job_status(job_id: str):
    """Get status of a specific job"""
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    return jobs[job_id]

@app.get("/api/jobs/{job_id}/result", response_model=DockingResult)
async def get_job_result(job_id: str):
    """Get detailed results of a completed job"""
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = jobs[job_id]
    if job.status != "completed":
        raise HTTPException(status_code=400, detail="Job not completed yet")
    
    # Load detailed results
    try:
        run_dirs = list(Path("runs").glob("*/"))
        for run_dir in run_dirs:
            results_file = run_dir / "results.json"
            if results_file.exists():
                with open(results_file) as f:
                    data = json.load(f)
                    if data.get("job_id") == job_id:
                        return DockingResult(**data)
        
        raise HTTPException(status_code=404, detail="Job results not found")
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error loading results: {e}")

@app.get("/api/jobs", response_model=List[JobStatus])
async def list_jobs():
    """List all jobs"""
    return list(jobs.values())

@app.post("/api/upload/protein")
async def upload_protein(file: UploadFile = File(...)):
    """Upload protein file"""
    upload_dir = Path("uploads/proteins")
    upload_dir.mkdir(parents=True, exist_ok=True)
    
    file_path = upload_dir / file.filename
    with open(file_path, "wb") as f:
        content = await file.read()
        f.write(content)
    
    return {"message": "Protein file uploaded", "file_path": str(file_path)}

@app.post("/api/upload/ligand")
async def upload_ligand(file: UploadFile = File(...)):
    """Upload ligand file"""
    upload_dir = Path("uploads/ligands")
    upload_dir.mkdir(parents=True, exist_ok=True)
    
    file_path = upload_dir / file.filename
    with open(file_path, "wb") as f:
        content = await file.read()
        f.write(content)
    
    return {"message": "Ligand file uploaded", "file_path": str(file_path)}

@app.get("/")
async def root():
    return {"message": "Dippy Docking API", "version": "1.0.0", "status": "running"}

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
