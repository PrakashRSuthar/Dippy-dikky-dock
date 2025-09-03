# backend/api/docking_api.py

import sys
from pathlib import Path
import traceback
import shutil

# Fix import path
backend_dir = Path(__file__).parent.parent
sys.path.insert(0, str(backend_dir))

print(f"✅ Added to Python path: {backend_dir}")

from fastapi import FastAPI, HTTPException, UploadFile, File, Form, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import StreamingResponse, FileResponse
from pydantic import BaseModel
from typing import Optional, List
import json
import uuid
from datetime import datetime
import threading
from contextlib import asynccontextmanager
import asyncio

# Import your docking modules
try:
    from modules.protein_fetcher import ProteinFetcher
    from modules.protein_preprocessor import prepare_protein, ProteinPreprocessor
    from modules.ligand_preparer import process_ligand, _convert_to_pdbqt
    from modules.pocket_identifier import identify_binding_site
    from modules.docking_engine import DockingEngine
    from modules.result_parser import parse_vina_results
    print("✅ All modules imported successfully")
except ImportError as e:
    print(f"❌ Module import error: {e}")
    sys.exit(1)

@asynccontextmanager
async def lifespan(app: FastAPI):
    print("🚀 Dippy Docking API Started Successfully!")
    yield
    print("🛑 Shutting down Dippy Docking API")

app = FastAPI(title="Dippy Docking API", version="1.0.0", lifespan=lifespan)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Data Models
class DockingRequest(BaseModel):
    protein_input: str
    ligand_input: str
    job_name: Optional[str] = None
    cleaning_policy: Optional[dict] = {
        "keep_waters": False,
        "keep_ions": True,
        "keep_solvents": False,
        "keep_cofactors": True
    }
    retention: Optional[str] = "keep7d"

class JobStatus(BaseModel):
    job_id: str
    status: str
    progress: int
    message: str
    created_at: str
    completed_at: Optional[str] = None

# Global job tracking
jobs = {}
job_logs = {}

class CompleteFileTracker:
    """File tracker from your original pipeline"""
    def __init__(self, workspace: Path):
        self.workspace = workspace
        self.external_files: List[Path] = []
        self.data_dirs_to_check = [
            Path("data/prepared_proteins"),
            Path("data/prepared_ligands"), 
            Path("data/proteins"),
            Path("data/ligands")
        ]
    
    def find_external_files_by_timestamp(self, timestamp: str):
        import re
        pattern = re.compile(rf".*{timestamp}.*\.(pdb|pdbqt|json|sdf|mol2)$", re.IGNORECASE)
        for data_dir in self.data_dirs_to_check:
            if data_dir.exists():
                for file_path in data_dir.rglob("*"):
                    if file_path.is_file() and pattern.match(file_path.name):
                        if not self.is_in_workspace(file_path):
                            self.external_files.append(file_path)
    
    def is_in_workspace(self, file_path: Path) -> bool:
        try:
            file_path.relative_to(self.workspace)
            return True
        except ValueError:
            return False
    
    def save_permanently(self, project_name: str, summary_data: dict) -> bool:
        try:
            permanent_path = Path(f"data/results/{project_name}")
            permanent_path.mkdir(parents=True, exist_ok=True)
            
            # Copy workspace
            workspace_dest = permanent_path / "workspace"
            shutil.copytree(self.workspace, workspace_dest)
            
            # Copy external files
            external_dest = permanent_path / "external_files"
            external_dest.mkdir(exist_ok=True)
            for ext_file in self.external_files:
                if ext_file.exists():
                    dest_file = external_dest / ext_file.name
                    shutil.copy2(ext_file, dest_file)
            
            # Save summary
            summary_data['permanent_location'] = str(permanent_path)
            summary_file = permanent_path / "summary.json"
            with open(summary_file, 'w') as f:
                json.dump(summary_data, f, indent=2, default=str)
            
            return True
        except Exception as e:
            print(f"❌ Error saving permanently: {e}")
            return False
    
    def cleanup_all_tracked_files(self):
        for file_path in self.external_files:
            try:
                if file_path.exists():
                    file_path.unlink()
            except Exception as e:
                print(f"Warning: Could not remove {file_path}: {e}")
        
        if self.workspace.exists():
            shutil.rmtree(self.workspace)

class DockingPipeline:
    def __init__(self):
        self.status_callback = None
        self.log_callback = None
    
    def set_status_callback(self, callback):
        self.status_callback = callback
    
    def set_log_callback(self, callback):
        self.log_callback = callback
    
    def update_progress(self, job_id: str, percent: int, message: str):
        bar_length = 30
        filled_length = int(bar_length * percent // 100)
        bar = '█' * filled_length + '░' * (bar_length - filled_length)
        
        progress_msg = f"Progress: |{bar}| {percent}% - {message}"
        print(f"\r{progress_msg}", end='', flush=True)
        
        if self.status_callback:
            self.status_callback(job_id, "running", percent, message)
        
        self.log_message(job_id, f"[{percent}%] {message}", "INFO")
    
    def log_message(self, job_id: str, message: str, level: str = "INFO"):
        timestamp = datetime.now().strftime("%H:%M:%S")
        log_entry = {
            "timestamp": timestamp,
            "level": level,
            "message": message
        }
        
        if job_id not in job_logs:
            job_logs[job_id] = []
        job_logs[job_id].append(log_entry)
        
        level_emoji = {"INFO": "ℹ️", "SUCCESS": "✅", "ERROR": "❌", "WARNING": "⚠️"}
        print(f"\n{level_emoji.get(level, 'ℹ️')} [{timestamp}] {message}")
        
        if self.log_callback:
            self.log_callback(job_id, log_entry)

    # ⭐ MISSING METHODS - Adding them here ⭐
    def resolve_protein_input(self, user_input: str, raw_dir: Path, job_id: str) -> str:
        """Resolve protein input: file path or ID"""
        if not user_input or user_input.strip() == "":
            raise Exception("Empty protein input")
            
        user_input = user_input.strip()
        p = Path(user_input)
        
        if p.exists() and p.is_file():
            dst = raw_dir / p.name
            if dst.resolve() != p.resolve():
                dst.write_bytes(p.read_bytes())
            return str(dst)

        # Fetch from database
        fetcher = ProteinFetcher(download_dir=str(raw_dir))
        if len(user_input) == 4 and user_input.isalnum():
            self.log_message(job_id, f"🧬 Fetching protein from PDB: {user_input}", "INFO")
            return fetcher.fetch_from_pdb(user_input)
        else:
            self.log_message(job_id, f"🧬 Fetching protein from AlphaFold: {user_input}", "INFO")
            return fetcher.fetch_from_alphafold(user_input)

    def resolve_ligand_input(self, user_input: str, raw_dir: Path, prepared_dir: Path, job_id: str) -> str:
        """Resolve ligand input: file path or PubChem name"""
        if not user_input or user_input.strip() == "":
            raise Exception("Empty ligand input")
            
        user_input = user_input.strip()
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
        self.log_message(job_id, f"💊 Fetching ligand from PubChem: {user_input}", "INFO")
        return process_ligand(user_input, ligand_dir=str(raw_dir), prepared_dir=str(prepared_dir))
    
    def prepare_protein_in_workspace(self, protein_pdb_path: str, workspace_prepared_dir: Path, 
                                   run_timestamp: str, cleaning_policy: dict, job_id: str) -> str:
        """Prepare protein with output directed to pipeline workspace"""
        workspace_prepared_dir.mkdir(parents=True, exist_ok=True)
        
        # Create a custom protein preprocessor that outputs to workspace
        processor = ProteinPreprocessor(Path(protein_pdb_path))
        
        # Override the output directory to use workspace
        processor.output_dir = workspace_prepared_dir
        
        # Update file paths to use workspace
        processor.base_stem = f"{Path(protein_pdb_path).stem}_{run_timestamp}"
        processor.tmp_cleaned_pdb = workspace_prepared_dir / f"{processor.base_stem}_clean.pdb"
        processor.clean_report = workspace_prepared_dir / f"{processor.base_stem}_clean_report.json"
        processor.output_pdbqt = workspace_prepared_dir / f"{processor.base_stem}_prepared.pdbqt"
        
        # Process with workspace output
        result = processor.process(interactive=False, default_policy=cleaning_policy)
        
        return result
    
    def display_binding_sites(self, pocket_info: dict, job_id: str):
        """Display top 5 binding sites"""
        self.log_message(job_id, "🎯 TOP 5 BINDING SITES:", "INFO")
        self.log_message(job_id, "="*60, "INFO")
        
        all_sites = []
        if "primary" in pocket_info:
            all_sites.append(pocket_info["primary"])
        if "modes" in pocket_info:
            all_sites.extend(pocket_info["modes"])
        
        for i, site in enumerate(all_sites[:5], 1):
            site_msg = f"Site {i}: ({site.get('center_x', 0):.1f}, {site.get('center_y', 0):.1f}, {site.get('center_z', 0):.1f}) Box: {site.get('size_x', 20):.0f}×{site.get('size_y', 20):.0f}×{site.get('size_z', 20):.0f}Å"
            self.log_message(job_id, site_msg, "INFO")
        
        if len(all_sites) > 5:
            self.log_message(job_id, f"... and {len(all_sites) - 5} more sites", "INFO")
        self.log_message(job_id, "="*60, "INFO")

    def handle_data_retention(self, job_id: str, retention: str, file_tracker: CompleteFileTracker, summary_data: dict, run_id: str):
        """Handle data retention based on user choice"""
        self.log_message(job_id, "💾 DATA MANAGEMENT", "INFO")
        
        # Find external files
        file_tracker.find_external_files_by_timestamp(run_id)
        
        if retention == "save":
            self.log_message(job_id, "💾 Saving data permanently...", "INFO")
            project_name = f"docking_{job_id[:8]}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
            if file_tracker.save_permanently(project_name, summary_data):
                self.log_message(job_id, f"✅ Data saved permanently to data/results/{project_name}", "SUCCESS")
            else:
                self.log_message(job_id, "❌ Failed to save data permanently", "ERROR")
                
        elif retention == "delete":
            self.log_message(job_id, "🗑️ Deleting all data...", "INFO")
            file_tracker.cleanup_all_tracked_files()
            self.log_message(job_id, "✅ All data deleted", "SUCCESS")
            
        elif retention == "keep7d":
            self.log_message(job_id, "⏳ Keeping data temporarily (7 days)", "INFO")
            self.log_message(job_id, "🗑️ Will be auto-deleted in 7 days", "INFO")

    def run_pipeline(self, job_id: str, protein_input: str, ligand_input: str, cleaning_policy: dict, retention: str):
        try:
            self.log_message(job_id, f"🚀 Docking Pipeline [{job_id[:8]}]", "INFO")
            self.log_message(job_id, f"   Protein: {protein_input}", "INFO")
            self.log_message(job_id, f"   Ligand: {ligand_input}", "INFO")
            self.log_message(job_id, f"   Retention Policy: {retention}", "INFO")
            
            # Setup workspace
            run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
            workspace = Path(f"temp_runs/{run_id}_{job_id[:8]}")
            run_dir = workspace / "run"
            raw_dir = run_dir / "raw"
            prepared_dir = run_dir / "prepared"
            
            for d in [workspace, run_dir, raw_dir, prepared_dir]:
                d.mkdir(parents=True, exist_ok=True)
            
            # Initialize file tracker
            file_tracker = CompleteFileTracker(workspace)
            
            self.update_progress(job_id, 15, "Resolving protein")
            raw_protein_path = self.resolve_protein_input(protein_input, raw_dir, job_id)
            
            self.update_progress(job_id, 25, "Resolving ligand")
            prepared_ligand_pdbqt = self.resolve_ligand_input(ligand_input, raw_dir, prepared_dir, job_id)
            
            self.update_progress(job_id, 40, "Analyzing protein structure")
            analyzer = ProteinPreprocessor(raw_protein_path)
            summary = analyzer.summarize_pdb(analyzer.pdb_path)
            
            analysis_msg = f"Found: {summary['waters']} waters, {summary['ions']} ions, {summary['solvents']} solvents, {summary['cofactors']} cofactors"
            self.update_progress(job_id, 42, f"Protein analysis: {analysis_msg}")
            
            self.update_progress(job_id, 45, "Preparing protein (workspace-only)")
            prepared_protein_pdbqt = self.prepare_protein_in_workspace(
                raw_protein_path, prepared_dir, run_id, cleaning_policy, job_id
            )
            
            if not prepared_protein_pdbqt:
                raise Exception("Protein preparation failed")
            
            self.update_progress(job_id, 55, "Finding binding sites")
            pocket_info = identify_binding_site(
                raw_protein_path, prepared_protein_pdbqt, prepared_ligand_pdbqt,
                use_validated=True, return_n=5, detailed=True
            )
            if not pocket_info:
                raise Exception("Pocket detection failed")
            
            center = pocket_info["primary"]
            self.display_binding_sites(pocket_info, job_id)

            self.update_progress(job_id, 65, "Running molecular docking")
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

            self.update_progress(job_id, 85, "Analyzing results")
            results_csv = docking_dir / "binding_scores.csv"
            df = parse_vina_results(docked_file, output_csv_path=results_csv, pocket_info=pocket_info)
            if df is None or df.empty:
                raise Exception("Results parsing failed")
            
            best_affinity = float(df["Binding Affinity (kcal/mol)"].min())
            total_poses = len(df)

            # Extract all pocket information
            all_pockets = []
            all_sites = []
            if "primary" in pocket_info:
                all_sites.append(pocket_info["primary"])
            if "modes" in pocket_info:
                all_sites.extend(pocket_info["modes"])
            
            for site in all_sites:
                all_pockets.append({
                    "center": [site.get("center_x", 0), site.get("center_y", 0), site.get("center_z", 0)],
                    "size": [site.get("size_x", 20), site.get("size_y", 20), site.get("size_z", 20)],
                    "confidence": site.get("confidence", "unknown"),
                    "method": site.get("method", "unknown")
                })

            self.update_progress(job_id, 95, "Generating summary")
            
            # Create comprehensive result summary
            summary = {
                "job_id": job_id,
                "status": "completed",
                "run_dir": str(run_dir.resolve()),
                "best_affinity": best_affinity,
                "total_poses": total_poses,
                "pocket_method": center["method"],
                "pocket_center": [center["center_x"], center["center_y"], center["center_z"]],
                "pocket_size": [center["size_x"], center["size_y"], center["size_z"]],
                "confidence": center["confidence"],
                "docked_file": docked_file,
                "results_csv": str(results_csv),
                "cleaning_policy": cleaning_policy,
                "protein_analysis": analysis_msg,
                "all_pockets": all_pockets,
                "files": {
                    "raw_protein": raw_protein_path,
                    "prepared_protein": prepared_protein_pdbqt,
                    "prepared_ligand": prepared_ligand_pdbqt,
                    "docked_poses": docked_file,
                    "results_csv": str(results_csv)
                },
                "retention": retention,
                "workspace": str(workspace)
            }

            # Save results
            (run_dir / "results.json").write_text(json.dumps(summary, indent=2, default=str))
            
            self.update_progress(job_id, 100, "Pipeline completed")
            
            # Handle data retention
            self.handle_data_retention(job_id, retention, file_tracker, summary, run_id)
            
            self.log_message(job_id, f"🎉 DOCKING COMPLETED!", "SUCCESS")
            self.log_message(job_id, f"   Best Affinity: {best_affinity:.2f} kcal/mol", "SUCCESS")
            self.log_message(job_id, f"   Total Poses: {total_poses}", "SUCCESS")
            self.log_message(job_id, f"   All Pockets Found: {len(all_pockets)}", "SUCCESS")
            
            if self.status_callback:
                self.status_callback(job_id, "completed", 100, f"Docking completed! Best affinity: {best_affinity:.2f} kcal/mol")
            
            return summary

        except Exception as e:
            error_msg = f"Pipeline failed: {str(e)}"
            self.log_message(job_id, f"❌ {error_msg}", "ERROR")
            print(f"\n{traceback.format_exc()}")
            if self.status_callback:
                self.status_callback(job_id, "failed", 0, error_msg)
            return {"job_id": job_id, "status": "failed", "error_message": error_msg}

def update_job_status(job_id: str, status: str, progress: int, message: str):
    if job_id in jobs:
        jobs[job_id].status = status
        jobs[job_id].progress = progress
        jobs[job_id].message = message
        if status in ["completed", "failed"]:
            jobs[job_id].completed_at = datetime.now().isoformat()

# API Endpoints
@app.post("/api/dock", response_model=dict)
async def start_docking(request: DockingRequest):
    job_id = str(uuid.uuid4())
    
    cleaning_policy = request.cleaning_policy or {
        "keep_waters": False,
        "keep_ions": True,
        "keep_solvents": False,
        "keep_cofactors": True
    }
    
    jobs[job_id] = JobStatus(
        job_id=job_id,
        status="queued",
        progress=0,
        message="Job queued for processing",
        created_at=datetime.now().isoformat()
    )
    
    job_logs[job_id] = []
    
    print(f"🚀 Starting docking job {job_id}")
    print(f"   Retention: {request.retention}")
    
    def run_job():
        pipeline = DockingPipeline()
        pipeline.set_status_callback(update_job_status)
        pipeline.run_pipeline(job_id, request.protein_input, request.ligand_input, cleaning_policy, request.retention)
    
    thread = threading.Thread(target=run_job, daemon=True)
    thread.start()
    
    return {"job_id": job_id, "message": "Docking job started", "status": "queued"}

@app.get("/api/jobs/{job_id}/status", response_model=JobStatus)
async def get_job_status(job_id: str):
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    return jobs[job_id]

@app.get("/api/jobs/{job_id}/logs")
async def stream_job_logs(job_id: str):
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    async def log_generator():
        sent_count = 0
        while True:
            current_logs = job_logs.get(job_id, [])
            
            for i in range(sent_count, len(current_logs)):
                log_entry = current_logs[i]
                yield f"data: {json.dumps(log_entry)}\n\n"
            
            sent_count = len(current_logs)
            
            job_status = jobs.get(job_id, {})
            if getattr(job_status, 'status', None) in ['completed', 'failed']:
                break
            
            await asyncio.sleep(0.5)
    
    return StreamingResponse(
        log_generator(),
        media_type="text/event-stream",
        headers={"Cache-Control": "no-cache", "Connection": "keep-alive"}
    )

@app.get("/api/jobs/{job_id}/result")
async def get_job_result(job_id: str):
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = jobs[job_id]
    if job.status != "completed":
        raise HTTPException(status_code=400, detail="Job not completed yet")
    
    try:
        search_paths = [
            Path("temp_runs").glob("**/results.json"),
            Path("data/results").glob("**/results.json")
        ]
        
        for pattern in search_paths:
            for results_file in pattern:
                if results_file.exists():
                    with open(results_file) as f:
                        data = json.load(f)
                        if data.get("job_id") == job_id:
                            return data
        
        raise HTTPException(status_code=404, detail="Job results not found")
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error loading results: {e}")

@app.get("/api/jobs/{job_id}/files")
async def list_job_files(job_id: str):
    try:
        result_data = await get_job_result(job_id)
        files = result_data.get("files", {})
        
        file_list = []
        for name, path in files.items():
            file_path = Path(path)
            if file_path.exists():
                file_list.append({
                    "name": name.replace('_', ' ').title(),
                    "filename": file_path.name,
                    "path": str(path),
                    "size": f"{file_path.stat().st_size / 1024:.1f} KB",
                    "download_url": f"/api/download?path={path}"
                })
        
        return {"files": file_list}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error listing files: {e}")

@app.get("/api/download")
async def download_file(path: str = Query(..., description="File path to download")):
    file_path = Path(path)
    
    if not file_path.exists() or not file_path.is_file():
        raise HTTPException(status_code=404, detail=f"File not found: {path}")
    
    return FileResponse(
        str(file_path),
        filename=file_path.name,
        media_type='application/octet-stream'
    )

@app.get("/api/jobs")
async def list_jobs():
    return list(jobs.values())

@app.post("/api/upload/protein")
async def upload_protein(file: UploadFile = File(...)):
    upload_dir = Path("uploads/proteins")
    upload_dir.mkdir(parents=True, exist_ok=True)
    
    file_path = upload_dir / file.filename
    with open(file_path, "wb") as f:
        content = await file.read()
        f.write(content)
    
    return {"message": "Protein file uploaded", "file_path": str(file_path)}

@app.post("/api/upload/ligand")
async def upload_ligand(file: UploadFile = File(...)):
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
    uvicorn.run("api.docking_api:app", host="0.0.0.0", port=8000, reload=True)
