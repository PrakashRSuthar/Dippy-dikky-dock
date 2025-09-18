# backend/api/docking_api.py
import sys
from pathlib import Path
import traceback
import shutil

# Fix import path
backend_dir = Path(__file__).parent.parent
sys.path.insert(0, str(backend_dir))
print(f"✅ Added to Python path: {backend_dir}")

from fastapi import FastAPI, HTTPException, UploadFile, File, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, StreamingResponse
from pydantic import BaseModel
from typing import Optional, List, Dict, Any
import json
import uuid
from datetime import datetime
import threading
from contextlib import asynccontextmanager
import asyncio

# ========= Pipeline imports =========
try:
    from modules.protein_fetcher import ProteinFetcher
    from modules.protein_preprocessor import prepare_protein, ProteinPreprocessor
    from modules.ligand_preparer import process_ligand, _convert_to_pdbqt
    from modules.pocket_identifier import identify_binding_site
    from modules.docking_engine import DockingEngine
    from modules.result_parser import parse_vina_results
    print("🧩 Modules: ✅ Loaded all pipeline modules")
except ImportError as e:
    print(f"❌ Module import error: {e}")
    sys.exit(1)

# ========= App =========
@asynccontextmanager
async def lifespan(app: FastAPI):
    print("🚀 Dippy Docking API booting...")
    yield
    print("🛑 Dippy Docking API shutting down")

app = FastAPI(title="Dippy Docking API", version="1.0.0", lifespan=lifespan)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
print("🌐 CORS: ✅ All origins allowed")

# ========= Models =========
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

# ========= State =========
jobs: Dict[str, JobStatus] = {}
job_logs: Dict[str, List[Dict[str, Any]]] = {}

# ========= Logging helpers =========
def log_console(level: str, msg: str):
    prefix = {"INFO":"ℹ️","SUCCESS":"✅","ERROR":"❌","WARNING":"⚠️"}.get(level,"ℹ️")
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"{prefix} [{ts}] {msg}")

def progress_bar(pct: int, msg: str):
    bar_len = 30
    filled = int(bar_len * pct // 100)
    bar = '█' * filled + '░' * (bar_len - filled)
    print(f"\r⏳ |{bar}| {pct:3d}% - {msg}", end='', flush=True)

# ========= File tracking =========
class CompleteFileTracker:
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
        pattern = re.compile(rf".*{timestamp}.*\.(pdb|pdbqt|json|sdf|mol2|csv)$", re.IGNORECASE)
        for data_dir in self.data_dirs_to_check:
            if data_dir.exists():
                for file_path in data_dir.rglob("*"):
                    if file_path.is_file() and pattern.match(file_path.name):
                        if not self.is_in_workspace(file_path):
                            self.external_files.append(file_path)
        log_console("INFO", f"📦 Tracked external files: {len(self.external_files)}")

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
            workspace_dest = permanent_path / "workspace"
            shutil.copytree(self.workspace, workspace_dest)
            external_dest = permanent_path / "external_files"
            external_dest.mkdir(exist_ok=True)
            for ext_file in self.external_files:
                if ext_file.exists():
                    shutil.copy2(ext_file, external_dest / ext_file.name)
            summary_data['permanent_location'] = str(permanent_path)
            with open(permanent_path / "summary.json", "w") as f:
                json.dump(summary_data, f, indent=2, default=str)
            log_console("SUCCESS", f"💾 Saved permanently to {permanent_path}")
            return True
        except Exception as e:
            log_console("ERROR", f"Save permanently failed: {e}")
            return False

    def cleanup_all_tracked_files(self):
        for file_path in self.external_files:
            try:
                if file_path.exists():
                    file_path.unlink()
            except Exception as e:
                log_console("WARNING", f"Cleanup skip {file_path}: {e}")
        if self.workspace.exists():
            shutil.rmtree(self.workspace)
        log_console("SUCCESS", "🧹 Cleanup done")

# ========= Pipeline =========
class DockingPipeline:
    def __init__(self):
        self.status_callback = None
        self.log_callback = None

    def set_status_callback(self, callback):
        self.status_callback = callback

    def set_log_callback(self, callback):
        self.log_callback = callback

    def _status(self, job_id: str, pct: int, message: str):
        progress_bar(pct, message)
        if self.status_callback:
            self.status_callback(job_id, "running", pct, message)
        self._log(job_id, f"[{pct}%] {message}", "INFO")

    def _log(self, job_id: str, message: str, level: str = "INFO"):
        ts = datetime.now().strftime("%H:%M:%S")
        entry = {"timestamp": ts, "level": level, "message": message}
        job_logs.setdefault(job_id, []).append(entry)
        log_console(level, message)
        if self.log_callback:
            self.log_callback(job_id, entry)

    # ----- Input resolution -----
    def resolve_protein_input(self, user_input: str, raw_dir: Path, job_id: str) -> str:
        if not user_input or user_input.strip() == "":
            raise Exception("Empty protein input")
        user_input = user_input.strip()
        p = Path(user_input)
        if p.exists() and p.is_file():
            dst = raw_dir / p.name
            if dst.resolve() != p.resolve():
                dst.write_bytes(p.read_bytes())
            self._log(job_id, f"📥 Protein from local path -> {dst}", "INFO")
            return str(dst)
        fetcher = ProteinFetcher(download_dir=str(raw_dir))
        if len(user_input) == 4 and user_input.isalnum():
            self._log(job_id, f"🧬 Fetching PDB ID: {user_input}", "INFO")
            return fetcher.fetch_from_pdb(user_input)
        else:
            self._log(job_id, f"🧬 Fetching AlphaFold ID: {user_input}", "INFO")
            return fetcher.fetch_from_alphafold(user_input)

    def resolve_ligand_input(self, user_input: str, raw_dir: Path, prepared_dir: Path, job_id: str) -> str:
        if not user_input or user_input.strip() == "":
            raise Exception("Empty ligand input")
        user_input = user_input.strip()
        p = Path(user_input)
        if p.exists() and p.is_file():
            dst = raw_dir / p.name
            if dst.resolve() != p.resolve():
                dst.write_bytes(p.read_bytes())
            ext = p.suffix.lower()
            self._log(job_id, f"📥 Ligand from local path ({ext}) -> {dst}", "INFO")
            if ext == ".pdbqt":
                return str(dst)
            elif ext in {".sdf", ".mol2", ".pdb"}:
                out = _convert_to_pdbqt(str(dst), str(prepared_dir))
                self._log(job_id, f"🔁 Converted ligand to PDBQT -> {out}", "INFO")
                return out
        self._log(job_id, f"💊 Fetching ligand from PubChem: {user_input}", "INFO")
        return process_ligand(user_input, ligand_dir=str(raw_dir), prepared_dir=str(prepared_dir))

    def prepare_protein_in_workspace(self, protein_pdb_path: str, workspace_prepared_dir: Path,
                                     run_timestamp: str, cleaning_policy: dict, job_id: str) -> str:
        workspace_prepared_dir.mkdir(parents=True, exist_ok=True)
        processor = ProteinPreprocessor(Path(protein_pdb_path))
        processor.output_dir = workspace_prepared_dir
        processor.base_stem = f"{Path(protein_pdb_path).stem}_{run_timestamp}"
        processor.tmp_cleaned_pdb = workspace_prepared_dir / f"{processor.base_stem}_clean.pdb"
        processor.clean_report = workspace_prepared_dir / f"{processor.base_stem}_clean_report.json"
        processor.output_pdbqt = workspace_prepared_dir / f"{processor.base_stem}_prepared.pdbqt"
        self._log(job_id, "🧹 Cleaning & preparing protein...", "INFO")
        result = processor.process(interactive=False, default_policy=cleaning_policy)
        self._log(job_id, f"✅ Protein prepared: {result}", "SUCCESS")
        return result

    def run_pipeline(self, job_id: str, protein_input: str, ligand_input: str, cleaning_policy: dict, retention: str):
        try:
            self._log(job_id, f"🚀 Starting Docking Pipeline [{job_id[:8]}]", "INFO")
            self._log(job_id, f"• Protein: {protein_input}", "INFO")
            self._log(job_id, f"• Ligand: {ligand_input}", "INFO")
            self._log(job_id, f"• Retention: {retention}", "INFO")

            # Workspace
            run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
            workspace = Path(f"temp_runs/{run_id}_{job_id[:8]}")
            run_dir = workspace / "run"
            raw_dir = run_dir / "raw"
            prepared_dir = run_dir / "prepared"
            for d in [workspace, run_dir, raw_dir, prepared_dir]:
                d.mkdir(parents=True, exist_ok=True)
            self._log(job_id, f"📂 Workspace: {workspace}", "INFO")

            # Resolve inputs
            self._status(job_id, 15, "Resolving protein")
            raw_protein_path = self.resolve_protein_input(protein_input, raw_dir, job_id)

            self._status(job_id, 25, "Resolving ligand")
            prepared_ligand_pdbqt = self.resolve_ligand_input(ligand_input, raw_dir, prepared_dir, job_id)

            # Protein analysis
            self._status(job_id, 40, "Analyzing protein")
            analyzer = ProteinPreprocessor(raw_protein_path)
            summary = analyzer.summarize_pdb(analyzer.pdb_path)
            analysis_msg = f"Waters={summary['waters']}, Ions={summary['ions']}, Solvents={summary['solvents']}, Cofactors={summary['cofactors']}"
            self._log(job_id, f"🧪 Analysis: {analysis_msg}", "INFO")

            # Detailed cleaning report banner
            self._log(job_id, "===== PROTEIN CLEANING REPORT =====", "INFO")
            try:
                # try reading a JSON cleaning report (either analyzer or later processor)
                report_path = getattr(analyzer, 'clean_report', None)
                if report_path and Path(report_path).exists():
                    rep = json.loads(Path(report_path).read_text())
                    kept = rep.get('kept', {})
                    removed = rep.get('removed', {})
                    self._log(job_id, f"Kept  -> Waters={kept.get('waters',0)}, Ions={kept.get('ions',0)}, Solvents={kept.get('solvents',0)}, Cofactors={kept.get('cofactors',0)}", "INFO")
                    self._log(job_id, f"Removed-> Waters={removed.get('waters',0)}, Ions={removed.get('ions',0)}, Solvents={removed.get('solvents',0)}, Cofactors={removed.get('cofactors',0)}", "INFO")
                else:
                    self._log(job_id, f"Counts -> Waters={summary['waters']} Ions={summary['ions']} Solvents={summary['solvents']} Cofactors={summary['cofactors']}", "INFO")
            except Exception as e:
                self._log(job_id, f"(cleaning report not available) {e}", "WARNING")
            self._log(job_id, "===================================", "INFO")

            # Prepare protein (workspace)
            self._status(job_id, 45, "Preparing protein")
            prepared_protein_pdbqt = self.prepare_protein_in_workspace(raw_protein_path, prepared_dir, run_id, cleaning_policy, job_id)
            if not prepared_protein_pdbqt:
                raise Exception("Protein preparation failed")

            # Binding sites
            self._status(job_id, 55, "Identifying binding sites")
            pocket_info = identify_binding_site(raw_protein_path, prepared_protein_pdbqt, prepared_ligand_pdbqt, use_validated=True, return_n=5, detailed=True)
            if not pocket_info:
                raise Exception("Pocket detection failed")
            center = pocket_info["primary"]
            self._log(job_id, "===== BINDING SITES FOUND =====", "INFO")
            sites_list = []
            if "primary" in pocket_info:
                sites_list.append(("Primary", pocket_info["primary"]))
            if "modes" in pocket_info:
                for i, s in enumerate(pocket_info["modes"], 1):
                    sites_list.append((f"Alt {i}", s))
            for tag, s in sites_list:
                self._log(job_id,
                          f"{tag}: center=({s.get('center_x',0):.2f},{s.get('center_y',0):.2f},{s.get('center_z',0):.2f}) "
                          f"size=({s.get('size_x',0):.1f},{s.get('size_y',0):.1f},{s.get('size_z',0):.1f}) "
                          f"method={s.get('method','?')} confidence={s.get('confidence','?')}",
                          "INFO")
            self._log(job_id, "================================", "INFO")

            # Docking
            self._status(job_id, 65, "Running docking")
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
            self._log(job_id, f"🧷 Docked poses file: {docked_file}", "SUCCESS")

            # Parse results
            self._status(job_id, 85, "Parsing results")
            results_csv = docking_dir / "binding_scores.csv"
            df = parse_vina_results(docked_file, output_csv_path=results_csv, pocket_info=pocket_info)
            if df is None or df.empty:
                raise Exception("Results parsing failed")

            # Summary numbers
            best_affinity = float(df["Binding Affinity (kcal/mol)"].min())
            total_poses = int(len(df))
            self._log(job_id, f"📈 Best affinity: {best_affinity:.2f} kcal/mol; Poses: {total_poses}", "SUCCESS")

            # Print each pose
            self._log(job_id, "===== VINA POSES =====", "INFO")
            try:
                for _, row in df.iterrows():
                    pose_no = int(row.get('Pose', row.get('pose', 0)))
                    aff = float(row.get('Binding Affinity (kcal/mol)', row.get('Affinity', 0.0)))
                    rlb = row.get('RMSD_lb', row.get('rmsd_lb', ''))
                    rub = row.get('RMSD_ub', row.get('rmsd_ub', ''))
                    self._log(job_id, f"Pose {pose_no:02d} • Affinity={aff:.2f} kcal/mol • RMSD_lb={rlb} • RMSD_ub={rub}", "INFO")
            except Exception as e:
                self._log(job_id, f"Could not list poses: {e}", "WARNING")
            self._log(job_id, "======================", "INFO")

            # Collate all pockets for UI
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

            # Summary JSON for frontend
            self._status(job_id, 95, "Writing summary")
            summary_data = {
                "job_id": job_id,
                "status": "completed",
                "best_affinity": best_affinity,
                "total_poses": total_poses,
                "pocket_center": [center["center_x"], center["center_y"], center["center_z"]],
                "pocket_size": [center["size_x"], center["size_y"], center["size_z"]],
                "method": center.get("method", "unknown"),
                "confidence": center.get("confidence", "unknown"),
                "docked_file": docked_file,
                "results_csv": str(results_csv),
                "files": {
                    "raw_protein": raw_protein_path,
                    "prepared_protein": prepared_protein_pdbqt,
                    "prepared_ligand": prepared_ligand_pdbqt,
                    "docked_poses": docked_file,
                    "results_csv": str(results_csv)
                },
                "protein_analysis": analysis_msg,
                "cleaning_policy": cleaning_policy,
                "retention": retention,
                "run_directory": str(run_dir.resolve()),
                "all_pockets": all_pockets
            }
            (run_dir / "results.json").write_text(json.dumps(summary_data, indent=2, default=str))

            self._status(job_id, 100, "Pipeline completed")
            self._log(job_id, "🎉 Docking complete!", "SUCCESS")
            if self.status_callback:
                self.status_callback(job_id, "completed", 100, f"Docking completed! Best affinity: {best_affinity:.2f} kcal/mol")
            return summary_data
        except Exception as e:
            err = f"Pipeline failed: {e}"
            log_console("ERROR", err)
            print(traceback.format_exc())
            if self.status_callback:
                self.status_callback(job_id, "failed", 0, err)
            return {"job_id": job_id, "status": "failed", "error_message": str(e)}

# ========= Job status updates =========
def update_job_status(job_id: str, status: str, progress: int, message: str):
    if job_id in jobs:
        jobs[job_id].status = status
        jobs[job_id].progress = progress
        jobs[job_id].message = message
        if status in ["completed", "failed"]:
            jobs[job_id].completed_at = datetime.now().isoformat()

# ========= API endpoints =========
@app.post("/api/dock", response_model=dict)
async def start_docking(request: DockingRequest):
    job_id = str(uuid.uuid4())
    cleaning_policy = request.cleaning_policy or {
        "keep_waters": False, "keep_ions": True, "keep_solvents": False, "keep_cofactors": True
    }
    jobs[job_id] = JobStatus(job_id=job_id, status="queued", progress=0, message="Job queued", created_at=datetime.now().isoformat())
    job_logs[job_id] = []
    log_console("INFO", f"🧵 Spawn job: {job_id} (retention={request.retention})")

    def run_job():
        pipeline = DockingPipeline()
        pipeline.set_status_callback(update_job_status)
        pipeline.run_pipeline(job_id, request.protein_input, request.ligand_input, cleaning_policy, request.retention)

    threading.Thread(target=run_job, daemon=True).start()
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

    async def gen():
        sent = 0
        while True:
            logs = job_logs.get(job_id, [])
            for i in range(sent, len(logs)):
                yield f"data: {json.dumps(logs[i])}\n\n"
                sent = len(logs)
            job_status = jobs.get(job_id, {})
            if getattr(job_status, 'status', None) in ['completed', 'failed']:
                break
            await asyncio.sleep(0.5)

    return StreamingResponse(gen(), media_type="text/event-stream", headers={"Cache-Control":"no-cache","Connection":"keep-alive"})

@app.get("/api/jobs/{job_id}/result")
async def get_job_result(job_id: str):
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    job = jobs[job_id]
    if job.status != "completed":
        raise HTTPException(status_code=400, detail="Job not completed yet")
    try:
        search_paths = [Path("temp_runs").glob("**/results.json"), Path("data/results").glob("**/results.json")]
        for pattern in search_paths:
            for results_file in pattern:
                if results_file.exists():
                    with open(results_file) as f:
                        data = json.load(f)
                    if data.get("job_id") == job_id:
                        log_console("SUCCESS", f"📦 Results loaded for {job_id}")
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
        log_console("INFO", f"🗂️ Files listed: {len(file_list)}")
        return {"files": file_list}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error listing files: {e}")

@app.get("/api/download")
async def download_file(path: str = Query(..., description="File path to download")):
    file_path = Path(path)
    if not file_path.exists() or not file_path.is_file():
        log_console("ERROR", f"Download missing: {path}")
        raise HTTPException(status_code=404, detail=f"File not found: {path}")
    log_console("INFO", f"⬇️ Download: {file_path.name}")
    return FileResponse(str(file_path), filename=file_path.name, media_type='application/octet-stream')

@app.get("/api/jobs")
async def list_jobs():
    return list(jobs.values())

@app.post("/api/upload/protein")
async def upload_protein(file: UploadFile = File(...)):
    upload_dir = Path("uploads/proteins")
    upload_dir.mkdir(parents=True, exist_ok=True)
    file_path = upload_dir / file.filename
    with open(file_path, "wb") as f:
        f.write(await file.read())
    log_console("SUCCESS", f"📤 Protein uploaded -> {file_path}")
    return {"message": "Protein file uploaded", "file_path": str(file_path), "download_url": f"/api/download?path={file_path}"}

@app.post("/api/upload/ligand")
async def upload_ligand(file: UploadFile = File(...)):
    upload_dir = Path("uploads/ligands")
    upload_dir.mkdir(parents=True, exist_ok=True)
    file_path = upload_dir / file.filename
    with open(file_path, "wb") as f:
        f.write(await file.read())
    log_console("SUCCESS", f"📤 Ligand uploaded -> {file_path}")
    return {"message": "Ligand file uploaded", "file_path": str(file_path), "download_url": f"/api/download?path={file_path}"}

@app.get("/")
async def root():
    return {"message": "Dippy Docking API", "version": "1.0.0", "status": "running"}

if __name__ == "__main__":
    import uvicorn
    log_console("INFO", "🔊 Uvicorn starting at 0.0.0.0:8000 (reload=True)")
    uvicorn.run("api.docking_api:app", host="0.0.0.0", port=8000, reload=True)
