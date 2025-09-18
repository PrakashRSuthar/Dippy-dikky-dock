# backend/api/docking_api.py
import sys
from pathlib import Path
import traceback
import shutil
import concurrent.futures

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

jobs: Dict[str, JobStatus] = {}
job_logs: Dict[str, List[Dict[str, Any]]] = {}

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
                    dest_file = external_dest / ext_file.name
                    shutil.copy2(ext_file, dest_file)
            summary_data['permanent_location'] = str(permanent_path)
            with open(permanent_path / "summary.json", "w") as f:
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
    def set_status_callback(self, callback): self.status_callback = callback
    def set_log_callback(self, callback): self.log_callback = callback

    def _status(self, job_id: str, pct: int, msg: str):
        bar_len = 30
        filled = int(bar_len * pct // 100)
        bar = '█' * filled + '░' * (bar_len - filled)
        print(f"\r⏳ |{bar}| {pct}% - {msg}", end='', flush=True)
        if self.status_callback: self.status_callback(job_id, "running", pct, msg)
        self._log(job_id, f"[{pct}%] {msg}", "INFO")

    def _log(self, job_id: str, msg: str, level: str = "INFO"):
        ts = datetime.now().strftime("%H:%M:%S")
        entry = {"timestamp": ts, "level": level, "message": msg}
        job_logs.setdefault(job_id, []).append(entry)
        level_emoji = {"INFO":"ℹ️","SUCCESS":"✅","ERROR":"❌","WARNING":"⚠️"}.get(level,"ℹ️")
        print(f"\n{level_emoji} [{ts}] {msg}")
        if self.log_callback: self.log_callback(job_id, entry)

    def _ensure_file(self, job_id: str, path_str: Optional[str], label: str):
        if not path_str:
            raise FileNotFoundError(f"{label} missing: got None")
        p = Path(path_str)
        if not p.exists() or not p.is_file():
            raise FileNotFoundError(f"{label} missing: {path_str}")
        if p.stat().st_size < 64:
            raise IOError(f"{label} appears empty: {path_str}")
        self._log(job_id, f"📄 {label}: {p.name} ({p.stat().st_size/1024:.1f} KB)", "INFO")
        return str(p)

    def _with_timeout(self, job_id: str, seconds: int, func, *args, **kwargs):
        self._log(job_id, f"⏱️ Timeout set: {seconds}s for {func.__name__}", "INFO")
        with concurrent.futures.ThreadPoolExecutor(max_workers=1) as ex:
            fut = ex.submit(func, *args, **kwargs)
            try:
                return fut.result(timeout=seconds)
            except concurrent.futures.TimeoutError:
                self._log(job_id, f"❌ Step timed out: {func.__name__}", "ERROR")
                raise TimeoutError(f"Timeout in {func.__name__}")

    # Input resolvers
    def resolve_protein_input(self, user_input: str, raw_dir: Path, job_id: str) -> str:
        if not user_input or user_input.strip() == "": raise Exception("Empty protein input")
        user_input = user_input.strip()
        p = Path(user_input)
        if p.exists() and p.is_file():
            dst = raw_dir / p.name
            if dst.resolve() != p.resolve(): dst.write_bytes(p.read_bytes())
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
        if not user_input or user_input.strip() == "": raise Exception("Empty ligand input")
        user_input = user_input.strip()
        p = Path(user_input)
        if p.exists() and p.is_file():
            dst = raw_dir / p.name
            if dst.resolve() != p.resolve(): dst.write_bytes(p.read_bytes())
            ext = p.suffix.lower()
            self._log(job_id, f"📥 Ligand from local path ({ext}) -> {dst}", "INFO")
            if ext == ".pdbqt": return str(dst)
            elif ext in {".sdf", ".mol2", ".pdb"}:
                out = _convert_to_pdbqt(str(dst), str(prepared_dir))
                self._log(job_id, f"🔁 Converted ligand to PDBQT -> {out}", "INFO")
                return out
        self._log(job_id, f"💊 Fetching ligand from PubChem: {user_input}", "INFO")
        out = process_ligand(user_input, ligand_dir=str(raw_dir), prepared_dir=str(prepared_dir))
        if isinstance(out, dict):
            out = out.get('prepared_pdbqt') or out.get('path') or out.get('pdbqt')
        if not out or not isinstance(out, str):
            raise Exception("Ligand preparation did not return a file path")
        return out

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

    # Added: retention handler
    def _handle_retention(self, job_id: str, retention: str, workspace: Path, summary: dict):
        """
        Post-completion retention policy:
        - 'delete': remove the workspace directory immediately.
        - otherwise: keep files (log location).
        """
        try:
            if retention == "delete":
                self._log(job_id, "🗑️ Deleting workspace per retention=delete", "INFO")
                shutil.rmtree(workspace, ignore_errors=True)
                self._log(job_id, "✅ Workspace deleted", "SUCCESS")
            else:
                self._log(job_id, f"📦 Retention '{retention}': keeping workspace at {workspace}", "INFO")
        except Exception as e:
            self._log(job_id, f"⚠️ Retention handling failed: {e}", "WARNING")

    # Protein bounds and clamp
    def _protein_bounds_ca(self, pdb_path: str):
        xs=[]; ys=[]; zs=[]
        with open(pdb_path) as f:
            for line in f:
                if line.startswith("ATOM") and line[12:16].strip()=="CA":
                    parts = line[30:54].split()
                    if len(parts)>=3:
                        x,y,z = map(float, parts[:3])
                        xs.append(x); ys.append(y); zs.append(z)
        if not xs: return None
        return (min(xs),max(xs)),(min(ys),max(ys)),(min(zs),max(zs))

    def _clamp_box_to_protein(self, raw_protein_path: str, center: Dict) -> Dict:
        bounds = self._protein_bounds_ca(raw_protein_path)
        if not bounds: return center
        (xmin,xmax),(ymin,ymax),(zmin,zmax) = bounds
        cx,cy,cz = center["center_x"], center["center_y"], center["center_z"]
        sx,sy,sz = center["size_x"], center["size_y"], center["size_z"]
        def clamp_center(c, s, lo, hi):
            half = s/2.0
            if c+half < lo: c = lo + half + 1.5
            if c-half > hi: c = hi - half - 1.5
            return c
        cx = clamp_center(cx, sx, xmin, xmax)
        cy = clamp_center(cy, sy, ymin, ymax)
        cz = clamp_center(cz, sz, zmin, zmax)
        sx = max(14.0, min(26.0, sx)); sy = max(14.0, min(26.0, sy)); sz = max(14.0, min(26.0, sz))
        center.update({"center_x":cx,"center_y":cy,"center_z":cz,"size_x":sx,"size_y":sy,"size_z":sz})
        return center

    def run_pipeline(self, job_id: str, protein_input: str, ligand_input: str, cleaning_policy: dict, retention: str):
        try:
            self._log(job_id, f"🚀 Starting Docking Pipeline [{job_id[:8]}]", "INFO")
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
            raw_protein = self.resolve_protein_input(protein_input, raw_dir, job_id)
            raw_protein = self._ensure_file(job_id, raw_protein, "Raw protein")

            self._status(job_id, 25, "Resolving ligand")
            prepared_lig = self.resolve_ligand_input(ligand_input, raw_dir, prepared_dir, job_id)
            prepared_lig = self._ensure_file(job_id, prepared_lig, "Prepared ligand PDBQT")

            # Analyze and prepare protein
            self._status(job_id, 40, "Analyzing protein")
            analyzer = ProteinPreprocessor(raw_protein)
            analysis = analyzer.summarize_pdb(analyzer.pdb_path)
            self._log(job_id, f"🧪 Found: Waters={analysis['waters']} Ions={analysis['ions']} Solvents={analysis['solvents']} Cofactors={analysis['cofactors']}", "INFO")

            self._status(job_id, 45, "Preparing protein")
            prep_path = self._with_timeout(job_id, 180, self.prepare_protein_in_workspace, raw_protein, prepared_dir, run_id, cleaning_policy, job_id)
            try:
                prep_path = self._ensure_file(job_id, prep_path, "Prepared protein PDBQT")
            except Exception as e:
                self._log(job_id, f"⚠️ Prepared PDBQT missing ({e}); trying cleaned PDB", "WARNING")
                cleaned_pdb = prepared_dir / f"{Path(raw_protein).stem}_{run_id}_clean.pdb"
                prep_path = self._ensure_file(job_id, str(cleaned_pdb), "Cleaned protein PDB")
            self._status(job_id, 50, "Protein prepared")

            # Pocket identification
            self._status(job_id, 55, "Identifying binding sites")
            def find_pocket():
                return identify_binding_site(raw_protein, prep_path, prepared_lig, use_validated=True, return_n=5, detailed=True)
            pocket_info = self._with_timeout(job_id, 180, find_pocket)
            if not pocket_info or "primary" not in pocket_info:
                self._log(job_id, "⚠️ Pocket finder returned nothing; building default box", "WARNING")
                default_center = {"center_x": 0.0, "center_y": 0.0, "center_z": 0.0, "size_x": 22.0, "size_y": 22.0, "size_z": 22.0, "method": "default", "confidence": "low"}
                pocket_info = {"primary": default_center, "modes": []}
            center = pocket_info["primary"]
            self._log(job_id, f"🎯 Initial center=({center['center_x']:.1f},{center['center_y']:.1f},{center['center_z']:.1f}) size=({center['size_x']:.1f},{center['size_y']:.1f},{center['size_z']:.1f})", "INFO")

            # Clamp to protein bounds
            center = self._clamp_box_to_protein(raw_protein, center)
            self._log(job_id, f"🛠️ Box clamped: center=({center['center_x']:.1f},{center['center_y']:.1f},{center['center_z']:.1f}) size=({center['size_x']:.1f},{center['size_y']:.1f},{center['size_z']:.1f})", "INFO")
            self._status(job_id, 60, "Pocket ready")

            # Optional quick probe of alternate pocket
            modes = pocket_info.get("modes", [])
            if len(modes) >= 2:
                self._status(job_id, 63, "Probing alternate pocket")
                docking_dir = run_dir / "docking"
                docking_dir.mkdir(exist_ok=True)
                def quick_probe(c):
                    try:
                        eng = DockingEngine(prep_path, prepared_lig)
                        return eng.run_docking(center_x=c["center_x"], center_y=c["center_y"], center_z=c["center_z"],
                                               box_size_x=c["size_x"], box_size_y=c["size_y"], box_size_z=c["size_z"],
                                               output_dir=str(docking_dir), quick=True)
                    except Exception:
                        return None
                q1 = quick_probe(center)
                q2 = quick_probe(modes[1])
                def best_aff(path):
                    if not path or not Path(path).exists(): return None
                    vals=[]
                    try:
                        with open(path) as f:
                            for line in f:
                                if "REMARK VINA RESULT" in line:
                                    parts = line.split()
                                    vals.append(float(parts[3]))
                    except Exception:
                        return None
                    return min(vals) if vals else None
                a1 = best_aff(q1); a2 = best_aff(q2)
                if a2 is not None and (a1 is None or a2 < a1):
                    center = modes[1]
                    center = self._clamp_box_to_protein(raw_protein, center)
                    self._log(job_id, f"🔁 Switched to alternate pocket (quick better: {a2:.2f} < {a1 if a1 is not None else float('inf')})", "INFO")

            # Docking
            self._status(job_id, 65, "Running docking")
            docking_dir = run_dir / "docking"
            docking_dir.mkdir(exist_ok=True)
            engine = DockingEngine(prep_path, prepared_lig)
            def run_dock():
                return engine.run_docking(
                    center_x=center["center_x"], center_y=center["center_y"], center_z=center["center_z"],
                    box_size_x=center["size_x"], box_size_y=center["size_y"], box_size_z=center["size_z"],
                    output_dir=str(docking_dir)
                )
            docked_file = self._with_timeout(job_id, 600, run_dock)
            docked_file = self._ensure_file(job_id, docked_file, "Docked poses")
            self._status(job_id, 75, "Docking finished")

            # Results parsing
            self._status(job_id, 85, "Parsing results")
            results_csv = docking_dir / "binding_scores.csv"
            df = parse_vina_results(docked_file, output_csv_path=results_csv, pocket_info=pocket_info)
            if df is None or df.empty:
                self._log(job_id, "⚠️ Parser returned no rows; writing minimal CSV", "WARNING")
                with open(results_csv, "w") as f:
                    f.write("Pose,Binding Affinity (kcal/mol),RMSD_lb,RMSD_ub\n1,0.0,0.0,0.0\n")
                best_affinity, total_poses = 0.0, 1
            else:
                best_affinity = float(df["Binding Affinity (kcal/mol)"].min())
                total_poses = int(len(df))
                self._log(job_id, f"📈 Parsed {total_poses} poses; best={best_affinity:.2f} kcal/mol", "SUCCESS")
            self._status(job_id, 92, "Results ready")

            # Pockets summary for UI
            all_pockets = [{
                "center": [center.get("center_x",0), center.get("center_y",0), center.get("center_z",0)],
                "size": [center.get("size_x",22), center.get("size_y",22), center.get("size_z",22)],
                "confidence": center.get("confidence","unknown"),
                "method": center.get("method","unknown")
            }]
            for s in pocket_info.get("modes", []):
                all_pockets.append({
                    "center": [s.get("center_x",0), s.get("center_y",0), s.get("center_z",0)],
                    "size": [s.get("size_x",22), s.get("size_y",22), s.get("size_z",22)],
                    "confidence": s.get("confidence","unknown"),
                    "method": s.get("method","unknown")
                })

            summary = {
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
                    "raw_protein": raw_protein,
                    "prepared_protein": prep_path,
                    "prepared_ligand": prepared_lig,
                    "docked_poses": docked_file,
                    "results_csv": str(results_csv)
                },
                "protein_analysis": f"Waters={analysis['waters']} Ions={analysis['ions']} Solvents={analysis['solvents']} Cofactors={analysis['cofactors']}",
                "cleaning_policy": cleaning_policy,
                "retention": retention,
                "run_directory": str(run_dir.resolve()),
                "all_pockets": all_pockets
            }
            (run_dir / "results.json").write_text(json.dumps(summary, indent=2, default=str))

            self._status(job_id, 100, "Completed")
            self._log(job_id, "🎉 Docking complete!", "SUCCESS")
            if self.status_callback:
                self.status_callback(job_id, "completed", 100, f"Docking completed! Best affinity: {best_affinity:.2f} kcal/mol")
            self._handle_retention(job_id, retention, workspace, summary)
            return summary
        except Exception as e:
            msg = f"Pipeline failed: {e}"
            self._log(job_id, msg, "ERROR")
            print(traceback.format_exc())
            if self.status_callback: self.status_callback(job_id, "failed", 0, msg)
            return {"job_id": job_id, "status": "failed", "error_message": str(e)}

def update_job_status(job_id: str, status: str, progress: int, message: str):
    if job_id in jobs:
        jobs[job_id].status = status
        jobs[job_id].progress = progress
        jobs[job_id].message = message
        if status in ["completed", "failed"]:
            jobs[job_id].completed_at = datetime.now().isoformat()

@app.post("/api/dock", response_model=dict)
async def start_docking(request: DockingRequest):
    job_id = str(uuid.uuid4())
    cleaning_policy = request.cleaning_policy or {
        "keep_waters": False, "keep_ions": True, "keep_solvents": False, "keep_cofactors": True
    }
    jobs[job_id] = JobStatus(job_id=job_id, status="queued", progress=0, message="Job queued for processing", created_at=datetime.now().isoformat())
    job_logs[job_id] = []
    def run_job():
        pipeline = DockingPipeline()
        pipeline.set_status_callback(update_job_status)
        pipeline.set_log_callback(lambda j, e: None)
        try:
            pipeline.run_pipeline(job_id, request.protein_input, request.ligand_input, cleaning_policy, request.retention)
        except Exception as e:
            update_job_status(job_id, "failed", 0, f"Unhandled: {e}")
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
    async def log_generator():
        sent = 0
        while True:
            current = job_logs.get(job_id, [])
            for i in range(sent, len(current)):
                yield f"data: {json.dumps(current[i])}\n\n"
                sent = len(current)
            job_status = jobs.get(job_id, {})
            if getattr(job_status, 'status', None) in ['completed', 'failed']:
                break
            await asyncio.sleep(0.5)
    return StreamingResponse(log_generator(), media_type="text/event-stream", headers={"Cache-Control": "no-cache", "Connection": "keep-alive"})

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
    return {
        "message": "Protein file uploaded",
        "file_path": str(file_path),
        "download_url": f"/api/download?path={file_path}"
    }

@app.post("/api/upload/ligand")
async def upload_ligand(file: UploadFile = File(...)):
    upload_dir = Path("uploads/ligands")
    upload_dir.mkdir(parents=True, exist_ok=True)
    file_path = upload_dir / file.filename
    with open(file_path, "wb") as f:
        f.write(await file.read())
    return {
        "message": "Ligand file uploaded",
        "file_path": str(file_path),
        "download_url": f"/api/download?path={file_path}"
    }

@app.get("/")
async def root():
    return {"message": "Dippy Docking API", "version": "1.0.0", "status": "running"}

if __name__ == "__main__":
    import uvicorn
    uvicorn.run("api.docking_api:app", host="0.0.0.0", port=8000, reload=True)
