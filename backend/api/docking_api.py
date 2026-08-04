# docking_api.py
"""
Dippy Docking API — FastAPI backend
Paper Section III-G: REST endpoints + SSE log streaming + batch support
"""
import sys, os, json, uuid, threading, shutil, traceback, asyncio
import concurrent.futures
from pathlib import Path
from datetime import datetime
from contextlib import asynccontextmanager
from typing import Optional, List, Dict, Any

from fastapi import FastAPI, HTTPException, UploadFile, File, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, StreamingResponse
from pydantic import BaseModel

backend_dir = Path(__file__).parent.parent
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(backend_dir))
# Anchor every relative path (temp_runs/, data/, uploads/, benchmark_results/)
# to the repository root so the API behaves identically regardless of CWD.
os.chdir(PROJECT_ROOT)

try:
    from modules.protein_fetcher      import ProteinFetcher
    from modules.protein_preprocessor import prepare_protein, ProteinPreprocessor
    from modules.ligand_preparer       import process_ligand, _convert_to_pdbqt
    from modules.pocket_identifier     import identify_binding_site
    from modules.docking_engine        import DockingEngine
    from modules.result_parser         import parse_vina_results
    print(" All modules imported")
except ImportError as e:
    print(f" Module import error: {e}")
    sys.exit(1)

#  Shared job registries 
jobs:     Dict[str, Any] = {}
job_logs: Dict[str, List[Dict]] = {}

#  Pydantic models 
class DockingRequest(BaseModel):
    protein_input:   str
    ligand_input:    str
    job_name:        Optional[str]  = None
    cleaning_policy: Optional[dict] = {"keep_waters":False,"keep_ions":True,
                                        "keep_solvents":False,"keep_cofactors":True}
    retention:       Optional[str]  = "keep7d"

class BatchItemModel(BaseModel):
    protein_input: str
    ligand_input:  str
    job_name:      Optional[str] = None

class BatchDockingRequest(BaseModel):
    items:           List[BatchItemModel]
    cleaning_policy: Optional[dict] = {"keep_waters":False,"keep_ions":True,
                                        "keep_solvents":False,"keep_cofactors":True}
    retention:       Optional[str]  = "keep7d"

class JobStatus(BaseModel):
    job_id:       str
    status:       str
    progress:     int
    message:      str
    created_at:   str
    completed_at: Optional[str] = None
    result:       Optional[Dict[str, Any]] = None
    job_name:     Optional[str] = None

#  App lifecycle 
@asynccontextmanager
async def lifespan(app: FastAPI):
    print(" Dippy Docking API started"); yield
    print(" Dippy Docking API stopped")

app = FastAPI(title="Dippy Docking API", version="2.0.0", lifespan=lifespan)
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_credentials=True,
                   allow_methods=["*"], allow_headers=["*"])

#  Status helper (used by batch_processor too) 
def update_job_status(job_id, status, progress, message):
    if job_id in jobs:
        jobs[job_id].status    = status
        jobs[job_id].progress  = progress
        jobs[job_id].message   = message
        if status in ("completed","failed"):
            jobs[job_id].completed_at = datetime.now().isoformat()

# 
# Docking Pipeline
# 
class DockingPipeline:
    def __init__(self):
        self._status_cb = None; self._log_cb = None

    def set_status_callback(self, cb): self._status_cb = cb
    def set_log_callback(self, cb):    self._log_cb    = cb

    def _status(self, jid, pct, msg):
        bar = "█"*int(30*pct//100) + "░"*(30-int(30*pct//100))
        print(f"\r⏳ |{bar}| {pct}% {msg}", end="", flush=True)
        if self._status_cb: self._status_cb(jid,"running",pct,msg)
        self._log(jid,f"[{pct}%] {msg}")

    def _log(self, jid, msg, level="INFO"):
        ts    = datetime.now().strftime("%H:%M:%S")
        emoji = {"INFO":"ℹ️","SUCCESS":"✅","ERROR":"❌","WARNING":"⚠️"}.get(level,"ℹ️")
        entry = {"timestamp":ts,"level":level,"message":msg}
        job_logs.setdefault(jid,[]).append(entry)
        print(f"\n{emoji} [{ts}] {msg}")
        if self._log_cb: self._log_cb(jid,entry)

    def _ensure_file(self, jid, path_str, label):
        if not path_str: raise FileNotFoundError(f"{label}: None")
        p = Path(path_str)
        if not p.exists() or not p.is_file(): raise FileNotFoundError(f"{label}: {path_str}")
        if p.stat().st_size < 64: raise IOError(f"{label} too small: {path_str}")
        self._log(jid, f" {label}: {p.name} ({p.stat().st_size/1024:.1f} KB)")
        return str(p)

    def _timeout(self, jid, secs, fn, *args, **kwargs):
        self._log(jid, f"⏱ Timeout: {secs}s for {fn.__name__}")
        with concurrent.futures.ThreadPoolExecutor(max_workers=1) as ex:
            fut = ex.submit(fn, *args, **kwargs)
            try: return fut.result(timeout=secs)
            except concurrent.futures.TimeoutError:
                self._log(jid, f"Timeout in {fn.__name__} ({secs}s)", "ERROR")
                raise TimeoutError(f"Timeout: {fn.__name__}")

    def _compute_docking_timeout(self, center: dict, ligand_pdbqt: str = None) -> int:
        """Compute an adaptive docking timeout based on box volume and ligand complexity.

        Base timeout: 600s (10 min) for a reference 20×20×20 Å box.
        Scales with box volume^0.5 and adds time for flexible ligands.
        Hard-capped at 1800s (30 min) to prevent runaway jobs.
        """
        from modules.docking_engine import _estimate_rotatable_bonds, REFERENCE_BOX_VOL
        bx = float(center.get("size_x", 20))
        by = float(center.get("size_y", 20))
        bz = float(center.get("size_z", 20))
        box_vol = bx * by * bz
        vol_factor = (box_vol / REFERENCE_BOX_VOL) ** 0.5

        n_rot = 6
        if ligand_pdbqt:
            try:
                n_rot = _estimate_rotatable_bonds(Path(ligand_pdbqt))
            except Exception:
                pass
        rot_factor = 1.0 + 0.08 * max(0, n_rot - 6)

        base_timeout = 600
        timeout = int(base_timeout * vol_factor * rot_factor)
        return max(300, min(1800, timeout))

    #  Input resolvers 
    def _resolve_protein(self, inp, raw_dir, jid):
        inp = inp.strip()
        p   = Path(inp)
        if p.exists() and p.is_file():
            dst = raw_dir/p.name
            if dst.resolve()!=p.resolve(): dst.write_bytes(p.read_bytes())
            return str(dst)
        fetcher = ProteinFetcher(download_dir=str(raw_dir))
        if len(inp)==4 and inp.isalnum():
            self._log(jid, f" Fetching PDB: {inp}")
            return fetcher.fetch_from_pdb(inp)
        self._log(jid, f" Fetching AlphaFold: {inp}")
        return fetcher.fetch_from_alphafold(inp)

    def _resolve_ligand(self, inp, raw_dir, prep_dir, jid):
        inp = inp.strip(); p = Path(inp)
        if p.exists() and p.is_file():
            dst = raw_dir/p.name
            if dst.resolve()!=p.resolve(): dst.write_bytes(p.read_bytes())
            if p.suffix.lower()==".pdbqt": return str(dst)
            out = _convert_to_pdbqt(str(dst), str(prep_dir))
            return out
        self._log(jid, f" PubChem: {inp}")
        out = process_ligand(inp, ligand_dir=str(raw_dir), prepared_dir=str(prep_dir))
        if isinstance(out,dict):
            out = out.get("prepared_pdbqt") or out.get("path") or out.get("pdbqt")
        if not out or not isinstance(out,str):
            raise Exception("Ligand preparation returned no file path")
        return out

    def _prep_protein(self, pdb, prep_dir, ts, policy, jid):
        prep_dir.mkdir(parents=True, exist_ok=True)
        proc = ProteinPreprocessor(pdb)
        proc.output_dir  = prep_dir
        proc.base_stem   = f"{Path(pdb).stem}_{ts}"
        proc.tmp_cleaned_pdb = prep_dir/f"{proc.base_stem}_clean.pdb"
        proc.clean_report    = prep_dir/f"{proc.base_stem}_clean_report.json"
        proc.output_pdbqt    = prep_dir/f"{proc.base_stem}_prepared.pdbqt"
        self._log(jid, " Cleaning protein...")
        result = proc.process(interactive=False, default_policy=policy)
        self._log(jid, f" Protein prepared: {result}", "SUCCESS")
        return result

    def _handle_retention(self, jid, retention, workspace, summary):
        try:
            if retention=="delete":
                self._log(jid, " Deleting workspace (retention=delete)")
                shutil.rmtree(workspace, ignore_errors=True)
            elif retention in ("save","permanent","keep"):
                self._log(jid, f" Saving permanently to data/results (retention='{retention}')")
                self._save_permanent(jid, workspace, summary)
            else:
                self._log(jid, f" Retention '{retention}': keeping {workspace}")
        except Exception as e:
            self._log(jid, f"Retention handling failed: {e}", "WARNING")

    def _save_permanent(self, jid, workspace, summary):
        """Copy a finished run into data/results/<run_name> for permanent storage."""
        results_root = Path("data/results")
        results_root.mkdir(parents=True, exist_ok=True)
        run_name = workspace.name
        dest = results_root / run_name
        if dest.exists():
            shutil.rmtree(dest, ignore_errors=True)
        shutil.copytree(workspace, dest, dirs_exist_ok=True)
        results_file = dest / "run" / "results.json"
        if results_file.exists():
            try:
                data = json.loads(results_file.read_text(encoding="utf-8"))
                data["permanent_location"] = str(dest)
                data["retention"] = "save"
                results_file.write_text(json.dumps(data, indent=2, default=str), encoding="utf-8")
            except Exception as e:
                self._log(jid, f"Could not update permanent results.json: {e}", "WARNING")
        self._log(jid, f" Permanent run saved → {dest}")

    #  Main pipeline 
    def run_pipeline(self, jid, protein_input, ligand_input, cleaning_policy, retention):
        try:
            self._log(jid, f" Pipeline [{jid[:8]}]")
            ts        = datetime.now().strftime("%Y%m%d_%H%M%S")
            workspace = Path(f"temp_runs/{ts}_{jid[:8]}")
            run_dir   = workspace/"run"
            raw_dir   = run_dir/"raw"
            prep_dir  = run_dir/"prepared"
            for d in [workspace,run_dir,raw_dir,prep_dir]: d.mkdir(parents=True, exist_ok=True)

            # 1) Resolve inputs
            self._status(jid,15,"Resolving protein")
            raw_prot = self._resolve_protein(protein_input, raw_dir, jid)
            raw_prot = self._ensure_file(jid, raw_prot, "Raw protein")

            self._status(jid,25,"Resolving ligand")
            prep_lig = self._resolve_ligand(ligand_input, raw_dir, prep_dir, jid)
            prep_lig = self._ensure_file(jid, prep_lig, "Ligand PDBQT")

            # 2) Analyse & prepare protein
            self._status(jid,40,"Analysing protein")
            analyser = ProteinPreprocessor(raw_prot)
            analysis = analyser.summarize_pdb(analyser.pdb_path)
            self._log(jid, f" Waters={analysis['waters']} Ions={analysis['ions']} "
                           f"Solvents={analysis['solvents']} Cofactors={analysis['cofactors']}")

            self._status(jid,45,"Preparing protein")
            prep_path = self._timeout(jid,180,self._prep_protein,
                                      raw_prot,prep_dir,ts,cleaning_policy,jid)
            try:
                prep_path = self._ensure_file(jid, prep_path, "Prepared PDBQT")
            except Exception:
                cleaned = prep_dir/f"{Path(raw_prot).stem}_{ts}_clean.pdb"
                prep_path = self._ensure_file(jid, str(cleaned), "Cleaned PDB")

            # 3) Pocket identification (3-method consensus)
            self._status(jid,55,"Identifying binding sites")
            pocket_info = self._timeout(jid,180,identify_binding_site,
                                        raw_prot,prep_path,prep_lig,True,5,True)
            if not pocket_info or "primary" not in pocket_info:
                self._log(jid," Pocket finder returned nothing; default box","WARNING")
                pocket_info = {"primary":{"center_x":0,"center_y":0,"center_z":0,
                                          "size_x":22,"size_y":22,"size_z":22,
                                          "method":"default","confidence":"low"},
                               "modes":[],"drift_events":0}
            center = pocket_info["primary"]
            self._log(jid, f" Center=({center['center_x']:.1f},{center['center_y']:.1f},"
                           f"{center['center_z']:.1f}) "
                           f"Box=({center['size_x']:.1f},{center['size_y']:.1f},{center['size_z']:.1f}) "
                           f"method={center.get('method','?')} drift_events={pocket_info.get('drift_events',0)}")
            self._status(jid,60,"Pocket identified")

            # 4) Quick-probe alternate pocket (paper Section III-F)
            modes = pocket_info.get("modes",[])
            if len(modes)>=2:
                self._status(jid,63,"Probing alternate pocket")
                dock_dir = run_dir/"docking"; dock_dir.mkdir(exist_ok=True)
                def _quick(c):
                    try:
                        eng = DockingEngine(prep_path,prep_lig)
                        return eng.run_docking(center_x=c["center_x"],center_y=c["center_y"],
                                               center_z=c["center_z"],box_size_x=c["size_x"],
                                               box_size_y=c["size_y"],box_size_z=c["size_z"],
                                               output_dir=str(dock_dir),quick=True)
                    except: return None
                def _best_aff(path):
                    if not path or not Path(path).exists(): return None
                    vals=[]
                    with open(path) as f:
                        for line in f:
                            if "REMARK VINA RESULT" in line:
                                try: vals.append(float(line.split()[3]))
                                except: pass
                    return min(vals) if vals else None
                a1=_best_aff(_quick(center)); a2=_best_aff(_quick(modes[1]))
                if a2 is not None and (a1 is None or a2<a1):
                    center=modes[1]
                    self._log(jid,f" Alternate pocket selected (Δ={a2:.2f} vs {a1})")

            # 5) Full docking run
            self._status(jid,65,"Running docking")
            dock_dir = run_dir/"docking"; dock_dir.mkdir(exist_ok=True)
            engine   = DockingEngine(prep_path,prep_lig)
            dock_timeout = self._compute_docking_timeout(center, prep_lig)
            self._log(jid, f"Docking timeout: {dock_timeout}s (box=({center['size_x']:.1f},{center['size_y']:.1f},{center['size_z']:.1f}) Å)")
            docked   = self._timeout(jid,dock_timeout,engine.run_docking,
                                     center_x=center["center_x"],center_y=center["center_y"],
                                     center_z=center["center_z"],box_size_x=center["size_x"],
                                     box_size_y=center["size_y"],box_size_z=center["size_z"],
                                     output_dir=str(dock_dir))
            docked = self._ensure_file(jid,docked,"Docked poses")
            self._status(jid,75,"Docking complete")

            # 6) Parse results — REAL RMSD values, no random data
            self._status(jid,85,"Parsing results")
            csv_path = dock_dir/"binding_scores.csv"
            df = parse_vina_results(docked, output_csv_path=csv_path, pocket_info=pocket_info)
            if df is None or df.empty:
                self._log(jid," No poses parsed; minimal CSV","WARNING")
                csv_path.write_text("Pose,Binding Affinity (kcal/mol),RMSD_lb (A),RMSD_ub (A)\n1,0.0,0.0,0.0\n")
                best_aff, total_poses = 0.0, 1
            else:
                best_aff   = float(df["Binding Affinity (kcal/mol)"].min())
                total_poses= int(len(df))
                self._log(jid, f" {total_poses} poses; best={best_aff:.2f} kcal/mol; "
                               f"RMSD_lb range={df['RMSD_lb (A)'].min():.2f}—{df['RMSD_ub (A)'].max():.2f} Å","SUCCESS")
            self._status(jid,92,"Results ready")

            # 7) Build summary
            all_pockets=[{"center":[center.get("center_x",0),center.get("center_y",0),center.get("center_z",0)],
                           "size":[center.get("size_x",22),center.get("size_y",22),center.get("size_z",22)],
                           "confidence":center.get("confidence","?"),"method":center.get("method","?")}]
            for s in pocket_info.get("modes",[]):
                all_pockets.append({"center":[s.get("center_x",0),s.get("center_y",0),s.get("center_z",0)],
                                    "size":[s.get("size_x",22),s.get("size_y",22),s.get("size_z",22)],
                                    "confidence":s.get("confidence","?"),"method":s.get("method","?")})

            # Pose data for frontend — REAL values from parsed DataFrame
            pose_rows = []
            if df is not None and not df.empty:
                for _,row in df.iterrows():
                    pose_rows.append({
                        "pose":   int(row["Pose"]),
                        "affinity": float(row["Binding Affinity (kcal/mol)"]),
                        "rmsd_lb":  float(row.get("RMSD_lb (A)",0.0)),
                        "rmsd_ub":  float(row.get("RMSD_ub (A)",0.0)),
                    })

            summary = {
                "job_id":jid,"status":"completed",
                "best_affinity":best_aff,"total_poses":total_poses,
                "pocket_center":[center["center_x"],center["center_y"],center["center_z"]],
                "pocket_size":[center["size_x"],center["size_y"],center["size_z"]],
                "method":center.get("method","?"),"confidence":center.get("confidence","?"),
                "drift_events":pocket_info.get("drift_events",0),
                "consensus_clusters":pocket_info.get("consensus_clusters",1),
                "methods_used":pocket_info.get("methods_used",[]),
                "docked_file":docked,"results_csv":str(csv_path),
                "pose_data":pose_rows,   # ← real RMSD data for frontend
                "files":{"raw_protein":raw_prot,"prepared_protein":prep_path,
                         "prepared_ligand":prep_lig,"docked_poses":docked,"results_csv":str(csv_path)},
                "protein_analysis":f"Waters={analysis['waters']} Ions={analysis['ions']} "
                                   f"Solvents={analysis['solvents']} Cofactors={analysis['cofactors']}",
                "cleaning_policy":cleaning_policy,"retention":retention,
                "run_directory":str(run_dir.resolve()),
                "all_pockets":all_pockets,
            }
            (run_dir/"results.json").write_text(json.dumps(summary,indent=2,default=str))

            self._status(jid,100,"Completed")
            self._log(jid,f" Done! Best affinity={best_aff:.2f} kcal/mol","SUCCESS")
            if self._status_cb:
                self._status_cb(jid,"completed",100,f"Done! Best affinity: {best_aff:.2f} kcal/mol")
            self._handle_retention(jid,retention,workspace,summary)
            return summary

        except Exception as e:
            msg = f"Pipeline failed: {e}"
            self._log(jid,msg,"ERROR"); print(traceback.format_exc())
            if self._status_cb: self._status_cb(jid,"failed",0,msg)
            return {"job_id":jid,"status":"failed","error_message":str(e)}

#  Routes 
@app.post("/api/dock", response_model=dict)
async def start_docking(req: DockingRequest):
    jid     = str(uuid.uuid4())
    policy  = req.cleaning_policy or {"keep_waters":False,"keep_ions":True,
                                       "keep_solvents":False,"keep_cofactors":True}
    jobs[jid] = JobStatus(job_id=jid,status="queued",progress=0,
                          message="Queued",created_at=datetime.now().isoformat(),
                          job_name=req.job_name or f"job_{jid[:8]}")
    job_logs[jid] = []
    def _run():
        p=DockingPipeline(); p.set_status_callback(update_job_status)
        try:
            result = p.run_pipeline(jid,req.protein_input,req.ligand_input,policy,req.retention)
            if result and result.get("status")=="completed":
                jobs[jid].result = result
                _auto_save_job(jid,"completed",
                               f"Best affinity: {result.get('best_affinity')} kcal/mol")
        except Exception as e:
            update_job_status(jid,"failed",0,str(e))
    threading.Thread(target=_run, daemon=True).start()
    return {"job_id":jid,"status":"queued","message":"Docking started"}

@app.post("/api/dock/batch", response_model=dict)
async def start_batch(req: BatchDockingRequest):
    if not req.items: raise HTTPException(400,"No items")
    policy   = req.cleaning_policy or {"keep_waters":False,"keep_ions":True,
                                        "keep_solvents":False,"keep_cofactors":True}
    batch_id = str(uuid.uuid4())
    results  = []
    for it in req.items:
        jid = str(uuid.uuid4())
        jobs[jid]=JobStatus(job_id=jid,status="queued",progress=0,
                            message=f"Batch {batch_id[:8]}: queued",
                            created_at=datetime.now().isoformat(),
                            job_name=it.job_name or f"job_{jid[:8]}")
        job_logs[jid]=[]
        def _run(j=jid, pi=it.protein_input, li=it.ligand_input):
            p=DockingPipeline(); p.set_status_callback(update_job_status)
            try:
                result = p.run_pipeline(j,pi,li,policy,req.retention or "keep7d")
                if result and result.get("status")=="completed":
                    jobs[j].result = result
                    _auto_save_job(j,"completed",
                                   f"Best affinity: {result.get('best_affinity')} kcal/mol")
            except Exception as e: update_job_status(j,"failed",0,str(e))
        threading.Thread(target=_run,daemon=True).start()
        results.append({"job_id":jid,"protein_input":it.protein_input,
                        "ligand_input":it.ligand_input,"job_name":it.job_name or f"job_{jid[:8]}"})
    return {"batch_id":batch_id,"count":len(results),"jobs":results}

@app.get("/api/jobs/{job_id}/status", response_model=JobStatus)
async def get_status(job_id: str):
    if job_id in jobs:
        return jobs[job_id]
    row = _get_job_row(job_id)
    if row:
        if row["status"] in ("completed","failed"):
            progress = 100 if row["status"]=="completed" else 0
        else:
            progress = 50
        return JobStatus(
            job_id=row["id"], status=row["status"], progress=progress,
            message=row["error_message"] or f"Job {row['status']}",
            created_at=row["created_at"] or "",
            completed_at=row["completed_at"],
            job_name=row["name"],
            result=None,
        )
    raise HTTPException(404,"Job not found")

@app.get("/api/jobs/{job_id}/logs")
async def stream_logs(job_id: str):
    if job_id not in jobs and not _get_job_row(job_id):
        raise HTTPException(404,"Job not found")
    async def _gen():
        sent=0
        last_state=None
        while True:
            cur=job_logs.get(job_id,[])
            for i in range(sent,len(cur)):
                yield f"data: {json.dumps(cur[i])}\n\n"
            sent=len(cur)
            st=None
            if job_id in jobs:
                st=jobs[job_id].status
            else:
                row=_get_job_row(job_id)
                st=row["status"] if row else None
            if st is None:
                break
            if st!=last_state:
                yield f"data: {json.dumps({'timestamp':datetime.now().strftime('%H:%M:%S'),'level':'INFO','message':f'Job {st}'})}\n\n"
                last_state=st
            if st in ("completed","failed"):
                break
            await asyncio.sleep(0.5)
    return StreamingResponse(_gen(), media_type="text/event-stream",
                             headers={"Cache-Control":"no-cache","Connection":"keep-alive"})

@app.get("/api/jobs/{job_id}/result")
async def get_result(job_id: str):
    if job_id in jobs:
        st=jobs[job_id].status
        if st=="completed" and jobs[job_id].result:
            return jobs[job_id].result
        if st=="failed":
            raise HTTPException(400, jobs[job_id].message or "Job failed")
        if st in ("running","queued"):
            raise HTTPException(400, f"Job still {st}")
    row=_get_job_row(job_id)
    if row:
        if row["status"]=="failed":
            raise HTTPException(400, row["error_message"] or "Job failed")
        if row["result_data"]:
            try:
                return json.loads(row["result_data"])
            except Exception:
                pass
    rf,data=_find_result_file(job_id)
    if rf is not None:
        return _normalize_result(data, job_id)
    raise HTTPException(404,"Results file not found")

@app.get("/api/jobs/{job_id}/files")
async def list_files(job_id: str):
    result = await get_result(job_id)
    file_list=[]
    for name,path in result.get("files",{}).items():
        fp=Path(path)
        if fp.exists():
            file_list.append({"name":name.replace("_"," ").title(),"filename":fp.name,
                              "path":str(path),"size":f"{fp.stat().st_size/1024:.1f} KB",
                              "download_url":f"/api/download?path={path}"})
    return {"files":file_list}

@app.get("/api/download")
async def download(path: str = Query(...)):
    fp=Path(path)
    if not fp.exists() or not fp.is_file(): raise HTTPException(404,f"Not found: {path}")
    return FileResponse(str(fp),filename=fp.name,media_type="application/octet-stream")

@app.get("/api/jobs")
async def list_jobs():
    merged = {jid: _job_to_dict(j) for jid,j in jobs.items()}
    conn = _get_db()
    for row in conn.execute("SELECT * FROM jobs").fetchall():
        if row["id"] not in merged:
            merged[row["id"]] = dict(row)
    conn.close()
    return sorted(merged.values(), key=lambda j: (j.get("created_at") or ""), reverse=True)

@app.post("/api/upload/protein")
async def upload_protein(file: UploadFile = File(...)):
    d=Path("uploads/proteins"); d.mkdir(parents=True, exist_ok=True)
    fp=d/file.filename; fp.write_bytes(await file.read())
    return {"message":"Uploaded","file_path":str(fp),"download_url":f"/api/download?path={fp}"}

@app.post("/api/upload/ligand")
async def upload_ligand(file: UploadFile = File(...)):
    d=Path("uploads/ligands"); d.mkdir(parents=True, exist_ok=True)
    fp=d/file.filename; fp.write_bytes(await file.read())
    return {"message":"Uploaded","file_path":str(fp),"download_url":f"/api/download?path={fp}"}

@app.get("/")
async def root(): return {"message":"Dippy Docking API","version":"2.0.0","status":"running"}

if __name__=="__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000, reload=False)


# ══════════════════════════════════════════════════════════════════════════════
#  SQLite Persistence & Desktop App Endpoints
# ══════════════════════════════════════════════════════════════════════════════

import sqlite3 as _sqlite3

_DB_PATH = Path("data/dippydock.db")


# ── Database helpers ────────────────────────────────────────────────────────

def _get_db():
    _DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    conn = _sqlite3.connect(str(_DB_PATH), timeout=10)
    conn.row_factory = _sqlite3.Row
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA busy_timeout=5000")
    return conn


def _init_db():
    conn = _get_db()
    conn.executescript("""
        CREATE TABLE IF NOT EXISTS jobs (
            id              TEXT PRIMARY KEY,
            name            TEXT DEFAULT '',
            status          TEXT DEFAULT 'queued',
            protein_input   TEXT,
            ligand_inputs   TEXT,
            created_at      TEXT,
            started_at      TEXT,
            completed_at    TEXT,
            duration_s      REAL,
            best_affinity   REAL,
            best_rmsd       REAL,
            error_message   TEXT,
            run_directory   TEXT,
            result_data     TEXT
        );
        CREATE TABLE IF NOT EXISTS settings (
            key   TEXT PRIMARY KEY,
            value TEXT
        );
        CREATE INDEX IF NOT EXISTS idx_jobs_status   ON jobs(status);
        CREATE INDEX IF NOT EXISTS idx_jobs_created  ON jobs(created_at DESC);
    """)
    conn.commit()
    conn.close()


def _human_size(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} TB"


def _format_duration(secs):
    if secs is None:
        return None
    secs = float(secs)
    if secs < 60:
        return f"{secs:.0f}s"
    if secs < 3600:
        return f"{int(secs//60)}m {int(secs%60)}s"
    return f"{int(secs//3600)}h {int((secs%3600)//60)}m"


def _job_to_dict(job):
    if isinstance(job, JobStatus):
        return job.model_dump()
    return dict(job)


def _get_job_row(job_id):
    conn = _get_db()
    row = conn.execute("SELECT * FROM jobs WHERE id = ?", (job_id,)).fetchone()
    conn.close()
    return row


def _find_result_file(job_id: str):
    """Locate a results.json / summary.json file for a job, on disk."""
    for root in [Path("temp_runs"), Path("data/results")]:
        if not root.exists():
            continue
        for pattern in ("**/results.json", "**/summary.json",
                        "results.json", "summary.json"):
            for rf in root.glob(pattern):
                try:
                    d = json.loads(rf.read_text(encoding="utf-8"))
                    if d.get("job_id") == job_id:
                        return rf, d
                    # Legacy summary.json runs carry pipeline_info.id (short id).
                    pid = (d.get("pipeline_info") or {}).get("id")
                    if pid and pid == job_id[:8]:
                        return rf, d
                except Exception:
                    continue
    return None, None


def _normalize_result(data: dict, job_id: str) -> dict:
    """Convert an old summary.json-style result into the frontend's expected shape."""
    if isinstance(data.get("status"), str) and data.get("status") in ("completed", "failed") \
       and "best_affinity" in data and "files" in data:
        return data
    res    = data.get("results", {}) or {}
    inp    = data.get("inputs", {}) or {}
    pip    = data.get("pipeline_info", {}) or {}
    files  = data.get("files", {}) or {}
    binding = data.get("binding_sites", {}) or {}
    primary = binding.get("primary", {}) or {}
    center = [primary.get("center_x", 0), primary.get("center_y", 0), primary.get("center_z", 0)]
    size   = [primary.get("size_x", 22), primary.get("size_y", 22), primary.get("size_z", 22)]
    all_pockets = []
    if primary:
        all_pockets.append({"center": center, "size": size,
                            "confidence": primary.get("confidence", "?"),
                            "method": primary.get("method", "?")})
    for s in binding.get("modes", []) or []:
        all_pockets.append({"center": [s.get("center_x", 0), s.get("center_y", 0), s.get("center_z", 0)],
                            "size": [s.get("size_x", 22), s.get("size_y", 22), s.get("size_z", 22)],
                            "confidence": s.get("confidence", "?"),
                            "method": s.get("method", "?")})
    docked = res.get("poses_file") or files.get("docked_poses", "")
    csv    = res.get("results_csv") or files.get("results_csv", "")
    run_dir = ""
    ws = pip.get("workspace")
    if ws:
        try:
            run_dir = str(Path(ws).resolve())
        except Exception:
            run_dir = str(ws)
    return {
        "job_id": job_id, "status": "completed",
        "best_affinity": res.get("best_affinity", 0.0),
        "total_poses": res.get("total_poses", 0),
        "pocket_center": center, "pocket_size": size,
        "method": primary.get("method", "?"), "confidence": primary.get("confidence", "?"),
        "drift_events": binding.get("drift_events", 0),
        "consensus_clusters": binding.get("consensus_clusters", 1),
        "methods_used": binding.get("methods_used", []),
        "docked_file": docked, "results_csv": csv,
        "pose_data": [],
        "files": files,
        "protein_analysis": f"Protein={inp.get('protein','?')} Ligand={inp.get('ligand','?')}",
        "cleaning_policy": {}, "retention": "keep7d",
        "run_directory": run_dir,
        "all_pockets": all_pockets,
    }


def _deep_merge(base, override):
    out = dict(base)
    for k, v in (override or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = v
    return out


_DEFAULT_SETTINGS = {
    "docking":   {"exhaustiveness": 8, "num_modes": 9, "energy_range": 3},
    "box":       {"auto_box": True, "min_axis": 20, "max_axis": 30},
    "cleaning":  {"keep_waters": False, "keep_ions": True, "keep_solvents": False, "keep_cofactors": True},
    "retention": {"retention": "permanent"},
}


def _scan_run_dirs(root: Path):
    """Scan a directory for run folders with frontend-friendly metadata."""
    if not root.exists() or not root.is_dir():
        return []
    runs = []
    for entry in sorted(root.iterdir(),
                        key=lambda e: e.stat().st_mtime, reverse=True):
        if not entry.is_dir():
            continue
        has_result = any(p.name in ("results.json", "summary.json")
                         for p in entry.rglob("*"))
        file_count = sum(1 for _ in entry.rglob("*") if _.is_file())
        total_size = sum(f.stat().st_size for f in entry.rglob("*") if f.is_file())
        runs.append({
            "name": entry.name,
            "date": datetime.fromtimestamp(entry.stat().st_mtime).isoformat(),
            "created_at": datetime.fromtimestamp(
                entry.stat().st_ctime).isoformat(),
            "file_count": file_count,
            "total_size": total_size,
            "total_size_bytes": total_size,
            "total_size_human": _human_size(total_size),
            "has_result": has_result,
        })
    return runs


def _backfill_from_results(conn, job_id):
    """Find a results.json matching *job_id* and back-fill run data."""
    for root in [Path("temp_runs"), Path("data/results")]:
        if not root.exists():
            continue
        for rf in root.glob("**/results.json"):
            try:
                data = json.loads(rf.read_text())
                if data.get("job_id") != job_id:
                    continue
                run_dir = str(rf.parent.resolve())
                result_json = json.dumps(data, default=str)
                conn.execute("""
                    UPDATE jobs SET
                        run_directory = COALESCE(run_directory, ?),
                        best_affinity = COALESCE(best_affinity, ?),
                        result_data   = ?
                    WHERE id = ?
                """, (run_dir, data.get("best_affinity"), result_json, job_id))
                return
            except Exception:
                pass
    # Fall back to old-style summary.json files
    for root in [Path("temp_runs"), Path("data/results")]:
        if not root.exists():
            continue
        for rf in root.glob("**/summary.json"):
            try:
                data = json.loads(rf.read_text())
                if data.get("pipeline_info", {}).get("id") != job_id:
                    continue
                run_dir = str(rf.parent.resolve())
                result_json = json.dumps(_normalize_result(data, job_id), default=str)
                conn.execute("""
                    UPDATE jobs SET
                        run_directory = COALESCE(run_directory, ?),
                        best_affinity = COALESCE(best_affinity, ?),
                        result_data   = ?
                    WHERE id = ?
                """, (run_dir, data.get("results", {}).get("best_affinity"), result_json, job_id))
                return
            except Exception:
                pass


def _migrate_existing_runs():
    """On startup, scan temp_runs/ + data/results/ and import results.json not already in SQLite."""
    conn = _get_db()
    imported = 0
    seen = set()
    for results_file in list(Path("temp_runs").glob("**/results.json")) if Path("temp_runs").exists() else []:
        try:
            data = json.loads(results_file.read_text())
            jid = data.get("job_id")
            if not jid:
                continue
            if jid in seen:
                continue
            seen.add(jid)
            existing = conn.execute("SELECT id FROM jobs WHERE id=?", (jid,)).fetchone()
            if existing:
                _backfill_from_results(conn, jid)
                continue
            run_dir = str(results_file.parent.resolve())
            result_json = json.dumps(data, default=str)
            conn.execute("""
                INSERT OR IGNORE INTO jobs
                (id, name, status, created_at, completed_at,
                 best_affinity, run_directory, result_data)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                jid,
                f"job_{jid[:8]}",
                data.get("status", "completed"),
                data.get("created_at"),
                data.get("completed_at"),
                data.get("best_affinity"),
                run_dir,
                result_json,
            ))
            imported += 1
        except Exception as exc:
            print(f"Migration error for {results_file}: {exc}")

    # Old-style summary.json runs (legacy CleanDockingPipeline output)
    for root in ([Path("temp_runs"), Path("data/results")]
                 if Path("data/results").exists() else [Path("temp_runs")]):
        if not root.exists():
            continue
        for results_file in root.glob("**/summary.json"):
            try:
                data = json.loads(results_file.read_text())
                jid = data.get("pipeline_info", {}).get("id")
                if not jid or jid in seen:
                    continue
                seen.add(jid)
                existing = conn.execute("SELECT id FROM jobs WHERE id=?", (jid,)).fetchone()
                if existing:
                    _backfill_from_results(conn, jid)
                    continue
                run_dir = str(results_file.parent.resolve())
                result_json = json.dumps(_normalize_result(data, jid), default=str)
                conn.execute("""
                    INSERT OR IGNORE INTO jobs
                    (id, name, status, created_at, completed_at,
                     best_affinity, run_directory, result_data)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """, (
                    jid,
                    f"job_{jid[:8]}",
                    data.get("status", "completed"),
                    data.get("pipeline_info", {}).get("timestamp"),
                    None,
                    data.get("results", {}).get("best_affinity"),
                    run_dir,
                    result_json,
                ))
                imported += 1
            except Exception as exc:
                print(f"Migration error for {results_file}: {exc}")
    conn.commit()
    conn.close()
    if imported:
        print(f"  Migrated {imported} existing run(s) to SQLite")


def _auto_save_job(job_id, status, message=""):
    """Insert or update a job record on every status change."""
    try:
        conn = _get_db()
        now = datetime.now().isoformat()
        existing = conn.execute(
            "SELECT * FROM jobs WHERE id = ?", (job_id,)
        ).fetchone()

        if not existing:
            job_obj = jobs.get(job_id)
            created_at = now
            if job_obj and getattr(job_obj, "created_at", None):
                created_at = job_obj.created_at
            name = getattr(job_obj, "job_name", None) or f"job_{job_id[:8]}"
            started_at = now if status == "running" else None
            conn.execute("""
                INSERT INTO jobs (id, name, status, created_at, started_at)
                VALUES (?, ?, ?, ?, ?)
            """, (job_id, name, status, created_at, started_at))
        else:
            sets, params = ["status = ?"], [status]
            if status == "running" and not existing["started_at"]:
                sets.append("started_at = ?")
                params.append(now)
            if status in ("completed", "failed"):
                sets.append("completed_at = ?")
                params.append(now)
                if existing["started_at"]:
                    try:
                        start = datetime.fromisoformat(existing["started_at"])
                        sets.append("duration_s = ?")
                        params.append((datetime.now() - start).total_seconds())
                    except Exception:
                        pass
            if status == "failed" and message:
                sets.append("error_message = ?")
                params.append(message)
            params.append(job_id)
            conn.execute(
                f"UPDATE jobs SET {', '.join(sets)} WHERE id = ?", params
            )

        if status == "completed":
            _backfill_from_results(conn, job_id)

        conn.commit()
        conn.close()
    except Exception as exc:
        print(f"  Auto-save error for {job_id}: {exc}")


# ── Wrap update_job_status for auto-save ────────────────────────────────────

_original_update_job_status = update_job_status

def _persisting_update_job_status(job_id, status, progress, message):
    _original_update_job_status(job_id, status, progress, message)
    _auto_save_job(job_id, status, message)

update_job_status = _persisting_update_job_status

# ── Run DB init + migration at import time ──────────────────────────────────
_init_db()
_migrate_existing_runs()


# ── Pydantic models for new endpoints ───────────────────────────────────────

class JobRename(BaseModel):
    name: str

class SettingsPayload(BaseModel):
    settings: Dict[str, Any]


# ── New API endpoints ───────────────────────────────────────────────────────

@app.get("/api/history")
async def list_history(
    status: Optional[str] = Query(None),
    limit:  int = Query(50, ge=1, le=500),
    offset: int = Query(0, ge=0),
):
    """List all persisted jobs, newest first."""
    conn = _get_db()
    where = ""
    params: list = []
    if status:
        where = " WHERE status = ?"
        params.append(status)
    rows = conn.execute(
        f"SELECT * FROM jobs{where} ORDER BY created_at DESC LIMIT ? OFFSET ?",
        params + [limit, offset],
    ).fetchall()
    total = conn.execute(
        f"SELECT COUNT(*) as cnt FROM jobs{where}", params
    ).fetchone()["cnt"]
    conn.close()
    jobs_out = []
    for r in rows:
        j = dict(r)
        j["duration"] = _format_duration(j.get("duration_s"))
        j["rmsd"] = j.get("best_rmsd")
        rd = {}
        if j.get("result_data"):
            try:
                rd = json.loads(j["result_data"])
            except Exception:
                rd = {}
        j["protein_name"] = j.get("protein_input") or rd.get("protein_name") \
            or (rd.get("inputs", {}) or {}).get("protein")
        j["ligand_name"] = rd.get("ligand_name") or (rd.get("inputs", {}) or {}).get("ligand")
        j["updated_at"] = j.get("completed_at") or j.get("created_at")
        if j.get("best_affinity") is None:
            j["best_affinity"] = rd.get("best_affinity")
        jobs_out.append(j)
    return {"total": total, "limit": limit, "offset": offset, "jobs": jobs_out}


@app.get("/api/history/{job_id}")
async def get_history_job(job_id: str):
    """Get a single persisted job by ID."""
    conn = _get_db()
    row = conn.execute("SELECT * FROM jobs WHERE id = ?", (job_id,)).fetchone()
    conn.close()
    if not row:
        raise HTTPException(404, "Job not found in history")
    result = dict(row)
    for field in ("result_data", "ligand_inputs"):
        if result.get(field):
            try:
                result[field] = json.loads(result[field])
            except Exception:
                pass
    return result


@app.delete("/api/history/{job_id}")
async def delete_history_job(job_id: str, remove_files: bool = Query(False)):
    """Delete a job from SQLite; optionally remove its run directory."""
    conn = _get_db()
    row = conn.execute("SELECT * FROM jobs WHERE id = ?", (job_id,)).fetchone()
    if not row:
        conn.close()
        raise HTTPException(404, "Job not found in history")
    run_dir = row["run_directory"]
    conn.execute("DELETE FROM jobs WHERE id = ?", (job_id,))
    conn.commit()
    conn.close()
    files_deleted = False
    if remove_files and run_dir:
        rd = Path(run_dir)
        if rd.exists():
            shutil.rmtree(rd, ignore_errors=True)
            files_deleted = True
    return {"deleted": job_id, "files_deleted": files_deleted}


@app.put("/api/history/{job_id}/name")
async def rename_history_job(job_id: str, body: JobRename):
    """Rename a persisted job."""
    conn = _get_db()
    row = conn.execute("SELECT id FROM jobs WHERE id = ?", (job_id,)).fetchone()
    if not row:
        conn.close()
        raise HTTPException(404, "Job not found in history")
    conn.execute("UPDATE jobs SET name = ? WHERE id = ?", (body.name, job_id))
    conn.commit()
    conn.close()
    if job_id in jobs and hasattr(jobs[job_id], "job_name"):
        jobs[job_id].job_name = body.name
    return {"id": job_id, "name": body.name}


@app.patch("/api/jobs/{job_id}")
async def rename_job(job_id: str, body: JobRename):
    """Frontend HistoryPage rename alias."""
    return await rename_history_job(job_id, body)


@app.delete("/api/jobs/{job_id}")
async def delete_job(job_id: str, remove_files: bool = Query(False)):
    """Frontend HistoryPage delete alias."""
    jobs.pop(job_id, None)
    job_logs.pop(job_id, None)
    return await delete_history_job(job_id, remove_files)


@app.get("/api/runs")
async def list_runs():
    """List run directories — split into temporary and permanent (saved) runs."""
    temp      = _scan_run_dirs(Path("temp_runs"))
    permanent = _scan_run_dirs(Path("data/results"))
    return {"runs": temp, "temporary_runs": temp, "permanent_runs": permanent}


@app.get("/api/runs/benchmark_results")
async def list_benchmark_results():
    """List benchmark report files (CSV / JSON / log)."""
    bdir = Path("benchmark_results")
    if not bdir.exists():
        return {"files": []}
    files = []
    for f in sorted(bdir.rglob("*")):
        if f.is_file() and f.suffix.lower() in (".csv", ".json", ".log"):
            st = f.stat()
            files.append({
                "name": f.name,
                "path": str(f),
                "size": st.st_size,
                "size_human": _human_size(st.st_size),
                "date": datetime.fromtimestamp(st.st_mtime).isoformat(),
                "file_count": 1,
                "total_size": st.st_size,
                "download_url": f"/api/download?path={f}",
            })
    return {"files": files}


@app.get("/api/runs/{run_name}/result.json")
async def get_run_result(run_name: str):
    """Return the result.json / summary.json payload for a run."""
    for root in (Path("temp_runs"), Path("data/results")):
        run_dir = root / run_name
        if not run_dir.exists():
            continue
        for p in ("run/results.json", "results.json", "summary.json"):
            fp = run_dir / p
            if fp.exists():
                return json.loads(fp.read_text(encoding="utf-8"))
    raise HTTPException(404, "Results file not found")


@app.get("/api/runs/{run_name}")
async def get_run(run_name: str):
    """List files inside a specific run directory (temp or permanent)."""
    run_dir = Path("temp_runs") / run_name
    if not run_dir.exists() or not run_dir.is_dir():
        run_dir = Path("data/results") / run_name
    if not run_dir.exists() or not run_dir.is_dir():
        raise HTTPException(404, "Run not found")
    files = []
    for f in sorted(run_dir.rglob("*")):
        if f.is_file():
            st = f.stat()
            files.append({
                "path": str(f.relative_to(run_dir)),
                "name": f.name,
                "size": st.st_size,
                "size_bytes": st.st_size,
                "size_human": _human_size(st.st_size),
                "modified": datetime.fromtimestamp(st.st_mtime).isoformat(),
                "download_url": f"/api/download?path={f}",
            })
    return {"name": run_name, "path": str(run_dir.resolve()), "files": files}


@app.delete("/api/runs/{run_name}")
async def delete_run(run_name: str):
    """Delete an entire run directory (temp or permanent)."""
    run_dir = Path("temp_runs") / run_name
    if not run_dir.exists() or not run_dir.is_dir():
        run_dir = Path("data/results") / run_name
    if not run_dir.exists() or not run_dir.is_dir():
        raise HTTPException(404, "Run not found")
    shutil.rmtree(run_dir, ignore_errors=True)
    return {"deleted": run_name}


@app.get("/api/settings")
async def get_settings():
    """Get all persisted user settings, merged with app defaults."""
    conn = _get_db()
    rows = conn.execute("SELECT key, value FROM settings").fetchall()
    conn.close()
    stored = {}
    for r in rows:
        try:
            stored[r["key"]] = json.loads(r["value"])
        except Exception:
            stored[r["key"]] = r["value"]
    merged = _deep_merge(_DEFAULT_SETTINGS, stored)
    return {"settings": merged}


@app.put("/api/settings")
async def update_settings(body: SettingsPayload):
    """Update settings (upsert key-value pairs). Accepts {settings:{...}} or a raw object."""
    payload = body.settings
    if not isinstance(payload, dict):
        raise HTTPException(400, "settings must be an object")
    if set(payload.keys()) == {"settings"} and isinstance(payload.get("settings"), dict):
        payload = payload["settings"]
    conn = _get_db()
    for key, value in payload.items():
        conn.execute(
            "INSERT INTO settings (key, value) VALUES (?, ?) "
            "ON CONFLICT(key) DO UPDATE SET value = excluded.value",
            (key, json.dumps(value, default=str)),
        )
    conn.commit()
    conn.close()
    return {"updated": list(payload.keys())}


@app.get("/api/stats")
async def get_stats():
    """Aggregate statistics across all persisted jobs."""
    conn = _get_db()
    row = conn.execute("""
        SELECT
            COUNT(*)                                                   AS total,
            SUM(CASE WHEN status = 'completed' THEN 1 ELSE 0 END)     AS completed,
            SUM(CASE WHEN status = 'failed'    THEN 1 ELSE 0 END)     AS failed,
            AVG(CASE WHEN status = 'completed' THEN best_affinity END) AS avg_affinity,
            AVG(CASE WHEN status = 'completed' THEN duration_s   END)  AS avg_duration,
            AVG(CASE WHEN status = 'completed' THEN best_rmsd    END)  AS avg_rmsd
        FROM jobs
    """).fetchone()
    conn.close()
    total     = row["total"] or 0
    completed = row["completed"] or 0
    failed    = row["failed"] or 0
    return {
        "total_jobs":     total,
        "completed":      completed,
        "failed":         failed,
        "completed_count": completed,
        "failed_count":    failed,
        "running_count":   max(0, total - completed - failed),
        "avg_rmsd":       round(row["avg_rmsd"] or 0.0, 3),
        "avg_duration":   round(row["avg_duration"] or 0.0, 2),
        "avg_affinity":   round(row["avg_affinity"] or 0.0, 3),
        "success_rate":   round(completed / total * 100, 1) if total > 0 else 0.0,
    }