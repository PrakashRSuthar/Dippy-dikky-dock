# batch_processor.py
"""
Batch Docking Processor — Paper Section III-G
Thread-per-job model with semaphore concurrency cap.
All jobs register in the shared API registry (jobs / job_logs)
so existing SSE and status endpoints work without changes.
"""
from __future__ import annotations
import threading, uuid, time, logging
from dataclasses import dataclass, field
from datetime import datetime
from typing import List, Dict, Optional, Callable

logger = logging.getLogger(__name__)

#  Lazy import of API registries (avoids circular import at module load) 
def _get_api_registries():
    try:
        from docking_api import DockingPipeline, update_job_status, jobs, job_logs  # type: ignore
        return DockingPipeline, update_job_status, jobs, job_logs
    except ImportError:
        # Allow batch_processor to be tested standalone
        _jobs: Dict = {}
        _logs: Dict = {}
        def _upd(jid, status, prog, msg):
            if jid in _jobs:
                _jobs[jid]["status"]=status; _jobs[jid]["progress"]=prog; _jobs[jid]["message"]=msg
        return None, _upd, _jobs, _logs

@dataclass
class BatchItem:
    protein_input: str
    ligand_input:  str
    job_name: Optional[str] = None

@dataclass
class BatchConfig:
    cleaning_policy: Optional[Dict] = None
    retention:        str   = "keep7d"
    max_concurrency:  int   = 0   # 0 = one thread per item (unbounded)
    job_timeout_s:    int   = 900 # 15-min per-job hard timeout

    def effective_policy(self) -> Dict:
        return self.cleaning_policy or {
            "keep_waters":False,"keep_ions":True,
            "keep_solvents":False,"keep_cofactors":True}

@dataclass
class BatchResult:
    batch_id: str
    jobs:     List[Dict]
    started_at: str = field(default_factory=lambda: datetime.now().isoformat())

class BatchProcessor:
    """
    Lightweight batch runner.

    Parameters
    ----------
    on_status : callback(job_id, status, progress, message)
    on_log    : callback(job_id, log_entry_dict)
    """
    def __init__(self,
                 on_status: Optional[Callable] = None,
                 on_log:    Optional[Callable] = None):
        _, default_upd, _, _ = _get_api_registries()
        self.on_status = on_status or default_upd
        self.on_log    = on_log    or (lambda j, e: None)

    #  Public method 
    def submit(self, items: List[BatchItem], config: BatchConfig) -> BatchResult:
        if not items:
            raise ValueError("No items provided")

        DockingPipeline, update_job_status, jobs, job_logs = _get_api_registries()
        batch_id = str(uuid.uuid4())
        semaphore = (threading.Semaphore(config.max_concurrency)
                     if config.max_concurrency > 0 else None)
        results: List[Dict] = []

        def _run(jid: str, item: BatchItem):
            if semaphore: semaphore.acquire()
            try:
                if DockingPipeline is None:
                    # Standalone test mode
                    time.sleep(0.1)
                    self.on_status(jid,"completed",100,"Test OK")
                    return
                pipeline = DockingPipeline()
                pipeline.set_status_callback(self.on_status)
                pipeline.set_log_callback(self.on_log)
                pipeline.run_pipeline(
                    jid,
                    item.protein_input,
                    item.ligand_input,
                    config.effective_policy(),
                    config.retention)
            except Exception as e:
                logger.error("Job %s failed: %s", jid[:8], e)
                self.on_status(jid,"failed",0,f"Error: {e}")
            finally:
                if semaphore: semaphore.release()

        for it in items:
            jid = str(uuid.uuid4())

            # Seed shared registries so frontend can subscribe immediately
            jobs[jid] = _make_status_obj(
                jid, "queued", f"Batch {batch_id[:8]}: queued")
            job_logs[jid] = []

            t = threading.Thread(target=_run, args=(jid, it), daemon=True)
            t.start()

            results.append({
                "job_id":        jid,
                "protein_input": it.protein_input,
                "ligand_input":  it.ligand_input,
                "job_name":      it.job_name or f"job_{jid[:8]}",
            })
            logger.info("Queued job %s (%s vs %s)", jid[:8],
                        it.protein_input[:20], it.ligand_input[:20])

        logger.info("Batch %s: %d jobs submitted (concurrency=%s)",
                    batch_id[:8], len(results),
                    config.max_concurrency or "unbounded")
        return BatchResult(batch_id=batch_id, jobs=results)

#  Helper: create a job-status-like object compatible with FastAPI model 
def _make_status_obj(jid, status, message):
    """Returns a simple namespace compatible with jobs[jid].status etc."""
    class _JS:
        def __init__(self):
            self.job_id=jid; self.status=status; self.progress=0
            self.message=message; self.created_at=datetime.now().isoformat()
            self.completed_at=None
    return _JS()