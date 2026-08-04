"""DippyDock backend — compatibility entry point.

Historically this module defined its own, divergent FastAPI application
(backed by ``docking_pipeline.CleanDockingPipeline``). That duplicated server
only implemented a subset of the endpoints the frontend calls, which is why
logs / runs / results pages broke whenever the backend was launched via
``main:app`` (e.g. ``start_desktop.bat``, ``start_desktop.ps1``).

The canonical application now lives in ``api.docking_api``. This module
simply re-exports it so that ``python -m uvicorn main:app ...`` and
``python -m uvicorn api.docking_api:app ...`` serve the exact same API.
"""
import sys
from pathlib import Path

# Force UTF-8 output: Windows defaults stdout to cp1252 ("charmap"), which cannot
# encode the emoji/box-drawing characters used in logs and crashes the pipeline.
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

sys.path.insert(0, str(Path(__file__).parent))

from api.docking_api import app, jobs  # noqa: E402

__all__ = ["app", "jobs"]


if __name__ == "__main__":
    import uvicorn

    uvicorn.run("main:app", host="0.0.0.0", port=8000, reload=False)
