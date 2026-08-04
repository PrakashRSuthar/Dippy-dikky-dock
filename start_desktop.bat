@echo off
setlocal EnableDelayedExpansion

echo.
echo  ============================================
echo   DippyDock - Molecular Docking Desktop App
echo  ============================================
echo.

REM ── Check prerequisites ──────────────────────────────────────────────
where python >nul 2>&1
if %errorlevel% neq 0 (
    echo  [ERROR] Python not found. Install Python 3.10+ and add it to PATH.
    echo           https://www.python.org/downloads/
    pause
    exit /b 1
)

where node >nul 2>&1
if %errorlevel% neq 0 (
    echo  [ERROR] Node.js not found. Install Node.js 18+ and add it to PATH.
    echo           https://nodejs.org/
    pause
    exit /b 1
)

REM ── Install dependencies ─────────────────────────────────────────────
echo  [1/5] Checking Python dependencies...
pip install -r "%~dp0requirements.txt" --quiet --disable-pip-version-check 2>nul
if %errorlevel% neq 0 (
    echo  [WARN] Some Python packages may be missing - backend may fail to start.
)

echo  [2/5] Checking frontend dependencies...
cd /d "%~dp0frontend"
call npm install --quiet 2>nul

echo  [3/5] Building frontend assets...
call npm run build 2>nul

echo  [4/5] Starting backend server...

REM ── Start backend in a visible window so we can track it ─────────────
REM     Using "start" with a window lets uvicorn run properly with --reload.
REM     The title lets us identify and kill it on exit.
start "DippyDock-Backend" /min cmd /c "cd /d "%~dp0backend" && python -m uvicorn api.docking_api:app --host 0.0.0.0 --port 8000"

REM ── Wait for backend to respond ──────────────────────────────────────
set /a "RETRIES=0"
set /a "MAX_RETRIES=30"

:wait_for_backend
timeout /t 1 /nobreak >nul
set /a "RETRIES+=1"

curl -s http://localhost:8000/ >nul 2>&1
if %errorlevel% equ 0 (
    goto backend_ready
)

if %RETRIES% geq %MAX_RETRIES% (
    echo.
    echo  [ERROR] Backend did not respond after %MAX_RETRIES% seconds.
    echo           Check the backend window for errors.
    taskkill /fi "WINDOWTITLE eq DippyDock-Backend*" /f >nul 2>&1
    pause
    exit /b 1
)

goto wait_for_backend

:backend_ready
echo  [4/5] Backend is ready on http://localhost:8000

REM ── Start frontend ───────────────────────────────────────────────────
echo  [5/5] Starting frontend dev server...
echo.
echo  ============================================
echo   DippyDock is running!
echo   Frontend:  http://localhost:5173
echo   Backend:   http://localhost:8000
echo   Press Ctrl+C to stop
echo  ============================================
echo.

cd /d "%~dp0frontend"
call npm run dev

REM ── Cleanup on exit ──────────────────────────────────────────────────
echo.
echo  Stopping DippyDock...
taskkill /fi "WINDOWTITLE eq DippyDock-Backend*" /f >nul 2>&1
echo  Done.
