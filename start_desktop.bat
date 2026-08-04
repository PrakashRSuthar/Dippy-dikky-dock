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
echo  [1/3] Checking Python dependencies...
pip install -r "%~dp0requirements.txt" --quiet --disable-pip-version-check 2>nul

echo  [2/3] Checking frontend dependencies...
cd /d "%~dp0frontend"
call npm install --quiet 2>nul

echo  [3/3] Launching DippyDock desktop app...
echo.
echo  The app window will open shortly. Closing the app stops the backend.
echo.

call npm run dev:electron

echo.
echo  DippyDock closed. Stopping backend service...
taskkill /fi "WINDOWTITLE eq DippyDock-Backend*" /f >nul 2>&1
for /f "tokens=5" %%a in ('netstat -ano ^| findstr ":8000" ^| findstr "LISTENING"') do taskkill /pid %%a /f >nul 2>&1
echo  Done.
echo.
pause
