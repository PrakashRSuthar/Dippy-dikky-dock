@echo off
echo ============================================
echo   DippyDock - First Time Setup
echo ============================================
echo.

REM Check Python
where python >nul 2>&1
if %errorlevel% neq 0 (
    echo [ERROR] Python not found. 
    echo Please install Python 3.10+ from https://www.python.org/downloads/
    echo IMPORTANT: Check "Add Python to PATH" during installation!
    pause
    exit /b 1
)

REM Check Node.js
where node >nul 2>&1
if %errorlevel% neq 0 (
    echo [ERROR] Node.js not found.
    echo Please install Node.js 18+ from https://nodejs.org/
    pause
    exit /b 1
)

echo [1/5] Checking Python version...
python --version

echo.
echo [2/5] Installing Python dependencies...
pip install -r requirements.txt
if %errorlevel% neq 0 (
    echo [WARNING] Some dependencies may have failed. Trying with --user flag...
    pip install -r requirements.txt --user
)

echo.
echo [3/5] Installing frontend dependencies...
cd /d "%~dp0frontend"
call npm install

echo.
echo [4/5] Building frontend...
call npm run build

echo.
echo [5/5] Setup complete!
echo.
echo ============================================
echo   Setup Complete! 
echo   
echo   To start DippyDock:
echo     double-click "start_desktop.bat"
echo   
echo   Or run manually:
echo     cd backend
echo     python -m uvicorn api.docking_api:app --port 8000
echo     
echo     cd frontend  
echo     npm run dev
echo   
echo   Then open: http://localhost:5173
echo ============================================
echo.
pause
