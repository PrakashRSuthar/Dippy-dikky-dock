@echo off
title Dippy-Dikky-Dock Application Launcher
rem ====================================================================
#  Startup script for the Dippy-Dikky-Dock Molecular Docking App
# ====================================================================

echo =================================================================
echo       🧬  Dippy-Dikky-Dock Application Launcher  🧬        
echo =================================================================
echo  This launcher will activate the environment and start both
echo  the backend FastAPI server and React frontend dev server.
echo =================================================================
echo.

rem --- 1. Find Conda executable ---
set "CONDA_BIN="
where conda >nul 2>nul
if %errorlevel% equ 0 (
    set "CONDA_BIN=conda"
) else (
    rem Check standard folders
    if exist "%USERPROFILE%\miniconda3\condabin\conda.bat" (
        set "CONDA_BIN=%USERPROFILE%\miniconda3\condabin\conda.bat"
    ) else if exist "%USERPROFILE%\Anaconda3\condabin\conda.bat" (
        set "CONDA_BIN=%USERPROFILE%\Anaconda3\condabin\conda.bat"
    ) else if exist "C:\miniconda3\condabin\conda.bat" (
        set "CONDA_BIN=C:\miniconda3\condabin\conda.bat"
    ) else if exist "C:\ProgramData\miniconda3\condabin\conda.bat" (
        set "CONDA_BIN=C:\ProgramData\miniconda3\condabin\conda.bat"
    )
)

if "%CONDA_BIN%"=="" (
    echo [ERROR] Conda not found! Please run setup_win.ps1 first.
    echo.
    pause
    exit /b 1
)

rem --- 2. Find npm executable ---
set "NPM_BIN="
where npm >nul 2>nul
if %errorlevel% equ 0 (
    set "NPM_BIN=npm"
) else (
    if exist "C:\Program Files\nodejs\npm.cmd" (
        set "NPM_BIN=C:\Program Files\nodejs\npm.cmd"
    )
)

if "%NPM_BIN%"=="" (
    echo [ERROR] Node.js/npm not found! Please run setup_win.ps1 first.
    echo.
    pause
    exit /b 1
)

echo ✅ Conda found: %CONDA_BIN%
echo ✅ npm found: %NPM_BIN%
echo.
echo 🚀 Launching API Backend Server (FastAPI)...
echo ----------------------------------------------------
start "Dippy Docking API Backend" cmd /c "call "%CONDA_BIN%" activate dippydock_env && python backend/api/docking_api.py"

echo 🚀 Launching React Web Frontend (Vite)...
echo ----------------------------------------------------
start "Dippy Docking Web Frontend" cmd /c "cd frontend && "%NPM_BIN%" run dev"

echo.
echo =================================================================
echo 🎉  BOTH SERVICES LAUNCHED SUCCESSFULLY!
echo =================================================================
echo  • Web Frontend: http://localhost:5173
echo  • Backend API:  http://localhost:8000
echo.
echo  Keep the separate console windows open while using the app.
echo  To shut down, simply close the backend and frontend windows.
echo =================================================================
echo.
pause
