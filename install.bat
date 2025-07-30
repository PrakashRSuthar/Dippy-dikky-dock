@echo off
REM ================================================================
REM  Installer for the Dippy-Dikky-Dock Project
REM ================================================================

REM --- Check if conda command is available ---
conda --version >nul 2>nul
if %errorlevel% neq 0 (
    echo.
    echo ERROR: Conda command not found.
    echo Please make sure you have installed Miniconda/Anaconda
    echo and that it's accessible from your terminal.
    echo.
    pause
    exit /b 1
)

REM --- Create the environment from the .yml file ---
echo.
echo Found Conda. Now creating the 'dippy' environment...
echo This might take several minutes. Please be patient.
echo.
conda env create -f environment.yml

REM --- Final Instructions ---
echo.
echo =================================================================
echo Setup is complete!
echo To activate the new environment, open a new terminal and run:
echo   conda activate dippy
echo =================================================================
echo.
pause