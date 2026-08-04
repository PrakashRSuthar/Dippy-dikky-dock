@echo off
SETLOCAL ENABLEDELAYEDEXPANSION
if "%1"=="" (
    set "BENCH_DIR=D:\Dippy-dikky-dock\temp_runs\astex_full_85"
) else (
    set "BENCH_DIR=%1"
)
echo === Astex Benchmark Monitor ===
echo Output: %BENCH_DIR%
echo.

REM Count completed
set COUNT=0
if exist "%BENCH_DIR%\results" (
    for /d %%d in ("%BENCH_DIR%\results\*") do (
        if exist "%%d\result.json" set /a COUNT+=1
    )
)
echo Completed targets: %COUNT%/85

REM Check latest
echo.
echo Latest results:
if exist "%BENCH_DIR%\results" (
    for /d %%d in ("%BENCH_DIR%\results\*") do (
        if exist "%%d\result.json" (
            for /f "tokens=2 delims=:," %%a in ('findstr "rmsd" "%%d\result.json"') do (
                set "RMSD=%%a"
            )
            for /f "tokens=2 delims=:," %%b in ('findstr "success" "%%d\result.json"') do (
                set "OK=%%b"
            )
            echo   %%~nxd -- RMSD=!RMSD:~0,7! success=!OK!
        )
    )
)

REM Check log tail
echo.
echo Recent log entries:
if exist "%BENCH_DIR%\run_output.log" (
    powershell -Command "Get-Content '%BENCH_DIR%\run_output.log' -Tail 5"
)

echo.
echo Use: Get-Content "%BENCH_DIR%\run_output.log" -Wait  (to follow)
