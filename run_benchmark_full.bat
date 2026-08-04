@echo off
set PYTHONIOENCODING=utf-8
"C:\Users\pappu\miniconda3\envs\dippydock_env\python.exe" "D:\Dippy-dikky-dock\backend\benchmark\benchmark_runner.py" --output "D:\Dippy-dikky-dock\temp_runs\astex_full_85" --exhaustiveness 4 --workers 4 > "D:\Dippy-dikky-dock\temp_runs\astex_full_85\run_output.log" 2>&1
echo Exit code: %ERRORLEVEL% >> "D:\Dippy-dikky-dock\temp_runs\astex_full_85\run_output.log"
