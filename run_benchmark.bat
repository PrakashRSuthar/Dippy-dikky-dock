@echo off
SETLOCAL ENABLEDELAYEDEXPANSION
echo.
echo === Dippy-Dikky-Dock — Astex Benchmark Runner ===
echo.
set PYTHONIOENCODING=utf-8
"C:\Users\pappu\miniconda3\envs\dippydock_env\python.exe" backend/benchmark/benchmark_runner.py %*
echo.
echo Done.
