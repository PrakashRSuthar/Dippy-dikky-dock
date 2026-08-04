$benchmarkDir = "D:\Dippy-dikky-dock\temp_runs\astex_full_85"
$logFile = Join-Path $benchmarkDir "run.log"
$errFile = Join-Path $benchmarkDir "run_err.log"
$python = "C:\Users\pappu\miniconda3\envs\dippydock_env\python.exe"
$script = "D:\Dippy-dikky-dock\backend\benchmark\benchmark_runner.py"

New-Item -ItemType Directory -Path $benchmarkDir -Force | Out-Null

$env:PYTHONIOENCODING = "utf-8"

$psi = New-Object System.Diagnostics.ProcessStartInfo
$psi.FileName = $python
$psi.Arguments = "D:\Dippy-dikky-dock\backend\benchmark\benchmark_runner.py --output $benchmarkDir --exhaustiveness 8"
$psi.RedirectStandardOutput = $true
$psi.RedirectStandardError = $true
$psi.UseShellExecute = $false
$psi.CreateNoWindow = $true
$psi.WorkingDirectory = "D:\Dippy-dikky-dock"
$p = New-Object System.Diagnostics.Process
$p.StartInfo = $psi
$p.Start() | Out-Null

$stdout = $p.StandardOutput.ReadToEnd()
$stderr = $p.StandardError.ReadToEnd()

$stdout | Out-File $logFile
$stderr | Out-File $errFile

Write-Output "PID: $($p.Id)"
