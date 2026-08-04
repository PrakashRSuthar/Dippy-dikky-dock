# setup_win.ps1
# ====================================================================
#  Intelligent Installer for the Dippy-Dikky-Dock Molecular Docking App
# ====================================================================

$ErrorActionPreference = "Stop"

Write-Host "=================================================================" -ForegroundColor Cyan
Write-Host "       🧬  Dippy-Dikky-Dock Dependency Setup Script  🧬        " -ForegroundColor Cyan
Write-Host "=================================================================" -ForegroundColor Cyan
Write-Host "  This script will automatically detect and install all necessary"
Write-Host "  components, including Miniconda, Node.js, Conda environments,  "
Write-Host "  and frontend libraries. Please wait...                         "
Write-Host "=================================================================" -ForegroundColor Cyan
Write-Host ""

# --------------------------------------------------------------------
# 1. Detect or Install Miniconda
# --------------------------------------------------------------------
Write-Host "[1/4] Checking for Conda installation..." -ForegroundColor Yellow

$condaBin = $null

# Check if conda is already in PATH
$condaInPath = Get-Command conda -ErrorAction SilentlyContinue
if ($condaInPath) {
    $condaBin = "conda"
    Write-Host "  ✅ Conda detected in system PATH." -ForegroundColor Green
} else {
    # Check standard install locations
    $standardPaths = @(
        "$env:USERPROFILE\miniconda3\condabin\conda.bat",
        "$env:USERPROFILE\Anaconda3\condabin\conda.bat",
        "C:\miniconda3\condabin\conda.bat",
        "C:\ProgramData\miniconda3\condabin\conda.bat"
    )

    foreach ($path in $standardPaths) {
        if (Test-Path $path) {
            $condaBin = $path
            Write-Host "  ✅ Conda detected at: $path" -ForegroundColor Green
            break
        }
    }
}

if ($null -eq $condaBin) {
    Write-Host "  ⚠️  Conda not found. Downloading Miniconda installer..." -ForegroundColor DarkYellow
    $minicondaUrl = "https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe"
    $installerPath = "$env:TEMP\Miniconda3-installer.exe"
    
    Write-Host "  📥 Downloading installer to temporary folder..."
    Invoke-WebRequest -Uri $minicondaUrl -OutFile $installerPath
    
    Write-Host "  ⚙️  Installing Miniconda silently to $env:USERPROFILE\miniconda3..." -ForegroundColor Blue
    $installArgs = @("/S", "/RegisterPython=0", "/D=$env:USERPROFILE\miniconda3")
    $process = Start-Process -FilePath $installerPath -ArgumentList $installArgs -Wait -NoNewWindow -PassThru
    
    if ($process.ExitCode -eq 0) {
        $condaBin = "$env:USERPROFILE\miniconda3\condabin\conda.bat"
        Write-Host "  ✅ Miniconda installed successfully!" -ForegroundColor Green
    } else {
        Write-Error "❌ Miniconda installation failed with exit code $($process.ExitCode)"
    }
    
    # Cleanup installer
    if (Test-Path $installerPath) { Remove-Item $installerPath }
}

# --------------------------------------------------------------------
# 2. Detect or Install Node.js (npm)
# --------------------------------------------------------------------
Write-Host ""
Write-Host "[2/4] Checking for Node.js and npm..." -ForegroundColor Yellow

$npmBin = $null

# Check if npm is in PATH
$npmInPath = Get-Command npm -ErrorAction SilentlyContinue
if ($npmInPath) {
    $npmBin = "npm"
    Write-Host "  ✅ npm detected in system PATH." -ForegroundColor Green
} else {
    $nodePath = "C:\Program Files\nodejs\npm.cmd"
    if (Test-Path $nodePath) {
        $npmBin = $nodePath
        Write-Host "  ✅ npm detected at: $nodePath" -ForegroundColor Green
    }
}

if ($null -eq $npmBin) {
    Write-Host "  ⚠️  Node.js/npm not found. Downloading Node.js LTS installer..." -ForegroundColor DarkYellow
    $nodeUrl = "https://nodejs.org/dist/v20.11.1/node-v20.11.1-x64.msi"
    $nodeInstaller = "$env:TEMP\node-installer.msi"
    
    Write-Host "  📥 Downloading Node.js installer..."
    Invoke-WebRequest -Uri $nodeUrl -OutFile $nodeInstaller
    
    Write-Host "  ⚙️  Installing Node.js silently..." -ForegroundColor Blue
    $nodeArgs = @("/i", $nodeInstaller, "/qn", "/norestart")
    $process = Start-Process -FilePath msiexec.exe -ArgumentList $nodeArgs -Wait -NoNewWindow -PassThru
    
    if ($process.ExitCode -eq 0 -or $process.ExitCode -eq 3010) {
        $npmBin = "C:\Program Files\nodejs\npm.cmd"
        Write-Host "  ✅ Node.js and npm installed successfully!" -ForegroundColor Green
    } else {
        Write-Error "❌ Node.js installation failed with exit code $($process.ExitCode)"
    }
    
    # Cleanup installer
    if (Test-Path $nodeInstaller) { Remove-Item $nodeInstaller }
}

# --------------------------------------------------------------------
# 3. Create or Update Conda Environment
# --------------------------------------------------------------------
Write-Host ""
Write-Host "[3/4] Setting up Python environment via Conda..." -ForegroundColor Yellow
Write-Host "  This includes scientific libraries (RDKit, OpenBabel, Meeko, NumPy, etc.)"
Write-Host "  It may take a few minutes. Please don't close this window..." -ForegroundColor DarkGray

# Check if dippydock_env already exists
$envList = & $condaBin env list
$envExists = $envList -match "dippydock_env"

if ($envExists) {
    Write-Host "  ♻️  Conda environment 'dippydock_env' already exists. Updating packages..." -ForegroundColor Blue
    & $condaBin env update -f environment.yml --prune
} else {
    Write-Host "  🚀 Creating Conda environment 'dippydock_env' from environment.yml..." -ForegroundColor Blue
    & $condaBin env create -f environment.yml
}

Write-Host "  📦 Installing extra PIP packages..." -ForegroundColor Blue
& $condaBin run -n dippydock_env pip install -r requirements.txt

Write-Host "  ✅ Python Conda environment is ready!" -ForegroundColor Green

# --------------------------------------------------------------------
# 4. Install Frontend NPM Packages
# --------------------------------------------------------------------
Write-Host ""
Write-Host "[4/4] Setting up React frontend libraries..." -ForegroundColor Yellow

if (Test-Path "frontend") {
    Write-Host "  📂 Navigating to frontend directory..."
    Push-Location frontend
    
    Write-Host "  ⚙️  Running npm install..." -ForegroundColor Blue
    & $npmBin install
    
    Pop-Location
    Write-Host "  ✅ Frontend dependencies installed successfully!" -ForegroundColor Green
} else {
    Write-Host "  ⚠️  Frontend directory not found! Skipping..." -ForegroundColor Red
}

# --------------------------------------------------------------------
# Setup complete!
# --------------------------------------------------------------------
Write-Host ""
Write-Host "=================================================================" -ForegroundColor Green
Write-Host "🎉  SETUP COMPLETE! Dippy-Dikky-Dock is ready to run!" -ForegroundColor Green
Write-Host "=================================================================" -ForegroundColor Green
Write-Host "  To launch the entire software (API + Web Frontend):"
Write-Host "  Simply run the startup script:"
Write-Host "     .\start_win.bat" -ForegroundColor Yellow
Write-Host "=================================================================" -ForegroundColor Green
Write-Host ""

Read-Host "Press Enter to exit..."
