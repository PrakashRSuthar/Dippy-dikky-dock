# DippyDock Desktop Launcher
# Run this script to start DippyDock

$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $MyInvocation.MyCommand.Path

function Write-Header {
    param([string]$Message)
    Write-Host ""
    Write-Host "============================================" -ForegroundColor Cyan
    Write-Host "  $Message" -ForegroundColor Cyan
    Write-Host "============================================" -ForegroundColor Cyan
    Write-Host ""
}

function Write-Step {
    param([string]$Message)
    Write-Host "[INFO] $Message" -ForegroundColor Yellow
}

function Write-Success {
    param([string]$Message)
    Write-Host "[OK]   $Message" -ForegroundColor Green
}

function Write-Error {
    param([string]$Message)
    Write-Host "[ERROR] $Message" -ForegroundColor Red
}

function Write-Warning {
    param([string]$Message)
    Write-Host "[WARN] $Message" -ForegroundColor Yellow
}

function Test-Command {
    param([string]$Command, [string]$ErrorMessage)
    try {
        & $Command 2>$null | Out-Null
        return $true
    } catch {
        return $false
    }
}

function Wait-BackendReady {
    param(
        [string]$Url = "http://localhost:8000/",
        [int]$TimeoutSeconds = 30,
        [int]$IntervalSeconds = 1
    )
    Write-Step "Waiting for backend to be ready..."
    $elapsed = 0
    while ($elapsed -lt $TimeoutSeconds) {
        try {
            $response = Invoke-RestMethod -Uri $Url -Method Get -TimeoutSec 2 -ErrorAction Stop
            if ($response.status -eq "running") {
                Write-Success "Backend is ready at $Url"
                return $true
            }
        } catch {
            # Backend not ready yet
        }
        Start-Sleep -Seconds $IntervalSeconds
        $elapsed += $IntervalSeconds
        Write-Host "." -NoNewline -ForegroundColor Gray
    }
    Write-Error "Backend did not become ready within $TimeoutSeconds seconds"
    return $false
}

Write-Header "DippyDock - Molecular Docking Desktop App"

# Check Python
Write-Step "Checking Python installation..."
try {
    $pythonVersion = python --version 2>&1
    Write-Success "Python: $pythonVersion"
} catch {
    Write-Error "Python not found. Install Python 3.10+ from https://python.org"
    exit 1
}

# Check Node.js
Write-Step "Checking Node.js installation..."
try {
    $nodeVersion = node --version 2>&1
    Write-Success "Node.js: $nodeVersion"
} catch {
    Write-Error "Node.js not found. Install Node.js 18+ from https://nodejs.org"
    exit 1
}

# Install Python dependencies
Write-Step "Installing Python dependencies..."
Push-Location "$root\backend"
try {
    pip install -r ..\requirements.txt --quiet
    Write-Success "Python dependencies installed"
} catch {
    Write-Warning "Some Python dependencies may have failed to install"
}
Pop-Location

# Install frontend dependencies
Write-Step "Installing frontend dependencies..."
Push-Location "$root\frontend"
try {
    npm install --quiet 2>$null
    Write-Success "Frontend dependencies installed"
} catch {
    Write-Warning "Some frontend dependencies may have failed to install"
}
Pop-Location

# Build frontend
Write-Step "Building frontend..."
Push-Location "$root\frontend"
try {
    npm run build 2>$null
    Write-Success "Frontend built successfully"
} catch {
    Write-Warning "Frontend build had warnings"
}
Pop-Location

# Start backend
Write-Step "Starting backend server..."
$backendJob = Start-Job -ScriptBlock {
    param($backendPath)
    Set-Location $backendPath
    python -m uvicorn api.docking_api:app --host 0.0.0.0 --port 8000
} -ArgumentList "$root\backend"

# Wait for backend to be ready
if (-not (Wait-BackendReady)) {
    Write-Error "Failed to start backend. Check backend logs for errors."
    Stop-Job $backendJob -ErrorAction SilentlyContinue
    exit 1
}

# Start frontend (blocking)
Write-Success "Backend is ready. Starting frontend..."
Write-Host ""
Write-Header "DippyDock is running!"
Write-Host "  Frontend: http://localhost:5173" -ForegroundColor Green
Write-Host "  Backend:  http://localhost:8000" -ForegroundColor Green
Write-Host "  Press Ctrl+C to stop both servers" -ForegroundColor Green
Write-Host "============================================" -ForegroundColor Green
Write-Host ""

Push-Location "$root\frontend"
try {
    npm run dev
} finally {
    Write-Host ""
    Write-Step "Shutting down..."
    Stop-Job $backendJob -ErrorAction SilentlyContinue
    Pop-Location
}