# DippyDock Distribution Packager
# Creates a zip file for student distribution

$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
$outputDir = "$root\dist"
$outputFile = "$outputDir\DippyDock-Student-Edition.zip"

Write-Host ""
Write-Host "Packaging DippyDock for distribution..." -ForegroundColor Cyan
Write-Host ""

# Create output directory
if (!(Test-Path $outputDir)) {
    New-Item -ItemType Directory -Path $outputDir -Force | Out-Null
}

# Files to include
$includeFiles = @(
    "README.md",
    "requirements.txt",
    "environment.yml",
    "setup.bat",
    "start_desktop.bat",
    "start_desktop.ps1",
    ".gitignore"
)

$includeDirs = @(
    "backend",
    "frontend\dist",
    "frontend\package.json",
    "frontend\vite.config.ts",
    "frontend\tsconfig.json",
    "data"
)

# Create temp staging area
$stagingDir = "$outputDir\DippyDock-staging"
if (Test-Path $stagingDir) {
    Remove-Item -Recurse -Force $stagingDir
}
New-Item -ItemType Directory -Path $stagingDir -Force | Out-Null

# Copy files
Write-Host "Copying files..." -ForegroundColor Yellow
foreach ($file in $includeFiles) {
    $src = Join-Path $root $file
    if (Test-Path $src) {
        Copy-Item $src -Destination $stagingDir -Force
        Write-Host "  + $file"
    }
}

# Copy directories
foreach ($dir in $includeDirs) {
    $src = Join-Path $root $dir
    $dst = Join-Path $stagingDir $dir
    if (Test-Path $src) {
        if ($dir -like "*\*") {
            # Nested directory - create parent first
            $parent = Split-Path $dst -Parent
            if (!(Test-Path $parent)) {
                New-Item -ItemType Directory -Path $parent -Force | Out-Null
            }
        }
        Copy-Item $src -Destination $dst -Recurse -Force
        Write-Host "  + $dir/"
    }
}

# Clean up unnecessary files from staging
Write-Host ""
Write-Host "Cleaning up..." -ForegroundColor Yellow
$cleanPatterns = @(
    "*.pyc",
    "__pycache__",
    "node_modules",
    ".git",
    "*.log",
    "temp_runs"
)

foreach ($pattern in $cleanPatterns) {
    Get-ChildItem -Path $stagingDir -Filter $pattern -Recurse -ErrorAction SilentlyContinue | 
        Remove-Item -Recurse -Force -ErrorAction SilentlyContinue
}

# Create zip
Write-Host ""
Write-Host "Creating zip archive..." -ForegroundColor Yellow
if (Test-Path $outputFile) {
    Remove-Item $outputFile -Force
}

Compress-Archive -Path "$stagingDir\*" -DestinationPath $outputFile -CompressionLevel Optimal

# Cleanup staging
Remove-Item -Recurse -Force $stagingDir

# Report
$size = (Get-Item $outputFile).Length / 1MB
Write-Host ""
Write-Host "============================================" -ForegroundColor Green
Write-Host "  Package created successfully!" -ForegroundColor Green
Write-Host "  Location: $outputFile" -ForegroundColor Green
Write-Host "  Size: $([math]::Round($size, 2)) MB" -ForegroundColor Green
Write-Host "============================================" -ForegroundColor Green
Write-Host ""
Write-Host "Students can:" -ForegroundColor Cyan
Write-Host "  1. Extract the zip" -ForegroundColor White
Write-Host "  2. Run setup.bat" -ForegroundColor White
Write-Host "  3. Run start_desktop.bat" -ForegroundColor White
Write-Host ""
