// electron/install.js
// Script to install Python dependencies when app is first launched
const { exec } = require('child_process');
const path = require('path');
const fs = require('fs');

const isDev = !app?.isPackaged;

function getPythonPath() {
  if (isDev) return 'python';
  return path.join(process.resourcesPath, 'python', 'python.exe');
}

function getBackendPath() {
  if (isDev) return path.join(__dirname, '..', '..', 'backend');
  return path.join(process.resourcesPath, 'backend');
}

function installDependencies(callback) {
  const pythonPath = getPythonPath();
  const requirementsPath = path.join(getBackendPath(), '..', 'requirements.txt');

  if (!fs.existsSync(requirementsPath)) {
    console.log('No requirements.txt found, skipping install');
    callback(null);
    return;
  }

  console.log('Installing Python dependencies...');
  exec(`"${pythonPath}" -m pip install -r "${requirementsPath}" --quiet`, (error, stdout, stderr) => {
    if (error) {
      console.error('Failed to install dependencies:', error);
      callback(error);
      return;
    }
    console.log('Dependencies installed successfully');
    callback(null);
  });
}

module.exports = { installDependencies };
