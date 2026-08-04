// electron/main.js
const { app, BrowserWindow, shell } = require('electron');
const path = require('path');
const { spawn } = require('child_process');

const isDev = !app.isPackaged;
let mainWindow;
let backendProcess;

function startBackend() {
  const backendPath = isDev
    ? path.join(__dirname, '..', '..', 'backend', 'main.py')
    : path.join(process.resourcesPath, 'backend', 'main.py');

  const pythonPath = isDev
    ? 'python'
    : path.join(process.resourcesPath, 'python', 'python.exe');

  console.log('Starting backend:', pythonPath, backendPath);

  backendProcess = spawn(pythonPath, [backendPath], {
    cwd: isDev
      ? path.join(__dirname, '..', '..', 'backend')
      : path.join(process.resourcesPath, 'backend'),
    env: {
      ...process.env,
      PYTHONPATH: isDev
        ? path.join(__dirname, '..', '..', 'backend')
        : path.join(process.resourcesPath, 'backend'),
    },
  });

  backendProcess.stdout.on('data', (data) => {
    console.log(`Backend: ${data.toString()}`);
  });

  backendProcess.stderr.on('data', (data) => {
    console.error(`Backend error: ${data.toString()}`);
  });

  backendProcess.on('close', (code) => {
    console.log(`Backend exited with code ${code}`);
  });
}

function createWindow() {
  mainWindow = new BrowserWindow({
    width: 1400,
    height: 900,
    minWidth: 1024,
    minHeight: 700,
    title: 'DippyDock',
    icon: path.join(__dirname, '..', 'public', 'icon.png'),
    webPreferences: {
      nodeIntegration: false,
      contextIsolation: true,
    },
    titleBarStyle: 'hiddenInset',
    backgroundColor: '#030712',
  });

  if (isDev) {
    mainWindow.loadURL('http://localhost:5173');
    mainWindow.webContents.openDevTools();
  } else {
    mainWindow.loadFile(path.join(__dirname, '..', 'dist', 'index.html'));
  }

  mainWindow.setWindowButtonVisibility?.(true);

  mainWindow.webContents.setWindowOpenHandler(({ url }) => {
    shell.openExternal(url);
    return { action: 'deny' };
  });

  mainWindow.on('closed', () => {
    mainWindow = null;
  });
}

app.whenReady().then(() => {
  startBackend();

  // Wait for backend to start
  setTimeout(() => {
    createWindow();
  }, 2000);

  app.on('activate', () => {
    if (BrowserWindow.getAllWindows().length === 0) {
      createWindow();
    }
  });
});

app.on('window-all-closed', () => {
  if (backendProcess) {
    backendProcess.kill();
  }
  if (process.platform !== 'darwin') {
    app.quit();
  }
});

app.on('before-quit', () => {
  if (backendProcess) {
    backendProcess.kill();
  }
});
