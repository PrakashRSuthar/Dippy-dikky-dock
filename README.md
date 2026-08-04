# DippyDock - Molecular Docking Desktop Application

A modern molecular docking tool with a web-based GUI for protein-ligand docking studies.

## Quick Start (Students)

### Prerequisites
- **Python 3.10+** - [Download](https://www.python.org/downloads/)
- **Node.js 18+** - [Download](https://nodejs.org/)
- **Git** - [Download](https://git-scm.com/)

### Installation

1. **Clone or download this repository**
   ```
   git clone <repository-url>
   cd Dippy-dikky-dock
   ```

2. **Run the setup script** (Windows)
   ```
   setup.bat
   ```

   Or manually:
   ```bash
   # Install Python dependencies
   pip install -r requirements.txt

   # Install frontend dependencies
   cd frontend
   npm install
   npm run build
   ```

3. **Start DippyDock**
   ```
   start_desktop.bat
   ```

4. **Open your browser** to `http://localhost:5173`

## Features

### Docking
- **Protein Input**: PDB ID lookup or file upload (PDB/PDBQT)
- **Ligand Input**: PubChem search, SMILES, or file upload (SDF/MOL/MOL2/PDBQT)
- **Multi-ligand**: Batch docking with multiple ligands
- **Auto-pocket detection**: Automatic binding site identification
- **Real-time progress**: Live logs and status updates

### Results
- **3D Visualization**: Interactive molecular viewer with 3Dmol.js
- **Pose Analysis**: Affinity rankings, RMSD values, binding poses
- **Export**: Download results in various formats (PDBQT, CSV, PNG)

### History & Management
- **Job History**: Browse all past docking runs
- **Run Manager**: Explore temporary and benchmark runs
- **Settings**: Configure docking parameters, cleaning policies, retention

### Desktop App (Optional)
For a native desktop experience with Electron:
```bash
cd frontend
npm run dev:electron
```

## Project Structure

```
Dippy-dikky-dock/
├── backend/                 # Python FastAPI backend
│   ├── main.py             # API server entry point
│   ├── modules/            # Core docking modules
│   │   ├── docking_api.py  # REST API endpoints
│   │   ├── docking_engine.py # AutoDock Vina wrapper
│   │   ├── protein_prep.py # Protein preparation
│   │   └── pocket_identifier.py # Binding site detection
│   └── benchmark/          # Benchmark tools
│       ├── benchmark_runner.py
│       └── paper_gen.py
├── frontend/               # React TypeScript frontend
│   ├── src/
│   │   ├── pages/         # Route pages
│   │   │   ├── HistoryPage.tsx
│   │   │   ├── RunManagerPage.tsx
│   │   │   ├── SettingsPage.tsx
│   │   │   └── ResultsPage.tsx
│   │   ├── components/    # UI components
│   │   │   ├── ProteinUpload.tsx
│   │   │   ├── JobProgress.tsx
│   │   │   ├── ResultsChart.tsx
│   │   │   ├── ResultsTable.tsx
│   │   │   ├── MoleculeVisualization.tsx
│   │   │   └── Sidebar.tsx
│   │   └── hooks/         # React hooks
│   └── electron/          # Electron desktop wrapper
├── data/                   # Reference data
├── temp_runs/             # Docking results (auto-created)
└── requirements.txt       # Python dependencies
```

## API Reference

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/dock` | POST | Start a docking job |
| `/api/dock/batch` | POST | Start batch docking |
| `/api/jobs/{id}/status` | GET | Get job status |
| `/api/jobs/{id}/result` | GET | Get job results |
| `/api/jobs/{id}/logs` | GET | Stream job logs (SSE) |
| `/api/history` | GET | List job history |
| `/api/runs` | GET | List temp runs |
| `/api/settings` | GET/PUT | Manage settings |
| `/api/stats` | GET | Dashboard statistics |

## Configuration

Settings are persisted in `data/dippydock.db` (SQLite):

- **Docking Engine**: exhaustiveness (1-32), num_modes, energy_range
- **Box Settings**: auto-box, min/max axis sizes
- **Cleaning Policy**: keep waters/ions/solvents/cofactors
- **Retention**: permanent, 7 days, or delete immediately

## Troubleshooting

### "Port already in use"
```bash
# Kill process on port 8000
netstat -ano | findstr :8000
taskkill /PID <PID> /F
```

### "Module not found"
```bash
pip install -r requirements.txt
```

### "Frontend build failed"
```bash
cd frontend
rm -rf node_modules
npm install
npm run build
```

## License

Academic use only. See LICENSE file for details.

## Acknowledgments

- [AutoDock Vina](http://vina.scripps.edu/) - Docking engine
- [3Dmol.js](https://3dmol.csb.pitt.edu/) - Molecular visualization
- [React](https://react.dev/) - UI framework
- [FastAPI](https://fastapi.tiangolo.com/) - Backend framework
