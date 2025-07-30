# backend/modules/docking_engine.py

import os
import subprocess
from pathlib import Path
import re

# ... (The DockingEngine class remains the same as before)
SCRIPT_DIR = Path(__file__).parent
VINA_EXECUTABLE = SCRIPT_DIR.parent / "tools" / "vina" / "vina.exe"

class DockingEngine:
    def __init__(self, protein_pdbqt_path: str, ligand_pdbqt_path: str):
        self.protein_path = Path(protein_pdbqt_path)
        self.ligand_path = Path(ligand_pdbqt_path)
        if not self.protein_path.exists() or not self.ligand_path.exists():
            raise FileNotFoundError("Prepared protein or ligand PDBQT file not found.")

    def run_docking(self, center_x: float, center_y: float, center_z: float,
                    box_size_x: float = 20.0, box_size_y: float = 20.0, box_size_z: float = 20.0,
                    output_dir: str = "data/docking_results"):
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        protein_name = self.protein_path.stem.replace('_prepared', '')
        ligand_name = self.ligand_path.stem.replace('_prepared', '')
        output_poses_path = output_path / f"{protein_name}_{ligand_name}_docked.pdbqt"
        print("[INFO] 🧬 Starting docking simulation via command line...")
        if not VINA_EXECUTABLE.is_file():
            print(f"❌ ERROR: vina.exe not found at {VINA_EXECUTABLE}")
            return None
        cmd = [
            str(VINA_EXECUTABLE), '--receptor', str(self.protein_path), '--ligand', str(self.ligand_path),
            '--out', str(output_poses_path), '--center_x', str(center_x), '--center_y', str(center_y),
            '--center_z', str(center_z), '--size_x', str(box_size_x), '--size_y', str(box_size_y),
            '--size_z', str(box_size_z), '--exhaustiveness', '8'
        ]
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            print(f"✅ Docking complete. Results saved to: {output_poses_path}")
            print("\n[INFO] Binding Affinity Scores (kcal/mol):")
            for line in result.stdout.splitlines():
                if re.match(r'^\s*\d+\s+', line):
                    parts = line.split()
                    affinity = float(parts[1])
                    print(f"   Mode {parts[0]}: {affinity:.4f}")
            return str(output_poses_path)
        except subprocess.CalledProcessError as e:
            print("❌ ERROR: AutoDock Vina failed to run.")
            print(f"Stderr: {e.stderr}")
            return None

if __name__ == '__main__':
    print("--- Running docking engine test ---")
    
    protein_file = input("Enter the path to the prepared protein PDBQT file: ").strip()
    ligand_file = input("Enter the path to the prepared ligand PDBQT file: ").strip()

    try:
        print("Enter the center coordinates for the docking box:")
        # --- SIMPLIFICATION ---
        # User input is immediately converted to a standard Python float
        center_x = float(input("Center X: ").strip())
        center_y = float(input("Center Y: ").strip())
        center_z = float(input("Center Z: ").strip())
        
        pocket_center = {"center_x": center_x, "center_y": center_y, "center_z": center_z}

        if not Path(protein_file).exists() or not Path(ligand_file).exists():
            print(f"❌ ERROR: One or both input files not found. Please check the paths.")
        else:
            engine = DockingEngine(protein_pdbqt_path=protein_file, ligand_pdbqt_path=ligand_file)
            engine.run_docking(**pocket_center)
    
    except ValueError:
        print("❌ ERROR: Invalid input. Please enter a valid number for the coordinates.")
    except Exception as e:
        print(f"❌ An unexpected error occurred: {e}")
        
    print("\n--- Test complete ---")
    