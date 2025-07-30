# backend/modules/pocket_identifier.py

import os
import subprocess
import pandas as pd
from pathlib import Path

# --- CONFIGURATION ---
SCRIPT_DIR = Path(__file__).parent
P2RANK_EXECUTABLE = SCRIPT_DIR.parent / "tools" / "p2rank_2.4.2" / "prank.bat"

class PocketIdentifier:
    # ... (this class remains the same)
    def __init__(self, protein_path: str):
        self.protein_path = Path(protein_path)
        if not self.protein_path.exists():
            raise FileNotFoundError(f"Protein file not found: {self.protein_path}")
        self.output_dir = Path("data/p2rank_results") / self.protein_path.stem
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def find_pockets(self):
        print(f"[INFO] Running P2Rank for {self.protein_path.name}...")
        if not P2RANK_EXECUTABLE.is_file():
            print(f"❌ ERROR: P2Rank executable not found at: {P2RANK_EXECUTABLE}")
            return None
        cmd = [str(P2RANK_EXECUTABLE), "predict", "-f", str(self.protein_path), "-o", str(self.output_dir)]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            predictions_csv = self.output_dir / f"{self.protein_path.name}_predictions.csv"
            if predictions_csv.exists():
                print(f"✅ P2Rank analysis complete for {self.protein_path.name}")
                return str(predictions_csv)
            return None
        except subprocess.CalledProcessError as e:
            print(f"❌ ERROR: P2Rank failed to run. Stderr: {e.stderr}")
            return None

def log_pockets_to_excel(predictions_csv_path: str, protein_id: str, output_excel_dir: str = "data/p2rank_analysis"):
    # ... (this function remains the same)
    if not predictions_csv_path or not Path(predictions_csv_path).exists():
        return None
    try:
        Path(output_excel_dir).mkdir(parents=True, exist_ok=True)
        output_excel_path = Path(output_excel_dir) / f"{protein_id}_pockets.xlsx"
        df = pd.read_csv(predictions_csv_path)
        df.columns = df.columns.str.strip()
        df['protein_id'] = protein_id
        top_pocket_df = df[df['rank'] == 1]
        top_5_pockets_df = df[df['rank'] <= 5]
        all_pockets_df = df
        with pd.ExcelWriter(output_excel_path, engine='openpyxl') as writer:
            top_pocket_df.to_excel(writer, sheet_name='Top_Pocket', index=False)
            top_5_pockets_df.to_excel(writer, sheet_name='Top_5_Pockets', index=False)
            all_pockets_df.to_excel(writer, sheet_name='All_Pockets', index=False)
        print(f"[INFO] 📊 Pocket analysis for '{protein_id}' saved to {output_excel_path}")
        
        # --- SIMPLIFICATION ---
        # Convert NumPy floats to standard Python floats before returning
        top_pocket_coords = top_pocket_df.iloc[0]
        return {
            'center_x': float(top_pocket_coords['center_x']),
            'center_y': float(top_pocket_coords['center_y']),
            'center_z': float(top_pocket_coords['center_z'])
        }
    except Exception as e:
        print(f"❌ ERROR: Could not parse or log pockets. Error: {e}")
        return None

if __name__ == "__main__":
    # ... (this part remains the same)
    protein_file = input("Enter path to a prepared protein PDBQT file: ").strip()
    if protein_file and Path(protein_file).exists():
        identifier = PocketIdentifier(protein_file)
        predictions_file = identifier.find_pockets()
        if predictions_file:
            protein_identifier = Path(protein_file).stem.replace('_prepared', '')
            top_pocket_coords = log_pockets_to_excel(predictions_file, protein_id=protein_identifier)
            if top_pocket_coords:
                print(f"\n[INFO] Coordinates to pass to docking engine: {top_pocket_coords}")
    else:
        print("File not found or no file provided. Exiting.")
    print("--- Test complete ---")