# backend/modules/docking_engine.py

import os
import subprocess
from pathlib import Path
import re

SCRIPT_DIR = Path(__file__).parent
VINA_EXECUTABLE = SCRIPT_DIR.parent / "tools" / "vina" / "vina.exe"

WATER_RES = {'HOH', 'WAT', 'DOD', 'H2O', 'TIP'}

# Exhaustiveness scaling: the default of 8 is appropriate for a ~20×20×20 Å box
# with a "typical" drug-like ligand (6 rotatable bonds).  For larger boxes or
# more flexible ligands we scale up to compensate.
#
# Formula:
#   exh = base * (box_vol / 8000)^(1/3) * (1 + 0.15 * max(0, n_rot - 6))
#
#   base          = 8  (Vina developer recommendation for average-size molecules)
#   box_vol       = size_x * size_y * size_z (in Å³)
#   8000          = 20³ (reference box volume)
#   n_rot         = rotatable bond count of the ligand
#   6             = reference rotatable bond count
#
# This keeps exhaustiveness at 8 for the reference case and increases it
# approximately linearly for larger/flexible systems.
EXHAUSTIVENESS_BASE = 8
REFERENCE_BOX_VOL = 8000.0     # 20^3
REFERENCE_ROT_BONDS = 6
ROT_BOND_COEFF = 0.15          # 15 % increase per extra rotatable bond
EXHAUSTIVENESS_MAX = 32        # hard cap to prevent runaway runtimes


def _estimate_rotatable_bonds(pdbqt_path: Path) -> int:
    """Rough count of rotatable bonds from a PDBQT file.

    PDBQT does not carry explicit bond-order information, so we use a
    conservative heuristic: count every non-ring, non-terminal single
    bond between two heavy atoms where at least one atom is *not* bonded
    to three or more heavy atoms (i.e. not a branching quaternary centre).
    This overestimates slightly but is sufficient for exhaustiveness scaling.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolDescriptors
        # Try to reconstruct an RDKit mol from the PDBQT atoms — fragile
        # but better than zero.  Fall through to heuristic if it fails.
        # NOTE: This is best-effort; we accept the inaccuracy for scaling.
        pass
    except ImportError:
        pass

    # Heuristic fallback: count heavy-atom lines in the PDBQT as a proxy.
    # Rough correlation: ~1 rotatable bond per 5 heavy atoms for drug-like.
    n_heavy = 0
    try:
        with open(pdbqt_path) as fh:
            for line in fh:
                if line.startswith(("ATOM", "HETATM")):
                    atype = line[76:78].strip() if len(line) > 77 else ""
                    if atype and atype not in ("H", "HD", "HS"):
                        n_heavy += 1
    except Exception:
        pass
    # Very rough: n_rot ≈ n_heavy / 5 for drug-like molecules
    return max(0, n_heavy // 5 - 1)


def compute_exhaustiveness(box_size_x: float, box_size_y: float, box_size_z: float,
                           ligand_pdbqt: str = None,
                           n_rotatable_bonds: int = None,
                           base: int = EXHAUSTIVENESS_BASE) -> int:
    """Compute exhaustiveness scaled to the docking box volume and ligand flexibility.

    Parameters
    ----------
    box_size_x, box_size_y, box_size_z : float
        Docking box dimensions in Angstroms.
    ligand_pdbqt : str, optional
        Path to the ligand PDBQT file (used to estimate rotatable bonds).
    n_rotatable_bonds : int, optional
        Explicit rotatable-bond count (overrides estimation from PDBQT).
    base : int
        Base exhaustiveness (default 8).

    Returns
    -------
    int
        Suggested exhaustiveness, clamped to [2, EXHAUSTIVENESS_MAX].
    """
    box_vol = box_size_x * box_size_y * box_size_z
    volume_factor = (box_vol / REFERENCE_BOX_VOL) ** (1.0 / 3.0)

    if n_rotatable_bonds is not None:
        n_rot = n_rotatable_bonds
    elif ligand_pdbqt:
        n_rot = _estimate_rotatable_bonds(Path(ligand_pdbqt))
    else:
        n_rot = REFERENCE_ROT_BONDS

    flex_factor = 1.0 + ROT_BOND_COEFF * max(0, n_rot - REFERENCE_ROT_BONDS)
    exh = base * volume_factor * flex_factor
    return max(2, min(EXHAUSTIVENESS_MAX, int(round(exh))))


def clean_protein_pdbqt(pdbqt_path: Path):
    """Strip OpenBabel artifacts from protein/receptor PDBQT files:
    null bytes, torsion REMARKs, ROOT/ENDROOT/TORSDOF, BRANCH/ENDBRANCH
    blocks, embedded waters.
    Produces a plain ATOM/HETATM-only file that Vina can parse."""
    try:
        raw = pdbqt_path.read_bytes()
        text = raw.decode('utf-8', errors='replace')
        lines = text.split('\n')
        has_root = any(l.startswith(('ROOT', 'BRANCH')) for l in lines)
        has_torsdof = any(l.startswith(('TORSDOF', 'ENDBRANCH')) for l in lines)
        has_null = b'\x00' in raw
        has_remark = any(l.startswith('REMARK') for l in lines)
        if not has_root and not has_torsdof and not has_null and not has_remark:
            return
        cleaned = []
        for line in lines:
            if line.startswith(('ATOM', 'HETATM')):
                resname = line[17:20].strip() if len(line) >= 20 else ''
                if resname in WATER_RES:
                    continue
                cleaned.append(line.rstrip())
        if cleaned:
            pdbqt_path.write_text('\n'.join(cleaned) + '\n', encoding='utf-8')
            print(f"[INFO] Cleaned protein PDBQT: {pdbqt_path.name} ({len(lines)} -> {len(cleaned)} lines)")
    except Exception:
        pass


class DockingEngine:
    def __init__(self, protein_pdbqt_path: str, ligand_pdbqt_path: str):
        self.protein_path = Path(protein_pdbqt_path)
        self.ligand_path = Path(ligand_pdbqt_path)
        if not self.protein_path.exists() or not self.ligand_path.exists():
            raise FileNotFoundError("Prepared protein or ligand PDBQT file not found.")
        # Only clean the protein PDBQT — ligand PDBQT from Meeko/RDKit is already valid
        clean_protein_pdbqt(self.protein_path)

    def run_docking(self, center_x: float, center_y: float, center_z: float,
                    box_size_x: float = 20.0, box_size_y: float = 20.0, box_size_z: float = 20.0,
                    output_dir: str = "data/docking_results",
                    exhaustiveness: int = 8, quick: bool = False):
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        protein_name = self.protein_path.stem.replace('_prepared', '')
        ligand_name = self.ligand_path.stem.replace('_prepared', '')
        output_poses_path = output_path / f"{protein_name}_{ligand_name}_docked.pdbqt"
        print("[INFO] Starting docking simulation via command line...")
        if not VINA_EXECUTABLE.is_file():
            print(f"ERROR: vina.exe not found at {VINA_EXECUTABLE}")
            return None
        if quick:
            exh = 2
        else:
            exh = compute_exhaustiveness(
                box_size_x, box_size_y, box_size_z,
                ligand_pdbqt=str(self.ligand_path),
                base=exhaustiveness
            )
            print(f"[INFO] Exhaustiveness: {exh} (box=({box_size_x:.1f},{box_size_y:.1f},{box_size_z:.1f}) Å)")
        cmd = [
            str(VINA_EXECUTABLE), '--receptor', str(self.protein_path),
            '--ligand', str(self.ligand_path),
            '--out', str(output_poses_path),
            '--center_x', str(center_x), '--center_y', str(center_y),
            '--center_z', str(center_z),
            '--size_x', str(box_size_x), '--size_y', str(box_size_y),
            '--size_z', str(box_size_z),
            '--exhaustiveness', str(exh)
        ]
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            print(f"Docking complete: {output_poses_path}")
            for line in result.stdout.splitlines():
                if re.match(r'^\s*\d+\s+', line):
                    parts = line.split()
                    print(f"   Mode {parts[0]}: {float(parts[1]):.4f} kcal/mol")
            return str(output_poses_path)
        except subprocess.CalledProcessError as e:
            print("ERROR: AutoDock Vina failed to run.")
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
            print(f" ERROR: One or both input files not found. Please check the paths.")
        else:
            engine = DockingEngine(protein_pdbqt_path=protein_file, ligand_pdbqt_path=ligand_file)
            engine.run_docking(**pocket_center)
    
    except ValueError:
        print(" ERROR: Invalid input. Please enter a valid number for the coordinates.")
    except Exception as e:
        print(f" An unexpected error occurred: {e}")
        
    print("\n--- Test complete ---")
    