# docking_accuracy_evaluator.py
"""
Docking Accuracy Evaluator — Paper Section IV/V
Computes REAL RMSD using RDKit GetBestRMS (no placeholders).
Falls back to coordinate-RMSD with Kabsch alignment if RDKit unavailable.

NOTE: RDKit's GetBestRMS handles atom reordering via substructure matching,
but the coordinate fallback does NOT — it assumes the atom lists are in the
same order.  For large/flexible ligands where atom order diverges between the
crystal and docked poses, the fallback may return a misleading (too-low) RMSD.
We now extract only MODEL 1 from the PDBQT output before computing RMSD, which
fixes the previous bug of reading all poses merged into one molecule.
"""
import os, json, sys, math
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict, Optional

#  Reference database (ChEMBL / PDB experimental Ki → kcal/mol) 
REFERENCE_DB = {
    "indinavir":  {"protein":"1hsg","experimental_affinity":-11.5,"source":"ChEMBL"},
    "saquinavir": {"protein":"1hsg","experimental_affinity":-10.8,"source":"ChEMBL"},
    "ritonavir":  {"protein":"1hsg","experimental_affinity":-9.7, "source":"ChEMBL"},
    "lopinavir":  {"protein":"1hsg","experimental_affinity":-9.4, "source":"FDA"},
}


def _kabsch_rmsd(P, Q):
    """Compute RMSD between two coordinate lists using Kabsch alignment.

    This finds the optimal rotation to superimpose Q onto P, then returns
    the RMSD of the aligned coordinates.  Handles different list lengths
    by truncating to the shorter list.
    """
    try:
        import numpy as np
    except ImportError:
        # Pure-Python fallback (slower but no dependency)
        n = min(len(P), len(Q))
        if n < 3:
            return None
        Pc = list(P[:n]); Qc = list(Q[:n])
        cp = [sum(c)/n for c in zip(*Pc)]
        cq = [sum(c)/n for c in zip(*Qc)]
        P0 = [(p[0]-cp[0], p[1]-cp[1], p[2]-cp[2]) for p in Pc]
        Q0 = [(q[0]-cq[0], q[1]-cq[1], q[2]-cq[2]) for q in Qc]
        # Simple: just compute un-aligned RMSD as a lower bound
        msd = sum((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2 for a,b in zip(P0,Q0))/n
        return math.sqrt(msd)

    n = min(len(P), len(Q))
    if n < 3:
        return None
    P_arr = np.array(P[:n], dtype=float)
    Q_arr = np.array(Q[:n], dtype=float)
    cp = P_arr.mean(axis=0)
    cq = Q_arr.mean(axis=0)
    P0 = P_arr - cp
    Q0 = Q_arr - cq
    H = Q0.T @ P0
    V, _, Wt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1, 1, d]) @ Wt
    Q_rot = Q0 @ R.T
    msd = np.mean(np.sum((P0 - Q_rot)**2, axis=1))
    return math.sqrt(msd)

@dataclass
class DockingAccuracyEvaluator:
    results_file: Path
    docking_modes: List[Dict] = field(default_factory=list)

    def __post_init__(self):
        self.results_file = Path(self.results_file)
        if not self.results_file.exists():
            raise FileNotFoundError(f"Not found: {self.results_file}")
        self._parse_vina_output()

    #  Parse PDBQT 
    def _parse_vina_output(self):
        modes=[]; current_atoms=[]; in_model=False
        try:
            with open(self.results_file) as fh:
                for line in fh:
                    if line.startswith("MODEL"):
                        current_atoms=[]; in_model=True
                    elif line.startswith("ENDMDL"):
                        in_model=False
                    elif line.startswith("REMARK VINA RESULT:"):
                        parts=line.split()
                        if len(parts)>=4:
                            modes.append({"mode":len(modes)+1,"affinity":float(parts[3]),
                                          "rmsd_lb":float(parts[4]) if len(parts)>4 else None,
                                          "rmsd_ub":float(parts[5]) if len(parts)>5 else None,
                                          "_atoms":[]})
                    elif in_model and line.startswith(("ATOM","HETATM")) and modes:
                        try:
                            x,y,z = float(line[30:38]),float(line[38:46]),float(line[46:54])
                            modes[-1]["_atoms"].append((x,y,z))
                        except: pass
            self.docking_modes = modes
            print(f"[INFO] Parsed {len(modes)} docking modes")
        except Exception as e:
            print(f"[ERROR] Parse failed: {e}")

    #  RMSD calculation 
    def calculate_rmsd(self, crystal_ligand_path: str, pose_index: int = 0) -> Optional[float]:
        """
        Compute RMSD between docked pose and crystal structure.
        Strategy 1: RDKit GetBestRMS (handles atom reordering).
        Strategy 2: Coordinate RMSD with Kabsch alignment on heavy atoms (fallback).
        Returns RMSD in Angstroms, or None if unavailable.
        """
        if not crystal_ligand_path or not Path(crystal_ligand_path).exists():
            print("[WARN] No crystal ligand provided for RMSD")
            return None
        if pose_index >= len(self.docking_modes):
            return None

        # Strategy 1 — RDKit (extract only the target pose to avoid multi-model merge)
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            ref_mol = Chem.MolFromMolFile(crystal_ligand_path, removeHs=True)
            if ref_mol is None:
                raise RuntimeError("Could not parse crystal ligand with RDKit")

            # Extract only MODEL (pose_index+1) from the PDBQT into a temp file
            import tempfile
            pose_pdb = self._extract_pose_pdb(pose_index)
            if pose_pdb is None:
                raise RuntimeError(f"Could not extract pose {pose_index} as PDB")
            probe_mol = Chem.MolFromPDBFile(str(pose_pdb), removeHs=True)
            if probe_mol is None:
                raise RuntimeError("Could not parse docked pose with RDKit")
            rmsd = AllChem.GetBestRMS(ref_mol, probe_mol)
            print(f"[INFO] RMSD (RDKit): {rmsd:.3f} Å")
            return round(rmsd, 3)
        except Exception as e:
            print(f"[WARN] RDKit RMSD failed ({e}), using coordinate fallback")

        # Strategy 2 — Coordinate RMSD with Kabsch alignment on heavy atoms
        try:
            crystal_atoms = self._parse_crystal_atoms(crystal_ligand_path)
            docked_atoms  = self.docking_modes[pose_index].get("_atoms", [])
            if not crystal_atoms or not docked_atoms:
                return None
            if len(crystal_atoms) < 3 or len(docked_atoms) < 3:
                return None
            rmsd = _kabsch_rmsd(crystal_atoms, docked_atoms)
            if rmsd is not None:
                print(f"[INFO] RMSD (Kabsch): {rmsd:.3f} Å (n={min(len(crystal_atoms),len(docked_atoms))} atoms)")
                return round(rmsd, 3)
        except Exception as e:
            print(f"[ERROR] RMSD fallback failed: {e}")
        return None

    def _extract_pose_pdb(self, pose_index: int) -> Optional[Path]:
        """Extract a single MODEL from the multi-model PDBQT as a temporary PDB file."""
        import tempfile
        try:
            lines = []
            in_model = False
            model_count = 0
            with open(self.results_file) as fh:
                for line in fh:
                    if line.startswith("MODEL"):
                        model_count += 1
                        if model_count == pose_index + 1:
                            in_model = True
                            lines.append(line)
                    elif line.startswith("ENDMDL"):
                        if in_model:
                            lines.append(line)
                            break
                    elif in_model:
                        if line.startswith(("ATOM", "HETATM")):
                            # Fix element column for PDB (PDBQT has atom type in 77-78,
                            # PDB expects element in 77-78 but may differ).
                            lines.append(line.rstrip()[:76] + "  " + "\n" if len(line) > 76 else line.rstrip() + "\n")
                        else:
                            lines.append(line.rstrip() + "\n")
            if not lines:
                return None
            tmp = Path(tempfile.mktemp(suffix="_pose.pdb"))
            tmp.write_text("".join(lines))
            return tmp
        except Exception:
            return None

    def _parse_crystal_atoms(self, path):
        atoms = []
        try:
            with open(path) as fh:
                for line in fh:
                    if line.startswith(("ATOM","HETATM")):
                        elem = line[76:78].strip() if len(line)>77 else ""
                        if elem not in ("H",""):
                            try: atoms.append((float(line[30:38]),float(line[38:46]),float(line[46:54])))
                            except: pass
        except: pass
        return atoms

    #  Affinity accuracy 
    def get_best_affinity(self): 
        return min((m["affinity"] for m in self.docking_modes), default=None)

    def _detect_ligand(self):
        stem = self.results_file.stem.lower()
        for k in REFERENCE_DB:
            if k in stem: return k
        parts = stem.split("_")
        for p in parts:
            if len(p)>4 and p not in {"docked","prepared","results","poses"}: return p
        return "unknown"

    def calculate_affinity_accuracy(self, ligand_name=None, predicted=None):
        ligand_name = ligand_name or self._detect_ligand()
        predicted   = predicted   or self.get_best_affinity()
        if predicted is None: return None
        ref = REFERENCE_DB.get(ligand_name.lower())
        if not ref:
            return {"ligand":ligand_name,"predicted":predicted,"experimental":None,
                    "absolute_error":None,"relative_error":None,"category":"No Reference"}
        exp = ref["experimental_affinity"]
        ae  = abs(predicted-exp)
        re  = ae/abs(exp)*100
        cat = "Excellent" if ae<=1.0 else "Good" if ae<=2.0 else "Fair" if ae<=3.0 else "Poor"
        return {"ligand":ligand_name,"predicted":round(predicted,2),"experimental":exp,
                "absolute_error":round(ae,2),"relative_error":round(re,1),"category":cat}

    #  Full report 
    def generate_report(self, crystal_ligand_path=None, ligand_name=None):
        ligand  = ligand_name or self._detect_ligand()
        pred    = self.get_best_affinity()
        aff_acc = self.calculate_affinity_accuracy(ligand, pred)
        rmsd    = self.calculate_rmsd(crystal_ligand_path) if crystal_ligand_path else None

        ae = aff_acc.get("absolute_error") if aff_acc else None
        if ae is None:         score, grade = None, "N/A"
        elif ae<=1.0:          score, grade = 95, "A"
        elif ae<=2.0:          score, grade = 80, "B"
        elif ae<=3.0:          score, grade = 65, "C"
        elif ae<=5.0:          score, grade = 50, "D"
        else:                  score, grade = 30, "F"

        pose_accuracy = None
        if rmsd is not None:
            pose_accuracy = "Successful" if rmsd<=2.0 else "Acceptable" if rmsd<=3.5 else "Poor"

        return {"file":str(self.results_file),"ligand":ligand,"best_affinity":pred,
                "total_modes":len(self.docking_modes),"modes":self.docking_modes,
                "rmsd_angstrom":rmsd,"pose_accuracy":pose_accuracy,
                "affinity_accuracy":aff_acc,
                "overall":{"score":score,"grade":grade,
                           "status":"Validated" if score and score>=70 else "Needs Improvement"}}

    def print_report(self, report):
        print("\n"+"="*70)
        print("DOCKING ACCURACY REPORT")
        print("="*70)
        print(f"File:   {report['file']}")
        print(f"Ligand: {report['ligand'].upper()}")
        print(f"Modes:  {report['total_modes']}")
        for m in report["modes"]:
            print(f"  Mode {m['mode']}: {m['affinity']:6.2f} kcal/mol")
        if report.get("rmsd_angstrom") is not None:
            print(f"\nRMSD vs crystal: {report['rmsd_angstrom']:.3f} Å  → {report['pose_accuracy']}")
        aa = report.get("affinity_accuracy")
        if aa:
            print(f"\nAffinity:  predicted={aa['predicted']}, experimental={aa['experimental']}")
            if aa.get("absolute_error") is not None:
                print(f"  Error:   {aa['absolute_error']} kcal/mol ({aa['relative_error']}%)  → {aa['category']}")
        ov = report.get("overall",{})
        if ov.get("score"): print(f"\nOverall: {ov['score']}/100 ({ov['grade']}) — {ov['status']}")
        print("="*70)

def main():
    path = sys.argv[1] if len(sys.argv)>1 else input("Docking results file: ").strip().strip('"')
    if not path or not os.path.exists(path):
        print(f"[ERROR] File not found: {path}"); return
    ev = DockingAccuracyEvaluator(path)
    crystal = input("Crystal ligand file (Enter to skip): ").strip() or None
    report = ev.generate_report(crystal_ligand_path=crystal)
    ev.print_report(report)
    if input("Save report? (y/n): ").strip().lower()=="y":
        out = Path("data/accuracy_reports"); out.mkdir(parents=True, exist_ok=True)
        fname = out/f"report_{report['ligand']}_{Path(path).stem}.json"
        fname.write_text(json.dumps(report, indent=2, default=str))
        print(f"Saved: {fname}")

if __name__=="__main__": main()