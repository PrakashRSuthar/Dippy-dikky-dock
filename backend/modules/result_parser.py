# result_parser.py
"""
Vina PDBQT result parser — Paper Section III-G
Extracts REAL binding affinities AND RMSD values from output file.
No synthetic/random data. All values are parsed from actual Vina output.
"""
import pandas as pd
from pathlib import Path
from typing import Dict, Optional

def parse_vina_results(docked_pdbqt_path: str,
                        output_csv_path: str = None,
                        pocket_info: Dict = None) -> pd.DataFrame:
    """
    Parse AutoDock Vina PDBQT output. Extracts per-pose:
      - Binding affinity (kcal/mol)
      - RMSD lower bound (Å)  — REAL value from Vina output
      - RMSD upper bound (Å)  — REAL value from Vina output
      - Pose centre of mass (X, Y, Z)
    """
    pdbqt = Path(docked_pdbqt_path)
    if not pdbqt.exists():
        print(f"[ERROR] File not found: {pdbqt}"); return pd.DataFrame()

    poses = []        # list of dicts, one per pose
    current = {}      # accumulates data for current MODEL block
    atom_buf = []     # heavy atoms for centre-of-mass

    try:
        with open(pdbqt) as fh:
            for line in fh:
                #  New pose block 
                if line.startswith("MODEL"):
                    current = {"pose": int(line.split()[-1]) if len(line.split())>1 else len(poses)+1}
                    atom_buf = []

                #  Vina result line: affinity + REAL RMSD values 
                elif line.startswith("REMARK VINA RESULT:"):
                    parts = line.split()
                    if len(parts) >= 4:
                        current["affinity"]  = float(parts[3])
                        current["rmsd_lb"]   = float(parts[4]) if len(parts) > 4 else 0.0
                        current["rmsd_ub"]   = float(parts[5]) if len(parts) > 5 else 0.0

                #  Atom coordinates (heavy atoms only) 
                elif line.startswith(("ATOM","HETATM")):
                    elem = line[76:78].strip() if len(line) > 77 else ""
                    if elem not in ("H","HD","HS"):
                        try:
                            x,y,z = float(line[30:38]),float(line[38:46]),float(line[46:54])
                            atom_buf.append((x,y,z))
                        except: pass

                #  End of pose block 
                elif line.startswith("ENDMDL"):
                    if atom_buf:
                        n = len(atom_buf)
                        current["pose_cx"] = round(sum(a[0] for a in atom_buf)/n, 3)
                        current["pose_cy"] = round(sum(a[1] for a in atom_buf)/n, 3)
                        current["pose_cz"] = round(sum(a[2] for a in atom_buf)/n, 3)
                    if "affinity" in current:
                        poses.append(dict(current))
                    current = {}; atom_buf = []

        # Handle last block if no ENDMDL at EOF
        if "affinity" in current:
            if atom_buf:
                n = len(atom_buf)
                current["pose_cx"] = round(sum(a[0] for a in atom_buf)/n, 3)
                current["pose_cy"] = round(sum(a[1] for a in atom_buf)/n, 3)
                current["pose_cz"] = round(sum(a[2] for a in atom_buf)/n, 3)
            poses.append(current)

    except Exception as e:
        print(f"[ERROR] Parsing failed: {e}"); return pd.DataFrame()

    if not poses:
        print(f"[WARN] No poses found in {pdbqt.name}"); return pd.DataFrame()

    #  Build DataFrame 
    data = {
        "Pose":                     [p["pose"] for p in poses],
        "Binding Affinity (kcal/mol)": [p["affinity"] for p in poses],
        "RMSD_lb (A)":              [p.get("rmsd_lb", 0.0) for p in poses],
        "RMSD_ub (A)":              [p.get("rmsd_ub", 0.0) for p in poses],
    }

    # Centre-of-mass columns (present when atoms were parsed)
    if any("pose_cx" in p for p in poses):
        data["Pose_CoM_X"] = [p.get("pose_cx", None) for p in poses]
        data["Pose_CoM_Y"] = [p.get("pose_cy", None) for p in poses]
        data["Pose_CoM_Z"] = [p.get("pose_cz", None) for p in poses]

    # Binding site columns from pocket_info
    if pocket_info and "primary" in pocket_info:
        pri = pocket_info["primary"]
        data["BS_Center_X"] = [pri.get("center_x",0)] * len(poses)
        data["BS_Center_Y"] = [pri.get("center_y",0)] * len(poses)
        data["BS_Center_Z"] = [pri.get("center_z",0)] * len(poses)
        data["Box_Size_X"]  = [pri.get("size_x",20)]  * len(poses)
        data["Box_Size_Y"]  = [pri.get("size_y",20)]  * len(poses)
        data["Box_Size_Z"]  = [pri.get("size_z",20)]  * len(poses)
        data["Pocket_Method"]  = [pri.get("method","?")] * len(poses)
        data["Final_Score"]    = [round(pri.get("final_score",0.0),4)] * len(poses)
        data["Consensus_Size"] = [pri.get("consensus_size",1)] * len(poses)

    df = pd.DataFrame(data)

    #  Summary output 
    best = df["Binding Affinity (kcal/mol)"].min()
    print(f"\n{'='*80}")
    print(f"DOCKING RESULTS — {pdbqt.name}")
    print(f"{'='*80}")
    print(f"  Poses parsed : {len(df)}")
    print(f"  Best affinity: {best:.2f} kcal/mol")
    print(f"  RMSD range   : {df['RMSD_lb (A)'].min():.2f} — {df['RMSD_ub (A)'].max():.2f} Å")
    if pocket_info:
        pri = pocket_info.get("primary",{})
        print(f"  Pocket method: {pri.get('method','?')}  |  consensus_size={pri.get('consensus_size',1)}")
        print(f"  Drift events : {pocket_info.get('drift_events',0)}")
    print(f"\n{df[['Pose','Binding Affinity (kcal/mol)','RMSD_lb (A)','RMSD_ub (A)']].to_string(index=False)}")
    print("="*80)

    #  Save CSV 
    if output_csv_path is None:
        output_csv_path = pdbqt.parent / f"{pdbqt.stem}_scores.csv"
    df.to_csv(output_csv_path, index=False)
    print(f"[INFO] Results saved: {output_csv_path}")
    return df

if __name__=="__main__":
    import sys
    path = sys.argv[1] if len(sys.argv)>1 else input("Docked PDBQT file: ").strip()
    if path and Path(path).exists():
        parse_vina_results(path)
    else:
        print("File not found")