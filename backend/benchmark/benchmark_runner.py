# backend/benchmark/benchmark_runner.py
# Run docking pipeline on Astex Diverse Set and calculate RMSD

import os, sys, json, time, math, csv, argparse, shutil, traceback, threading, concurrent.futures
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Tuple

import requests

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolDescriptors
    HAS_RDKIT = True
except ImportError:
    Chem = None
    AllChem = None
    rdMolDescriptors = None
    HAS_RDKIT = False

sys.path.insert(0, str(Path(__file__).parent.parent))
from benchmark.astex_reference import ASTEX_DIVERSE_SET
from modules.protein_preprocessor import ProteinPreprocessor
from modules.ligand_preparer import _convert_to_pdbqt
from modules.docking_engine import DockingEngine
from modules.pocket_identifier import identify_binding_site, get_pocket_analysis_config
from modules.result_parser import parse_vina_results
from modules.docking_accuracy_evaluator import DockingAccuracyEvaluator

COMMON_IONS = {"[Na+]","[K+]","[Ca+2]","[Mg+2]","[Mn+2]","[Zn+2]","[Fe+2]","[Fe+3]","[Cu+2]","[Ni+2]","[Co+2]","[Sr+2]","[Cd+2]","[Cs+]","[Rb+]","[Hg+2]","[Cl-]"}
COMMON_BUFFERS = {"CC(=O)[O-]","O=S(=O)([O-])[O-]","O=P([O-])([O-])[O-]","OCCO","OCC(O)CO","CS(C)=O","CC(C)O","CCO","N","O","[O-][SH](O)O","O=C([O-])[O-]","[PH](O)(O)=O"}


def _classify_failure(error_msg: str) -> str:
    """Classify a failure message into a semantic category.

    Categories
    ----------
    ligand_prep    — RDKit/Meeko/OpenBabel could not embed or convert the ligand.
    protein_prep   — Protein preprocessor (cleaning/PDBQT) failed.
    docking_timeout — Vina subprocess was killed or timed out.
    docking_failed  — Vina ran but returned no valid poses.
    pocket_not_found — No binding site could be identified.
    rmsd_error     — Docking succeeded but RMSD computation failed.
    download_error — PDB / UniProt fetch failed.
    unknown        — Could not determine category.
    """
    msg = error_msg.lower()
    if any(kw in msg for kw in ("sdf", "smiles", "embed", "conformer", "ligand prepar",
                                  "convert_to_pdbqt", "pdbqt", "meeko", "rdkit")):
        if "timeout" in msg:
            return "docking_timeout"
        return "ligand_prep"
    if any(kw in msg for kw in ("protein prepar", "preprocessor", "clean", "atom mapping",
                                  "missing atoms", "disulfide", "protonat")):
        return "protein_prep"
    if any(kw in msg for kw in ("timeout", "timed out", "time")):
        return "docking_timeout"
    if any(kw in msg for kw in ("no poses", "no docking results", "empty", "vina failed")):
        return "docking_failed"
    if any(kw in msg for kw in ("no binding", "pocket", "site not found")):
        return "pocket_not_found"
    if any(kw in msg for kw in ("rmsd", "rdkit", "alignment")):
        return "rmsd_error"
    if any(kw in msg for kw in ("download", "fetch", "rcsb", "http", "404", "503")):
        return "download_error"
    return "unknown"

def is_druglike_smiles(smiles: str) -> bool:
    s = smiles.strip()
    if s in COMMON_IONS or s in COMMON_BUFFERS:
        return False
    if len(s) < 8:
        return False
    heavy = sum(1 for c in s if c.isupper() and c not in "HBCNOSPFClBrI")
    return heavy >= 6

def pick_primary_ligand(smiles_list: List[str]) -> str:
    candidates = [s for s in smiles_list if is_druglike_smiles(s)]
    if not candidates:
        return max(smiles_list, key=lambda s: len(s))
    return max(candidates, key=lambda s: len(s))

_log_lock = threading.Lock()
def log(msg: str):
    ts = datetime.now().strftime("%H:%M:%S")
    with _log_lock:
        print(f"[{ts}] {msg}")

def download_pdb(pdb_id: str, dest_dir: Path) -> Path:
    dest = dest_dir / f"{pdb_id}.pdb"
    if dest.exists():
        return dest
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    log(f"Downloading {pdb_id}...")
    resp = requests.get(url, timeout=60)
    if resp.status_code != 200:
        raise RuntimeError(f"RCSB fetch failed for {pdb_id}: {resp.status_code}")
    dest.write_text(resp.text, encoding="utf-8")
    return dest

def extract_ligand_from_pdb(pdb_path: Path, ligand_code: str, output_dir: Path) -> Optional[Path]:
    seen_atoms = set()
    hetatm_lines = []
    for line in open(pdb_path):
        if line.startswith("HETATM") and line[17:20].strip() == ligand_code:
            alt = line[16:17]
            if alt and alt != " " and alt != "A":
                continue
            atom_key = line[12:16].strip()
            if atom_key in seen_atoms:
                continue
            seen_atoms.add(atom_key)
            hetatm_lines.append(line)
    if not hetatm_lines:
        log(f"  WARN: No HETATM for ligand {ligand_code}")
        return None
    lig_path = output_dir / f"{pdb_path.stem}_{ligand_code}_xtal.pdb"
    lig_path.write_text("".join(hetatm_lines) + "END\n")
    log(f"  Crystal ligand: {len(hetatm_lines)} unique atoms")
    return lig_path

def extract_protein(pdb_path: Path, ligand_code: str, output_dir: Path) -> Path:
    lines = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                lines.append(line)
            elif line.startswith("HETATM"):
                if line[17:20].strip() != ligand_code:
                    lines.append(line)
            elif line.startswith(("TER","END")):
                continue
            else:
                lines.append(line)
    prot_path = output_dir / f"{pdb_path.stem}_receptor.pdb"
    prot_path.write_text("".join(lines) + "END\n")
    return prot_path

def get_ligand_center(pdb_path: Path, ligand_pdbqt: str = None) -> Optional[Dict]:
    xs, ys, zs = [], [], []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM","HETATM")):
                try:
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                    xs.append(x); ys.append(y); zs.append(z)
                except:
                    pass
    if not xs:
        return None
    cx, cy, cz = sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs)
    from modules.pocket_identifier import _ligand_axis_max
    box_max = _ligand_axis_max(ligand_pdbqt) if ligand_pdbqt else 26.0
    sx = max(14.0, min(box_max, (max(xs)-min(xs)) + 8.0))
    sy = max(14.0, min(box_max, (max(ys)-min(ys)) + 8.0))
    sz = max(14.0, min(box_max, (max(zs)-min(zs)) + 8.0))
    return {"center_x": round(cx, 3), "center_y": round(cy, 3), "center_z": round(cz, 3),
            "size_x": round(sx, 1), "size_y": round(sy, 1), "size_z": round(sz, 1)}

def smiles_to_sdf(smiles: str, output_path: Path) -> Optional[Path]:
    if not HAS_RDKIT:
        return _smiles_to_sdf_obabel(smiles, output_path)
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            log(f"  WARN: RDKit could not parse SMILES: {smiles[:60]}")
            return None
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        status = AllChem.EmbedMolecule(mol, params)
        if status != 0:
            status = AllChem.EmbedMolecule(mol, randomSeed=42)
        try:
            AllChem.MMFFOptimizeMolecule(mol)
        except:
            pass
        w = Chem.SDWriter(str(output_path))
        w.write(mol)
        w.close()
        return output_path
    except Exception as e:
        log(f"  ERROR creating SDF: {e}")
        return None

def _smiles_to_sdf_obabel(smiles: str, output_path: Path) -> Optional[Path]:
    import subprocess, tempfile
    try:
        tmp_sdf = Path(tempfile.mktemp(suffix='.sdf'))
        cmd = ["obabel", "-:" + smiles, "-O", str(tmp_sdf), "--gen3d", "-h"]
        subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=120)
        if tmp_sdf.exists() and tmp_sdf.stat().st_size > 0:
            shutil.move(str(tmp_sdf), str(output_path))
            return output_path
    except Exception as e:
        log(f"  WARN: OpenBabel SMILES->SDF failed: {e}")
    return None

def pdbqt_to_pdb(pdbqt_path: str, output_pdb: Path):
    atype_map = {"C":"C","A":"C","N":"N","NA":"N","NS":"N","OA":"O","OS":"O","F":"F",
                 "MG":"MG","P":"P","SA":"S","S":"S","CL":"CL","CA":"CA","MN":"MN",
                 "FE":"FE","ZN":"ZN","BR":"BR","I":"I","H":"H","HD":"H","HS":"H"}
    lines = []
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith(("ATOM","HETATM")):
                atype = line[77:79].strip()
                elem = atype_map.get(atype, atype[:1])
                new_line = line[:76] + f"{elem:>2}" + "\n"
                lines.append(new_line)
            elif line.startswith("MODEL"):
                lines.append(line)
            elif line.startswith("ENDMDL"):
                lines.append(line)
                break
    if not lines:
        return None
    output_pdb.write_text("".join(lines))
    return output_pdb

def _kabsch_align(P, Q):
    cp = [sum(c)/len(c) for c in zip(*P)]
    cq = [sum(c)/len(c) for c in zip(*Q)]
    Pc = [(p[0]-cp[0], p[1]-cp[1], p[2]-cp[2]) for p in P]
    Qc = [(q[0]-cq[0], q[1]-cq[1], q[2]-cq[2]) for q in Q]
    n = min(len(Pc), len(Qc))
    Pc, Qc = Pc[:n], Qc[:n]
    import numpy as np
    A = np.dot(np.array(Qc).T, np.array(Pc))
    V, _, Wt = np.linalg.svd(A)
    d = np.sign(np.linalg.det(np.dot(V, Wt)))
    R = np.dot(V, np.dot(np.diag([1,1,d]), Wt))
    Q_rot = np.dot(np.array(Qc), R.T)
    msd = np.mean(np.sum((np.array(Pc) - Q_rot)**2, axis=1))
    return math.sqrt(msd)

def calculate_rmsd_rdkit(xtal_pdb: Path, docked_pdbqt: str) -> Optional[float]:
    docked_pdb = docked_pdbqt.replace(".pdbqt", "_pose1.pdb")
    pdbqt_to_pdb(docked_pdbqt, Path(docked_pdb))
    try:
        ref = Chem.MolFromPDBFile(str(xtal_pdb), removeHs=True)
        probe = Chem.MolFromPDBFile(docked_pdb, removeHs=True)
        if ref and probe:
            try:
                rmsd = AllChem.GetBestRMS(ref, probe)
                return round(rmsd, 3)
            except:
                pass
    except:
        pass
    try:
        xtal_coords = []
        with open(xtal_pdb) as f:
            for line in f:
                if line.startswith(("ATOM","HETATM")):
                    elem = line[76:78].strip()
                    if elem not in ("H",""):
                        try:
                            xtal_coords.append((float(line[30:38]),float(line[38:46]),float(line[46:54])))
                        except: pass
        docked_coords = []
        with open(docked_pdbqt) as f:
            in_first = False
            for line in f:
                if line.startswith("MODEL"):
                    in_first = True; docked_coords = []
                elif line.startswith("ENDMDL"):
                    break
                elif in_first and line.startswith(("ATOM","HETATM")):
                    elem = line[76:78].strip()
                    if elem not in ("H","HD","HS",""):
                        try:
                            docked_coords.append((float(line[30:38]),float(line[38:46]),float(line[46:54])))
                        except: pass
        if len(xtal_coords) < 3 or len(docked_coords) < 3:
            return None
        return round(float(_kabsch_align(xtal_coords, docked_coords)), 3)
    except:
        return None

def run_single_benchmark(pdb_id: str, ligand_code: str, smiles_list: List[str],
                         benchmark_dir: Path, skip_existing: bool = True,
                         exhaustiveness: int = 8, known_site: bool = True) -> Dict:
    entry_dir = benchmark_dir / f"{pdb_id}_{ligand_code}"
    entry_dir.mkdir(parents=True, exist_ok=True)
    result_file = entry_dir / "result.json"

    if skip_existing and result_file.exists():
        try:
            return json.loads(result_file.read_text())
        except:
            pass

    result = {"pdb": pdb_id, "ligand_code": ligand_code, "status": "failed",
              "rmsd": None, "best_affinity": None, "error": None, "duration_s": 0, "success": False,
              "failure_category": None, "ligand_heavy_atoms": None, "ligand_rotatable_bonds": None,
              "ligand_mw": None}

    t0 = time.time()
    primary_smiles = pick_primary_ligand(smiles_list) if smiles_list else None

    try:
        raw_dir = entry_dir / "raw"; prep_dir = entry_dir / "prepared"; docking_dir = entry_dir / "docking"
        for d in [raw_dir, prep_dir, docking_dir]:
            d.mkdir(parents=True, exist_ok=True)

        complex_pdb = download_pdb(pdb_id, raw_dir)
        xtal_lig = extract_ligand_from_pdb(complex_pdb, ligand_code, raw_dir)
        protein_pdb = extract_protein(complex_pdb, ligand_code, raw_dir)

        lig_sdf = smiles_to_sdf(primary_smiles, prep_dir / f"{pdb_id}_{ligand_code}.sdf")
        if lig_sdf is None:
            raise RuntimeError("Could not create SDF from SMILES")

        prepared_lig = _convert_to_pdbqt(str(lig_sdf), str(prep_dir))
        if not prepared_lig or not Path(prepared_lig).exists():
            raise RuntimeError("Ligand PDBQT preparation failed")

        processor = ProteinPreprocessor(str(protein_pdb))
        processor.output_dir = prep_dir
        processor.base_stem = f"{pdb_id}_{ligand_code}"
        processor.tmp_cleaned_pdb = prep_dir / f"{processor.base_stem}_clean.pdb"
        processor.clean_report = prep_dir / f"{processor.base_stem}_clean_report.json"
        processor.output_pdbqt = prep_dir / f"{processor.base_stem}_prepared.pdbqt"
        prepared_prot = processor.process(interactive=False, default_policy={
            "keep_waters": False, "keep_ions": True, "keep_solvents": False, "keep_cofactors": True
        })
        if not prepared_prot:
            raise RuntimeError("Protein preparation failed")

        if known_site and xtal_lig:
            center = get_ligand_center(xtal_lig, prepared_lig)
            if center is None:
                raise RuntimeError("Could not compute ligand center")
            log(f"  Using crystal ligand site: ({center['center_x']:.1f}, {center['center_y']:.1f}, {center['center_z']:.1f})")
        else:
            pocket_info = identify_binding_site(
                str(protein_pdb), prepared_prot, prepared_lig,
                **get_pocket_analysis_config("production")
            )
            if not pocket_info or "primary" not in pocket_info:
                raise RuntimeError("No binding sites found")
            center = pocket_info["primary"]
            log(f"  Using predicted site: method={center.get('method','?')} score={center.get('final_score','?'):.3f}")

        engine = DockingEngine(prepared_prot, prepared_lig)
        docked_file = engine.run_docking(
            center_x=center["center_x"], center_y=center["center_y"],
            center_z=center["center_z"],
            box_size_x=center["size_x"], box_size_y=center["size_y"],
            box_size_z=center["size_z"],
            output_dir=str(docking_dir),
            exhaustiveness=exhaustiveness
        )
        if not docked_file:
            raise RuntimeError("Docking failed")

        results_csv = docking_dir / "binding_scores.csv"
        df = parse_vina_results(docked_file, output_csv_path=results_csv)
        if df is None or df.empty:
            raise RuntimeError("No docking results parsed")

        best_affinity = float(df["Binding Affinity (kcal/mol)"].min())

        rmsd = None
        if xtal_lig and xtal_lig.exists():
            rmsd = calculate_rmsd_rdkit(xtal_lig, docked_file)

        result["status"] = "completed"
        result["rmsd"] = rmsd
        result["best_affinity"] = round(best_affinity, 2)
        result["success"] = rmsd is not None and rmsd <= 2.0

        # Attach ligand properties for post-hoc analysis
        try:
            if HAS_RDKIT:
                mol = Chem.MolFromSmiles(primary_smiles)
                if mol:
                    result["ligand_heavy_atoms"] = mol.GetNumHeavyAtoms()
                    result["ligand_rotatable_bonds"] = rdMolDescriptors.CalcNumRotatableBonds(mol)
                    result["ligand_mw"] = round(Chem.Descriptors.MolWt(mol), 1)
        except Exception:
            pass

    except Exception as e:
        result["error"] = str(e)[:200]
        result["failure_category"] = _classify_failure(str(e))
        log(f"  FAILED [{result['failure_category']}]: {e}")

    result["duration_s"] = round(time.time() - t0, 1)
    result_file.write_text(json.dumps(result, indent=2))
    return result

def print_summary_table(results: List[Dict]):
    completed = [r for r in results if r["status"] == "completed"]
    n_success = sum(1 for r in completed if r.get("success"))
    rmsd_vals = [r["rmsd"] for r in completed if r["rmsd"] is not None]

    print("\n" + "=" * 80)
    print("BENCHMARK SUMMARY - Astex Diverse Set")
    print("=" * 80)
    print(f"  Total targets:     {len(results)}")
    print(f"  Completed:         {len(completed)}")
    print(f"  Failed:            {len(results) - len(completed)}")
    if rmsd_vals:
        mean_rmsd = sum(rmsd_vals) / len(rmsd_vals)
        sorted_rmsd = sorted(rmsd_vals)
        median_rmsd = sorted_rmsd[len(sorted_rmsd) // 2]
        print(f"  Success (RMSD<=2A): {n_success}/{len(completed)} ({n_success/max(len(completed),1)*100:.1f}%)")
        print(f"  Mean RMSD:         {mean_rmsd:.3f} A")
        print(f"  Median RMSD:       {median_rmsd:.3f} A")
        print(f"  Min RMSD:          {min(rmsd_vals):.3f} A")
        print(f"  Max RMSD:          {max(rmsd_vals):.3f} A")
    print("=" * 80)

    hdr = f"{'PDB':<6} {'Lig':<4} {'RMSD(A)':<10} {'Affinity':<10} {'Time(s)':<8} {'Status':<10}"
    print(f"\n{hdr}")
    print("-" * len(hdr))
    for r in results:
        pdb = r["pdb"]; lig = r["ligand_code"]; st = r["status"]
        if st == "completed":
            rmsd_s = f'{r["rmsd"] or "?":<8}' if r.get("rmsd") else "N/A      "
            aff_s = f'{r["best_affinity"]:<8}'
        else:
            rmsd_s = "N/A      "; aff_s = "N/A      "
        dur = r["duration_s"]
        marker = "" if r.get("success") else ("" if st == "completed" else "!")
        print(f'{pdb:<6}{lig:<5} {rmsd_s:<10} {aff_s:<10} {dur:<8} {st:<10} {marker}')

def generate_report(results: List[Dict], output_dir: Path):
    report_path = output_dir / "benchmark_report.json"
    summary_csv = output_dir / "benchmark_results.csv"

    completed = [r for r in results if r["status"] == "completed"]
    rmsd_vals = [r["rmsd"] for r in completed if r["rmsd"] is not None]
    n_success = sum(1 for r in completed if r.get("success"))

    report = {
        "benchmark": "Astex Diverse Set",
        "date": datetime.now().isoformat(),
        "total_targets": len(results),
        "completed": len(completed),
        "failed": len(results) - len(completed),
        "n_success_rmsd2": n_success,
        "success_rate": round(n_success / max(len(completed), 1) * 100, 1),
        "mean_rmsd": round(sum(rmsd_vals) / len(rmsd_vals), 3) if rmsd_vals else None,
        "median_rmsd": round(sorted(rmsd_vals)[len(rmsd_vals) // 2], 3) if rmsd_vals else None,
        "min_rmsd": round(min(rmsd_vals), 3) if rmsd_vals else None,
        "max_rmsd": round(max(rmsd_vals), 3) if rmsd_vals else None,
        "results": results
    }
    report_path.write_text(json.dumps(report, indent=2))
    log(f"Report: {report_path}")

    with open(summary_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["PDB","Ligand","Status","RMSD","BestAffinity","Success_RMSD2","Duration_s","Error"])
        for r in results:
            w.writerow([r["pdb"],r["ligand_code"],r["status"],r.get("rmsd",""),
                       r.get("best_affinity",""),r.get("success",""),r["duration_s"],r.get("error","")])
    log(f"CSV: {summary_csv}")

def main():
    parser = argparse.ArgumentParser(description="Astex Diverse Set Benchmark")
    parser.add_argument("--targets", "-t", nargs="+", help="Specific PDB IDs")
    parser.add_argument("--limit", "-n", type=int, help="Limit to first N targets")
    parser.add_argument("--output", "-o", default=f"temp_runs/astex_benchmark_{datetime.now().strftime('%Y%m%d_%H%M%S')}", help="Output directory")
    parser.add_argument("--no-skip", action="store_true", help="Re-run all targets")
    parser.add_argument("--exhaustiveness", type=int, default=8, help="Vina exhaustiveness")
    parser.add_argument("--blind", action="store_true", help="Blind docking (use predicted pockets)")
    parser.add_argument("--workers", "-w", type=int, default=1, help="Parallel workers (default: 1)")
    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    results_dir = output_dir / "results"
    results_dir.mkdir(exist_ok=True)

    targets = ASTEX_DIVERSE_SET
    if args.targets:
        targets = [t for t in targets if t["pdb"] in args.targets]
    if args.limit:
        targets = targets[:args.limit]

    log(f"Benchmark: {len(targets)} targets -> {output_dir}")
    log(f"Mode: {'blind docking' if args.blind else 'self-docking (known site)'}")

    results = [None] * len(targets)
    n_workers = min(args.workers or 1, len(targets))

    def run_one(idx_entry):
        idx, entry = idx_entry
        pdb, lig_code = entry["pdb"], entry["ligand_code"]
        log(f"[{idx+1}/{len(targets)}] {pdb}_{lig_code}")
        res = run_single_benchmark(pdb, lig_code, entry["smiles"], results_dir,
                                   skip_existing=not args.no_skip,
                                   exhaustiveness=args.exhaustiveness,
                                   known_site=not args.blind)
        status = res["status"]
        rmsd_s = f" RMSD={res['rmsd']:.3f}" if res.get("rmsd") else ""
        log(f"  [{idx+1}/{len(targets)}] -> {status}{rmsd_s} ({res['duration_s']}s)")
        return idx, res

    if n_workers > 1:
        log(f"Using {n_workers} parallel workers")
        with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as ex:
            for idx, res in ex.map(run_one, enumerate(targets)):
                results[idx] = res
    else:
        for idx, res in map(run_one, enumerate(targets)):
            results[idx] = res

    print_summary_table(results)
    generate_report(results, output_dir)
    log("Done!")

if __name__ == "__main__":
    main()
