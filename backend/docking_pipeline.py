# backend/docking_pipeline.py

import os
from pathlib import Path
from datetime import datetime
import json

# Import your existing modules
from modules.protein_fetcher import ProteinFetcher
from modules.protein_preprocessor import prepare_protein
from modules.ligand_preparer import process_ligand, _convert_to_pdbqt
from modules.pocket_identifier import identify_binding_site
from modules.docking_engine import DockingEngine
from modules.result_parser import parse_vina_results
from modules.docking_accuracy_evaluator import evaluate_docking_accuracy_from_file


def json_dump(obj) -> str:
    return json.dumps(obj, indent=2, default=str)


def resolve_protein_input(user_input: str, raw_dir: Path) -> str:
    """
    Resolve protein input:
    - If 'user_input' is an existing file path, copy it under raw_dir and return the path.
    - Else treat as an ID: try PDB ID via RCSB, fallback to AlphaFold (UniProt ID).
    Returns a path to a raw PDB file.
    """
    raw_dir.mkdir(parents=True, exist_ok=True)
    p = Path(user_input)

    if p.exists() and p.is_file():
        dst = raw_dir / p.name
        if dst.resolve() != p.resolve():
            dst.write_bytes(p.read_bytes())
        print(f"[INFO] Protein file resolved: {dst}")
        return str(dst)

    # Not a file path: fetch
    fetcher = ProteinFetcher(download_dir=str(raw_dir))
    try:
        if len(user_input) == 4 and user_input.isalnum():
            # PDB ID heuristic
            pdb_path = fetcher.fetch_from_pdb(user_input)
            return pdb_path
        else:
            # Try AlphaFold with UniProt-like ID
            pdb_path = fetcher.fetch_from_alphafold(user_input)
            return pdb_path
    except Exception as e:
        raise FileNotFoundError(f"Failed to resolve protein input '{user_input}': {e}")


def resolve_ligand_input(user_input: str, raw_dir: Path, prepared_dir: Path) -> str:
    """
    Resolve ligand input:
    - If 'user_input' is an existing file, convert to PDBQT if needed.
    - Else treat as PubChem name/CID: fetch + prepare via process_ligand().
    Returns path to prepared ligand PDBQT.
    """
    raw_dir.mkdir(parents=True, exist_ok=True)
    prepared_dir.mkdir(parents=True, exist_ok=True)
    p = Path(user_input)

    if p.exists() and p.is_file():
        dst = raw_dir / p.name
        if dst.resolve() != p.resolve():
            dst.write_bytes(p.read_bytes())

        ext = p.suffix.lower()
        if ext == ".pdbqt":
            print(f"[INFO] Ligand already in PDBQT format: {dst}")
            return str(dst)
        elif ext in {".sdf", ".mol2", ".pdb"}:
            pdbqt_path = _convert_to_pdbqt(str(dst), str(prepared_dir))
            if not pdbqt_path:
                raise RuntimeError(f"Ligand conversion to PDBQT failed for {dst}")
            return pdbqt_path
        else:
            raise ValueError(f"Unsupported ligand format: {ext}")

    # Not a file path: fetch from PubChem
    pdbqt_path = process_ligand(user_input, ligand_dir=str(raw_dir), prepared_dir=str(prepared_dir))
    if not pdbqt_path:
        raise FileNotFoundError(f"Failed to fetch/prepare ligand '{user_input}' from PubChem.")
    return pdbqt_path


def docking_pipeline(protein_input: str, ligand_input: str):
    """
    Full docking pipeline (IDs or file paths):
    1. Resolve and prepare protein (PDBQT)
    2. Resolve and prepare ligand (PDBQT)
    3. Identify binding pocket (P2Rank + ligand-aware)
    4. Run docking (AutoDock Vina)
    5. Parse results & optional accuracy evaluation
    """

    # --- Step 0: Setup output folder ---
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = Path(f"runs/{run_id}")
    raw_dir = run_dir / "raw"
    prepared_dir = run_dir / "prepared"
    pocket_dir = run_dir / "pocket"
    docking_dir = run_dir / "docking"
    for d in (raw_dir, prepared_dir, pocket_dir, docking_dir):
        d.mkdir(parents=True, exist_ok=True)

    # --- Resolve Protein (path or ID) ---
    raw_protein_path = resolve_protein_input(protein_input, raw_dir)

    # --- Prepare receptor PDBQT ---
    prepared_protein_pdbqt = prepare_protein(raw_protein_path)
    if not prepared_protein_pdbqt:
        raise RuntimeError("Protein preparation failed")

    # --- Resolve Ligand (path or PubChem name/CID) ---
    prepared_ligand_pdbqt = resolve_ligand_input(ligand_input, raw_dir, prepared_dir)

    # --- Identify pocket ---
    # Prefer original raw protein (PDB) for P2Rank; provide prepared receptor & ligand PDBQT for ligand-aware sizing
    pocket_info = identify_binding_site(
        raw_protein_path,               # protein_pdb (optional but preferred)
        prepared_protein_pdbqt,         # protein_pdbqt (required for docking)
        prepared_ligand_pdbqt,          # ligand_pdbqt (required and used for ligand-aware sizing)
        use_validated=True,
        return_n=1
    )
    if not pocket_info:
        raise RuntimeError("Pocket detection failed")

    center = pocket_info["primary"]
    # Persist a copy of pocket info into run_dir for convenience
    (pocket_dir / "pocket_selection.json").write_text(json_dump(pocket_info))

    # --- Docking ---
    engine = DockingEngine(
        protein_pdbqt_path=prepared_protein_pdbqt,
        ligand_pdbqt_path=prepared_ligand_pdbqt
    )
    docked_file = engine.run_docking(
        center_x=center["center_x"], center_y=center["center_y"], center_z=center["center_z"],
        box_size_x=center["size_x"], box_size_y=center["size_y"], box_size_z=center["size_z"],
        output_dir=str(docking_dir)
    )
    if not docked_file:
        raise RuntimeError("Docking failed.")

    # --- Parse results ---
    results_csv = docking_dir / "binding_scores.csv"
    df = parse_vina_results(docked_file, output_csv_path=results_csv)
    if df is None or df.empty:
        raise RuntimeError("No docking results parsed.")

    # --- Optional accuracy evaluation ---
    accuracy_report = evaluate_docking_accuracy_from_file(docked_file)

    # --- Summarize and save ---
    best_affinity = float(df["Binding Affinity (kcal/mol)"].min())
    summary = {
        "run_dir": str(run_dir.resolve()),
        "raw_protein": raw_protein_path,
        "prepared_protein_pdbqt": prepared_protein_pdbqt,
        "prepared_ligand_pdbqt": prepared_ligand_pdbqt,
        "pocket_method": center["method"],
        "pocket_center": [center["center_x"], center["center_y"], center["center_z"]],
        "pocket_size": [center["size_x"], center["size_y"], center["size_z"]],
        "best_affinity": best_affinity,
        "docked_file": docked_file,
        "results_csv": str(results_csv),
        "accuracy_report": accuracy_report
    }
    (run_dir / "results.json").write_text(json_dump(summary))

    print("\n=== FINAL SUMMARY ===")
    print(f"Run folder: {summary['run_dir']}")
    print(f"Best binding affinity: {best_affinity} kcal/mol")
    if accuracy_report:
        print("Accuracy evaluation complete (see detailed output above).")

    return summary


if __name__ == "__main__":
    print("=== Automated Docking Pipeline (IDs or file paths) ===")
    prot = input("Enter protein file path or ID (e.g., 4lde or path/to/file.pdb): ").strip().strip('"')
    lig = input("Enter ligand file path or PubChem name/CID (e.g., darunavir or path/to/file.sdf): ").strip().strip('"')

    results = docking_pipeline(prot, lig)

    print("\n[PIPELINE COMPLETED]")
    print(f"Run directory: {results['run_dir']}")
    print(f"Best affinity: {results['best_affinity']} kcal/mol")
