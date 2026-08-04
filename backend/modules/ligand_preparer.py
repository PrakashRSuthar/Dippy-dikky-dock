# backend/modules/ligand_preparer.py

import os
import requests

try:
    from openbabel import openbabel
except ImportError:
    openbabel = None

def _fetch_from_pubchem(identifier: str, output_dir: str):
    """Bulletproof PubChem fetcher"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Try multiple PubChem endpoints
    endpoints = [
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{identifier}/SDF",
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{identifier}/SDF"
    ]
    
    sdf_path = os.path.join(output_dir, f"{identifier}.sdf")
    print(f"Downloading ligand '{identifier}' from PubChem...")
    
    for url in endpoints:
        try:
            response = requests.get(url, timeout=30)
            if response.status_code == 200:
                with open(sdf_path, 'w', encoding='utf-8') as f:
                    f.write(response.text)
                print(f" Successfully saved SDF file to {sdf_path}")
                return sdf_path
        except:
            continue
    
    raise Exception(f"Could not fetch '{identifier}' from any PubChem endpoint")

def _convert_to_pdbqt(sdf_path: str, output_dir: str):
    """SDF to PDBQT conversion using RDKit/Meeko or OpenBabel fallback"""
    # 1. Try RDKit + Meeko (Gold Standard for Vina, handles complex ligands best)
    try:
        print("[INFO] Attempting RDKit + Meeko ligand preparation...")
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from meeko import MoleculePreparation
        
        os.makedirs(output_dir, exist_ok=True)
        base_name = os.path.splitext(os.path.basename(sdf_path))[0]
        output_pdbqt_path = os.path.join(output_dir, f"{base_name}.pdbqt")
        
        # Read SDF
        suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
        mols = [m for m in suppl if m is not None]
        if not mols:
            raise Exception("RDKit failed to read molecule from SDF")
        mol = mols[0]
        
        # Check if molecule has 3D coordinates.
        # If not, or if 3D coordinates are invalid/flat, generate them.
        has_3d = False
        if mol.GetNumConformers() > 0:
            conf = mol.GetConformer()
            z_coords = [conf.GetAtomPosition(i).z for i in range(mol.GetNumAtoms())]
            if any(abs(z) > 1e-3 for z in z_coords):
                has_3d = True
                
        if not has_3d:
            print("[INFO] RDKit: No valid 3D coordinates found. Generating 3D conformation via ETKDG...")
            mol = Chem.AddHs(mol)
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            embed_status = AllChem.EmbedMolecule(mol, params)
            if embed_status != 0:
                print("[WARN] RDKit ETKDGv3 embedding failed. Trying standard embedding...")
                embed_status = AllChem.EmbedMolecule(mol, randomSeed=42)
                if embed_status != 0:
                    raise Exception("RDKit failed to embed 3D coordinates for the ligand")
        else:
            mol = Chem.AddHs(mol, addCoords=True)
            
        # Energy minimization
        try:
            print("[INFO] RDKit: Attempting MMFF94 energy minimization...")
            ff_setup = AllChem.MMFFOptimizeMolecule(mol)
            if ff_setup == 0:
                print("[INFO]  RDKit MMFF94 energy minimization completed")
            else:
                print("[WARN] RDKit MMFF94 minimization returned non-zero code; proceeding anyway")
        except Exception as e:
            print(f"[WARN] RDKit MMFF94 minimization failed: {e} (proceeding anyway)")
            
        # Meeko preparation
        print("[INFO] Meeko: Preparing ligand for AutoDock Vina...")
        preparer = MoleculePreparation(
            hydrate=False,
            flexible_amides=True
        )
        setups = preparer.prepare(mol)
        from meeko import PDBQTWriterLegacy
        for setup in setups:
            pdbqt_string, is_ok, err_msg = PDBQTWriterLegacy.write_string(setup)
            if is_ok:
                with open(output_pdbqt_path, "w") as f:
                    f.write(pdbqt_string)
                break
        print(f"[OK] Ligand prepared using RDKit/Meeko: {output_pdbqt_path}")
        return output_pdbqt_path

    except Exception as e:
        print(f"[WARN] RDKit + Meeko ligand preparation failed: {e}")
        print("[INFO] Falling back to OpenBabel CLI...")

        os.makedirs(output_dir, exist_ok=True)
        base_name = os.path.splitext(os.path.basename(sdf_path))[0]
        output_pdbqt_path = os.path.join(output_dir, f"{base_name}.pdbqt")

        import subprocess
        for strategy in [
            ["obabel", sdf_path, "-O", output_pdbqt_path, "-xh", "--partialcharge", "gasteiger"],
            ["obabel", sdf_path, "-O", output_pdbqt_path, "-h", "--partialcharge", "gasteiger"],
            ["obabel", sdf_path, "-O", output_pdbqt_path]
        ]:
            try:
                r = subprocess.run(strategy, capture_output=True, text=True, timeout=120)
                if r.returncode == 0 and os.path.exists(output_pdbqt_path) and os.path.getsize(output_pdbqt_path) > 100:
                    print(f"[OK] Ligand prepared via OpenBabel CLI: {output_pdbqt_path}")
                    return output_pdbqt_path
            except:
                continue
        
        raise Exception("All ligand preparation methods failed")

def process_ligand(identifier: str, ligand_dir="data/ligands", prepared_dir="data/prepared_ligands"):
    """Ultimate robust ligand processing"""
    try:
        downloaded_sdf_path = _fetch_from_pubchem(identifier, ligand_dir)
        if downloaded_sdf_path:
            final_pdbqt_path = _convert_to_pdbqt(downloaded_sdf_path, prepared_dir)
            return final_pdbqt_path
        return None
    except Exception as e:
        print(f"[ERROR] Ligand processing failed: {e}")
        return None
