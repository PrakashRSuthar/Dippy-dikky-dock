# backend/modules/ligand_preparer.py

import os
import requests
from openbabel import openbabel

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
                print(f"✅ Successfully saved SDF file to {sdf_path}")
                return sdf_path
        except:
            continue
    
    raise Exception(f"Could not fetch '{identifier}' from any PubChem endpoint")

def _convert_to_pdbqt(sdf_path: str, output_dir: str):
    """Bulletproof SDF to PDBQT conversion"""
    os.makedirs(output_dir, exist_ok=True)
    
    base_name = os.path.splitext(os.path.basename(sdf_path))[0]
    output_pdbqt_path = os.path.join(output_dir, f"{base_name}.pdbqt")
    
    ob_conversion = openbabel.OBConversion()
    ob_conversion.SetInAndOutFormats("sdf", "pdbqt")
    
    mol = openbabel.OBMol()
    if not ob_conversion.ReadFile(mol, sdf_path):
        raise Exception("Could not read SDF file")
    
    # Add hydrogens
    mol.AddHydrogens()
    
    # Try energy minimization (optional - don't fail if this breaks)
    try:
        print("[INFO] Attempting energy minimization...")
        ff = openbabel.OBForceField.FindForceField("mmff94")
        if ff and ff.Setup(mol):
            ff.ConjugateGradients(500, 1.0e-4)
            ff.GetCoordinates(mol)
            print("[INFO] ✅ Energy minimization completed")
        else:
            print("[WARN] MMFF94 not available; skipping minimization (this is fine)")
    except Exception as e:
        print(f"[WARN] Minimization failed: {e} (continuing anyway)")
    
    # Add Gasteiger charges (critical for docking)
    try:
        charge_model = openbabel.OBChargeModel.FindType("gasteiger")
        if charge_model:
            charge_model.ComputeCharges(mol)
            print("[INFO] ✅ Gasteiger charges added")
        else:
            print("[WARN] Could not add Gasteiger charges")
    except Exception as e:
        print(f"[WARN] Charge calculation failed: {e}")
    
    # Write PDBQT
    if not ob_conversion.WriteFile(mol, output_pdbqt_path):
        raise Exception("Could not write PDBQT file")
    
    print(f"✅ Successfully converted and prepared {sdf_path} to {output_pdbqt_path}")
    return output_pdbqt_path

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
