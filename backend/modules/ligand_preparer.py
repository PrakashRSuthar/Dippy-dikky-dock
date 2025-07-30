# backend/modules/ligand_preparer.py

import os
import requests
from openbabel import openbabel

def _fetch_from_pubchem(identifier: str, output_dir: str):
    """(Internal function) Fetches SDF file from PubChem."""
    os.makedirs(output_dir, exist_ok=True)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{identifier}/SDF"
    sdf_path = os.path.join(output_dir, f"{identifier}.sdf")

    print(f"Downloading ligand '{identifier}' from PubChem...")
    try:
        response = requests.get(url, timeout=15)
        if response.status_code == 404:
            raise requests.exceptions.HTTPError(f"Compound '{identifier}' not found on PubChem.")
        response.raise_for_status()
        with open(sdf_path, 'w', encoding='utf-8') as f:
            f.write(response.text)
        print(f"✅ Successfully saved SDF file to {sdf_path}")
        return sdf_path
    except requests.exceptions.RequestException as e:
        print(f"❌ Error downloading '{identifier}': {e}")
        return None

def _convert_to_pdbqt(sdf_path: str, output_dir: str):
    """(Internal function) Converts SDF to PDBQT, adding charges and performing energy minimization."""
    os.makedirs(output_dir, exist_ok=True)
    base_name = os.path.splitext(os.path.basename(sdf_path))[0]
    output_pdbqt_path = os.path.join(output_dir, f"{base_name}.pdbqt")

    ob_conversion = openbabel.OBConversion()
    ob_conversion.SetInAndOutFormats("sdf", "pdbqt")
    mol = openbabel.OBMol()
    if not ob_conversion.ReadFile(mol, sdf_path):
        return None

    # --- ADDED: Energy Minimization Step ---
    print("[INFO] Performing energy minimization...")
    mol.AddHydrogens() # Hydrogens are needed for forcefield calculations
    ff = openbabel.OBForceField.FindForceField("mmff94") # MMFF94 is a good general-purpose forcefield
    if ff:
        ff.Setup(mol)
        ff.ConjugateGradients(500, 1.0e-4) # Run for 500 steps
        ff.GetCoordinates(mol)
    else:
        print("[WARN] Could not find MMFF94 forcefield. Skipping energy minimization.")
    
    # --- Add Gasteiger charges (standard for Vina) ---
    charge_model = openbabel.OBChargeModel.FindType("gasteiger")
    if charge_model:
        charge_model.ComputeCharges(mol)
    
    # --- Write the final PDBQT file ---
    if not ob_conversion.WriteFile(mol, output_pdbqt_path):
        return None

    print(f"✅ Successfully converted and prepared {sdf_path} to {output_pdbqt_path}")
    return output_pdbqt_path

def process_ligand(identifier: str, ligand_dir="data/ligands", prepared_dir="data/prepared_ligands"):
    """Main function to fetch a ligand from PubChem and prepare it for docking."""
    downloaded_sdf_path = _fetch_from_pubchem(identifier, ligand_dir)
    if downloaded_sdf_path:
        final_pdbqt_path = _convert_to_pdbqt(downloaded_sdf_path, prepared_dir)
        return final_pdbqt_path
    return None

# This block allows you to run this file directly as a script for testing
if __name__ == "__main__":
    # --- UPDATED: Now asks for user input ---
    ligand_name = input("Enter ligand name or PubChem CID (e.g., caffeine, aspirin): ").strip()
    
    if not ligand_name:
        print("No ligand name provided. Exiting.")
    else:
        print(f"\n--- Running full fetch and preparation for '{ligand_name}' ---")
        final_path = process_ligand(ligand_name)
        if final_path:
            print(f"\n✅ Pipeline successful. Final file at: {final_path}")
        else:
            print(f"\n❌ Pipeline failed for '{ligand_name}'.")
    print("\n--- Test complete ---")