# backend/modules/protein_preprocessor.py

import os
from pathlib import Path
from Bio.PDB import PDBParser, PDBIO, Select
from meeko import MoleculePreparation
from rdkit import Chem

class ChainSelect(Select):
    """Selects only a single, specified chain from the protein."""
    def __init__(self, chain_id='A'):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.get_id() == self.chain_id

class ProteinPreprocessor:
    def __init__(self, pdb_path):
        self.pdb_path = Path(pdb_path)
        if not self.pdb_path.exists():
            raise FileNotFoundError(f"Input PDB file not found: {self.pdb_path}")

        # Define output directories
        self.output_cleaned = Path("data/cleaned_proteins")
        self.output_prepared = Path("data/prepared_proteins")
        self.output_cleaned.mkdir(parents=True, exist_ok=True)
        self.output_prepared.mkdir(parents=True, exist_ok=True)
        
        # Define file paths
        self.cleaned_pdb_path = self.output_cleaned / f"{self.pdb_path.stem}_chainA.pdb"
        self.prepared_pdbqt_path = self.output_prepared / f"{self.pdb_path.stem}_prepared.pdbqt"

    def process_and_prepare(self, chain_id_to_keep='A'):
        """
        Runs the full, robust preparation pipeline.
        1. Extracts a single chain using BioPython.
        2. Prepares the single chain for docking using RDKit and Meeko.
        """
        try:
            # --- STAGE 1: Extract a single chain using BioPython ---
            print(f"[INFO] Extracting chain '{chain_id_to_keep}' from {self.pdb_path.name}...")
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("protein", self.pdb_path)
            io = PDBIO()
            io.set_structure(structure)
            io.save(str(self.cleaned_pdb_path), ChainSelect(chain_id=chain_id_to_keep))
            print(f"[INFO] ✅ Saved single chain to: {self.cleaned_pdb_path}")

            # --- STAGE 2: Prepare the single-chain PDB with RDKit and Meeko ---
            print(f"[INFO] Preparing {self.cleaned_pdb_path.name} for Vina using Meeko...")
            mol = Chem.MolFromPDBFile(str(self.cleaned_pdb_path), removeHs=False)
            if mol is None:
                raise ValueError("RDKit could not read the cleaned PDB file.")
            
            # Check for fragments again, just in case. This should now pass.
            if len(Chem.GetMolFrags(mol)) > 1:
                raise ValueError("Molecule still has multiple fragments after cleaning.")

            mol = Chem.AddHs(mol, addCoords=True)
            
            preparator = MoleculePreparation()
            preparator.prepare(mol)
            preparator.write_pdbqt_file(str(self.prepared_pdbqt_path))
            
            print(f"[INFO] ✅ Prepared protein saved to: {self.prepared_pdbqt_path}")
            return str(self.prepared_pdbqt_path)

        except Exception as e:
            print(f"❌ ERROR: Failed to prepare the protein. Reason: {e}")
            return None

if __name__ == "__main__":
    pdb_input_path = input("Enter path to the ORIGINAL downloaded PDB file (e.g., data/input_proteins/1HSG_...pdb): ").strip()
    
    if not os.path.isfile(pdb_input_path):
        print("[ERROR] ❌ File not found. Check path and try again.")
    else:
        processor = ProteinPreprocessor(pdb_input_path)
        processor.process_and_prepare()