# backend/modules/protein_preprocessor.py

import os
import subprocess
from pathlib import Path
from datetime import datetime

class ProteinPreprocessor:
    def __init__(self, pdb_path):
        self.pdb_path = Path(pdb_path)
        if not self.pdb_path.exists():
            raise FileNotFoundError(f"PDB file not found: {self.pdb_path}")
        
        # Simple output structure
        self.output_dir = Path("data/prepared_proteins")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Generate unique output filename
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.output_pdbqt = self.output_dir / f"{self.pdb_path.stem}_{timestamp}_prepared.pdbqt"
    
    def prepare_with_openbabel(self):
        """Simple OpenBabel preparation - most reliable method"""
        try:
            print(f"[INFO] 🧬 Preparing protein receptor: {self.pdb_path.name}")
            
            cmd = [
                "obabel", 
                str(self.pdb_path),
                "-O", str(self.output_pdbqt),
                "-xr",  # Rigid receptor flag
                "-p", "7.4",  # Add hydrogens at pH 7.4
                "--partialcharge", "gasteiger"  # Add charges
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Validate output
            if self.output_pdbqt.exists() and self.output_pdbqt.stat().st_size > 0:
                print(f"[INFO] ✅ Protein prepared successfully: {self.output_pdbqt}")
                return str(self.output_pdbqt)
            else:
                raise Exception("Output file is empty or missing")
                
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] ❌ OpenBabel failed: {e.stderr}")
            return None
        except FileNotFoundError:
            print("[ERROR] ❌ OpenBabel not found. Install with: conda install -c conda-forge openbabel")
            return None
        except Exception as e:
            print(f"[ERROR] ❌ Preparation failed: {e}")
            return None
    
    def prepare_manual_fallback(self):
        """Simple manual preparation without any complex dependencies"""
        try:
            print(f"[INFO] 🔧 Manual preparation fallback...")
            
            # Read original PDB and create basic PDBQT
            with open(self.pdb_path, 'r') as infile, open(self.output_pdbqt, 'w') as outfile:
                for line in infile:
                    if line.startswith(('ATOM', 'HETATM')):
                        # Extract coordinates and basic info
                        if len(line) >= 54:
                            # Simple atom type mapping
                            atom_name = line[12:16].strip()
                            element = atom_name[0] if atom_name else 'C'
                            
                            # Basic AutoDock atom type
                            autodock_type = self.get_simple_atom_type(element)
                            
                            # Rebuild line with AutoDock atom type
                            pdbqt_line = line[:78] + f"  {autodock_type:>2}\n"
                            outfile.write(pdbqt_line)
            
            if self.output_pdbqt.exists() and self.output_pdbqt.stat().st_size > 0:
                print(f"[INFO] ✅ Manual preparation completed: {self.output_pdbqt}")
                return str(self.output_pdbqt)
            else:
                raise Exception("Manual preparation failed")
                
        except Exception as e:
            print(f"[ERROR] ❌ Manual preparation failed: {e}")
            return None
    
    def get_simple_atom_type(self, element):
        """Simple atom type mapping for AutoDock"""
        mapping = {
            'C': 'C', 'N': 'N', 'O': 'O', 'S': 'S', 'P': 'P',
            'H': 'H', 'F': 'F', 'I': 'I', 'B': 'B'
        }
        return mapping.get(element.upper(), 'C')
    
    def process(self):
        """Main processing method - simple interface for integration"""
        print(f"[INFO] 🚀 Starting protein preparation...")
        
        # Try OpenBabel first
        result = self.prepare_with_openbabel()
        if result:
            return result
        
        # Fallback to manual method
        print("[INFO] 🔄 Trying fallback method...")
        result = self.prepare_manual_fallback()
        if result:
            return result
        
        print("[ERROR] ❌ All preparation methods failed")
        return None
    
    def get_info(self):
        """Return basic info for UI integration"""
        return {
            'input_file': str(self.pdb_path),
            'output_file': str(self.output_pdbqt),
            'status': 'ready' if self.output_pdbqt.exists() else 'pending'
        }

# Simple function interface for easy integration
def prepare_protein(pdb_file_path):
    """Simple function interface - easy to call from main pipeline"""
    try:
        processor = ProteinPreprocessor(pdb_file_path)
        result = processor.process()
        return result
    except Exception as e:
        print(f"[ERROR] ❌ Protein preparation error: {e}")
        return None

# Test function
def test_protein_preparation():
    """Test with a sample file"""
    test_file = input("Enter PDB file path: ").strip()
    result = prepare_protein(test_file)
    
    if result:
        print(f"[SUCCESS] 🎉 Ready for docking: {result}")
        print(f"[INFO] Output file: {result}")
    else:
        print("[FAILED] ❌ Preparation unsuccessful")

if __name__ == "__main__":
    test_protein_preparation()
