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
        
        self.output_dir = Path("data/prepared_proteins")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.output_pdbqt = self.output_dir / f"{self.pdb_path.stem}_{timestamp}_prepared.pdbqt"

    def _is_valid_pdbqt(self, pdbqt_path: Path) -> bool:
        """Validate PDBQT has proper AutoDock atom types"""
        if not pdbqt_path.exists() or pdbqt_path.stat().st_size < 100:
            return False
        
        valid_types = {
            "C", "A", "N", "NA", "NS", "OA", "OS", "F", "Mg", "MG", "P", "SA", "S", "Cl", "CL", 
            "Ca", "CA", "Mn", "MN", "Fe", "FE", "Zn", "ZN", "Br", "BR", "I", "H", "HD", "HS"
        }
        
        try:
            total_atoms = 0
            valid_atoms = 0
            
            with pdbqt_path.open() as f:
                for line in f:
                    if line.startswith(("ATOM", "HETATM")):
                        total_atoms += 1
                        if len(line) >= 79:
                            atom_type = line[77:79].strip()
                            if atom_type in valid_types:
                                valid_atoms += 1
            
            return total_atoms > 0 and (valid_atoms / total_atoms) >= 0.9
        except:
            return False

    def prepare_with_mgltools(self):
        """Try MGLTools prepare_receptor4.py (most reliable)"""
        try:
            # Check for MGLTools/ADFR
            mgl_utils = os.environ.get("MGLTOOLS_UTILS")
            if not mgl_utils:
                # Try common Windows locations
                common_paths = [
                    "C:/MGLTools/MGLToolsPckgs/AutoDockTools/Utilities24",
                    "C:/ADFR/bin",
                    "C:/Program Files/MGLTools/MGLToolsPckgs/AutoDockTools/Utilities24"
                ]
                for path in common_paths:
                    if Path(path, "prepare_receptor4.py").exists():
                        mgl_utils = path
                        break
            
            if not mgl_utils:
                return None
                
            script = Path(mgl_utils) / "prepare_receptor4.py"
            if not script.exists():
                return None
            
            # Try pythonsh first, then python
            for interpreter in ["pythonsh", "python"]:
                try:
                    cmd = [interpreter, str(script), "-r", str(self.pdb_path), "-o", str(self.output_pdbqt)]
                    subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=300)
                    
                    if self._is_valid_pdbqt(self.output_pdbqt):
                        print(f"[INFO] ✅ Protein prepared via MGLTools: {self.output_pdbqt}")
                        return str(self.output_pdbqt)
                    break
                except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
                    continue
                    
            return None
        except Exception:
            return None

    def prepare_with_openbabel_advanced(self):
        """Multiple OpenBabel strategies"""
        strategies = [
            # Strategy 1: Standard receptor preparation
            ["obabel", str(self.pdb_path), "-O", str(self.output_pdbqt), "-xr", "-p", "7.4", "--partialcharge", "gasteiger"],
            # Strategy 2: Alternative approach
            ["obabel", str(self.pdb_path), "-O", str(self.output_pdbqt), "-h", "--partialcharge", "gasteiger"],
            # Strategy 3: Simple conversion
            ["obabel", str(self.pdb_path), "-O", str(self.output_pdbqt)]
        ]
        
        for i, cmd in enumerate(strategies, 1):
            try:
                print(f"[INFO] 🧬 Trying OpenBabel strategy {i}...")
                subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=300)
                
                if self._is_valid_pdbqt(self.output_pdbqt):
                    print(f"[INFO] ✅ Protein prepared via OpenBabel (strategy {i}): {self.output_pdbqt}")
                    return str(self.output_pdbqt)
                    
            except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
                continue
                
        return None

    def prepare_with_pdb2pqr_fallback(self):
        """Use PDB2PQR + manual PDBQT conversion as fallback"""
        try:
            # First, try to clean PDB and add hydrogens
            temp_pqr = self.output_pdbqt.with_suffix('.pqr')
            
            # Simple hydrogen addition
            cmd = ["obabel", str(self.pdb_path), "-O", str(temp_pqr), "-h", "-p", "7.4"]
            subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=300)
            
            if temp_pqr.exists():
                # Convert PQR to PDBQT
                cmd2 = ["obabel", str(temp_pqr), "-O", str(self.output_pdbqt)]
                subprocess.run(cmd2, check=True, capture_output=True, text=True, timeout=300)
                
                if self._is_valid_pdbqt(self.output_pdbqt):
                    print(f"[INFO] ✅ Protein prepared via PDB2PQR fallback: {self.output_pdbqt}")
                    temp_pqr.unlink()  # cleanup
                    return str(self.output_pdbqt)
                    
        except Exception:
            pass
            
        return None

    def process(self):
        """Ultimate robust processing with multiple fallbacks"""
        print(f"[INFO] 🚀 Starting protein preparation...")
        
        # Method 1: MGLTools (most reliable for Vina)
        result = self.prepare_with_mgltools()
        if result:
            return result
        
        # Method 2: Advanced OpenBabel
        result = self.prepare_with_openbabel_advanced()
        if result:
            return result
        
        # Method 3: PDB2PQR fallback
        result = self.prepare_with_pdb2pqr_fallback()
        if result:
            return result
        
        # STOP HERE - don't proceed with invalid files
        print("[ERROR] ❌ All PDBQT preparation methods failed")
        print("[SOLUTION] 💡 Install MGLTools or ADFR Suite for reliable receptor preparation:")
        print("  1. Download MGLTools: http://mgltools.scripps.edu/")
        print("  2. Or ADFR Suite: https://ccsb.scripps.edu/adfr/downloads/")
        print("  3. Set MGLTOOLS_UTILS environment variable")
        print("[ALTERNATIVE] 🔧 Fix OpenBabel installation:")
        print("  conda uninstall openbabel")
        print("  conda install -c conda-forge openbabel")
        
        return None

def prepare_protein(pdb_file_path):
    """Bulletproof protein preparation"""
    try:
        processor = ProteinPreprocessor(pdb_file_path)
        return processor.process()
    except Exception as e:
        print(f"[ERROR] ❌ Protein preparation error: {e}")
        return None
