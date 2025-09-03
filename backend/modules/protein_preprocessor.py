# backend/modules/protein_preprocessor.py
# Dual-mode protein preprocessor: CLI and API support
# - Cleans PDB with built-in cleaner (summary + optional prompts)
# - Converts to PDBQT using MGLTools -> OpenBabel -> PDB2PQR fallback
# - Automatic fallback to original PDB if cleaned conversion fails
# - Emits JSON cleaning report for UI
# - Works both as CLI tool and imported module

import os
import json
import sys
import argparse
import subprocess
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Optional

# Removable residue classes
REM_WATERS = {"HOH","WAT","DOD","H2O","W","TIP"}
REM_SOLVENTS = {
    "MPD","PEG","PG4","PG","BME","DMS","GOL","IPA","EOH","ACT","ACE",
    "FMT","ACY","TRS","TES","HEP","MES","MOPS","PO4","SO4",
    "NAG","NDG","MAN","GLC","GAL","BMA","FUC"
}
REM_IONS = {"NA","K","CL","MG","MN","CA","ZN","CU","FE","FE2","CD","NI","CO","SR","YB","CS","RB"}
COFACTORS = {"HEM","HEA","HEC","FAD","FMN","NAD","NAP","NDP","SAM","SAH","PLP","ATP","ADP","AMP","GTP","GDP","GNP","COA"}

class ProteinPreprocessor:
    def __init__(self, pdb_path: str):
        self.pdb_path = Path(pdb_path)
        if not self.pdb_path.exists():
            raise FileNotFoundError(f"PDB file not found: {self.pdb_path}")

        self.output_dir = Path("data/prepared_proteins")
        self.output_dir.mkdir(parents=True, exist_ok=True)

        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.base_stem = f"{self.pdb_path.stem}_{ts}"
        self.tmp_cleaned_pdb = self.output_dir / f"{self.base_stem}_clean.pdb"
        self.clean_report = self.output_dir / f"{self.base_stem}_clean_report.json"
        self.output_pdbqt = self.output_dir / f"{self.base_stem}_prepared.pdbqt"

    def _parse_atom_line(self, line: str) -> Tuple[str,str,str,str,str]:
        """Parse PDB line using fixed-width columns."""
        if len(line) < 78:
            return "", "", "", "", ""
        rec = line[:6].strip()
        resn = line[17:20].strip()
        chain = line[21:22].strip()
        alt = line[16:17].strip()
        elem = line[76:78].strip().upper()
        return rec, resn, chain, alt, elem

    def _is_valid_pdbqt(self, pdbqt_path: Path) -> bool:
        """Validate PDBQT has proper AutoDock/Vina atom types."""
        if not pdbqt_path.exists() or pdbqt_path.stat().st_size < 100:
            return False
        valid_types = {
            "C","A","N","NA","NS","OA","OS","F","Mg","MG","P","SA","S","Cl","CL",
            "Ca","CA","Mn","MN","Fe","FE","Zn","ZN","Br","BR","I","H","HD","HS"
        }
        try:
            total, valid = 0, 0
            with pdbqt_path.open() as f:
                for line in f:
                    if line.startswith(("ATOM","HETATM")):
                        total += 1
                        if len(line) >= 79 and line[77:79].strip() in valid_types:
                            valid += 1
            return total > 0 and (valid / total) >= 0.9
        except:
            return False

    def summarize_pdb(self, pdb_file: Path) -> Dict[str,int]:
        """Scan PDB and count removable classes."""
        counts = {"waters": 0, "solvents": 0, "ions": 0, "cofactors": 0, "hetatm_total": 0}
        with pdb_file.open() as f:
            for line in f:
                if not line.startswith(("ATOM","HETATM")):
                    continue
                rec, resn, chain, alt, elem = self._parse_atom_line(line)
                if rec != "HETATM":
                    continue
                counts["hetatm_total"] += 1
                if resn in REM_WATERS:
                    counts["waters"] += 1
                elif resn in REM_IONS or elem in REM_IONS:
                    counts["ions"] += 1
                elif resn in COFACTORS:
                    counts["cofactors"] += 1
                elif resn in REM_SOLVENTS:
                    counts["solvents"] += 1
        return counts

    def clean_pdb(self, pdb_file: Path, policy: Dict[str,bool]) -> Path:
        """Clean PDB according to policy and write report."""
        keep_w = policy.get("keep_waters", False)
        keep_i = policy.get("keep_ions", False)
        keep_s = policy.get("keep_solvents", False)
        keep_c = policy.get("keep_cofactors", False)

        seen_alt = set()
        kept_lines = []
        report = {
            "policy": policy,
            "removed": {"waters": 0, "solvents": 0, "ions": 0, "cofactors": 0, "altloc": 0},
            "kept": {"protein": 0, "hetatm": 0}
        }

        with pdb_file.open() as f:
            for line in f:
                if line.startswith(("ATOM","HETATM")):
                    rec, resn, chain, alt, elem = self._parse_atom_line(line)
                    
                    # Handle altLoc duplicates
                    atom_key = (line[6:11], line[12:16], line[17:26])
                    if alt and atom_key in seen_alt:
                        report["removed"]["altloc"] += 1
                        continue
                    if alt:
                        seen_alt.add(atom_key)

                    is_water = resn in REM_WATERS
                    is_ion = resn in REM_IONS or elem in REM_IONS
                    is_solvent = resn in REM_SOLVENTS
                    is_cofactor = resn in COFACTORS

                    if rec == "HETATM":
                        if is_water and not keep_w:
                            report["removed"]["waters"] += 1
                            continue
                        if is_ion and not keep_i:
                            report["removed"]["ions"] += 1
                            continue
                        if is_solvent and not keep_s and not is_cofactor:
                            report["removed"]["solvents"] += 1
                            continue
                        if is_cofactor and not keep_c:
                            report["removed"]["cofactors"] += 1
                            continue
                        report["kept"]["hetatm"] += 1
                    else:
                        report["kept"]["protein"] += 1

                    kept_lines.append(line)
                elif line.startswith(("MODEL","TER","END")):
                    continue
                else:
                    kept_lines.append(line)

        self.tmp_cleaned_pdb.write_text("".join(kept_lines) + "\nEND\n")
        
        # Write detailed report
        full_report = {
            "summary": self.summarize_pdb(pdb_file),
            "final_report": report
        }
        self.clean_report.write_text(json.dumps(full_report, indent=2))
        
        print(f"[INFO] 🧼 Cleaned PDB: {self.tmp_cleaned_pdb}")
        print(f"[INFO] 🧾 Clean report: {self.clean_report}")
        return self.tmp_cleaned_pdb

    def prepare_with_mgltools(self, receptor_pdb: Path) -> Optional[str]:
        """Try MGLTools/ADFR prepare_receptor."""
        try:
            mgl_utils = os.environ.get("MGLTOOLS_UTILS")
            if not mgl_utils:
                # Try common Windows ADFR locations
                for path in [
                    "C:/Program Files/ADFR/bin",
                    "C:/Program Files (x86)/ADFR/bin", 
                    "C:/ADFR/bin"
                ]:
                    if Path(path).exists():
                        mgl_utils = path
                        break
            
            if not mgl_utils:
                return None
                
            # Try prepare_receptor (ADFR) or prepare_receptor4.py (MGLTools)
            for script_name in ["prepare_receptor", "prepare_receptor4.py"]:
                script = Path(mgl_utils) / script_name
                if script.exists():
                    for interpreter in ["python", "pythonsh"]:
                        try:
                            cmd = [interpreter, str(script), "-r", str(receptor_pdb), "-o", str(self.output_pdbqt)]
                            subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=300)
                            if self._is_valid_pdbqt(self.output_pdbqt):
                                print(f"[INFO] ✅ Prepared via {script_name}: {self.output_pdbqt}")
                                return str(self.output_pdbqt)
                        except Exception as e:
                            print(f"[WARN] {script_name} failed with {interpreter}: {str(e)}")
            return None
        except Exception:
            return None

    def prepare_with_openbabel(self, receptor_pdb: Path) -> Optional[str]:
        """Try multiple OpenBabel strategies."""
        strategies = [
            ["obabel", str(receptor_pdb), "-O", str(self.output_pdbqt), "-xr", "-p", "7.4", "--partialcharge", "gasteiger"],
            ["obabel", str(receptor_pdb), "-O", str(self.output_pdbqt), "-h", "-p", "7.4", "--partialcharge", "gasteiger"],
            ["obabel", str(receptor_pdb), "-O", str(self.output_pdbqt)]
        ]
        
        for i, cmd in enumerate(strategies, 1):
            try:
                print(f"[INFO] 🧬 Trying OpenBabel strategy {i}...")
                subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=300)
                if self._is_valid_pdbqt(self.output_pdbqt):
                    print(f"[INFO] ✅ Prepared via OpenBabel strategy {i}: {self.output_pdbqt}")
                    return str(self.output_pdbqt)
            except Exception as e:
                print(f"[WARN] OpenBabel strategy {i} failed: {str(e)}")
        return None

    def prepare_with_pdb2pqr_fallback(self, receptor_pdb: Path) -> Optional[str]:
        """Fallback: PDB -> PQR -> PDBQT."""
        try:
            temp_pqr = self.output_pdbqt.with_suffix('.pqr')
            
            # Step 1: PDB to PQR with hydrogens
            cmd1 = ["obabel", str(receptor_pdb), "-O", str(temp_pqr), "-h", "-p", "7.4"]
            subprocess.run(cmd1, check=True, capture_output=True, text=True, timeout=300)
            
            if temp_pqr.exists():
                # Step 2: PQR to PDBQT
                cmd2 = ["obabel", str(temp_pqr), "-O", str(self.output_pdbqt)]
                subprocess.run(cmd2, check=True, capture_output=True, text=True, timeout=300)
                
                if self._is_valid_pdbqt(self.output_pdbqt):
                    print(f"[INFO] ✅ Prepared via PDB2PQR fallback: {self.output_pdbqt}")
                    temp_pqr.unlink(missing_ok=True)
                    return str(self.output_pdbqt)
        except Exception as e:
            print(f"[WARN] PDB2PQR fallback failed: {str(e)}")
        return None

    def process(self, interactive: bool = False, default_policy: Optional[Dict[str,bool]] = None) -> Optional[str]:
        """Main processing workflow with detailed analysis."""
        print("[INFO] 🚀 Starting protein preparation...")
        
        # 1) Detailed analysis
        summary = self.summarize_pdb(self.pdb_path)
        print(f"\n📊 DETAILED ANALYSIS of {self.pdb_path.name}:")
        print(f"   💧 Waters: {summary['waters']} molecules")
        print(f"   🧪 Solvents/Buffers: {summary['solvents']} molecules") 
        print(f"   ⚛️  Ions/Metals: {summary['ions']} atoms")
        print(f"   🔋 Cofactors: {summary['cofactors']} molecules")
        print(f"   📝 Total HETATM: {summary['hetatm_total']} records")
        
        # 2) Determine cleaning policy
        policy = default_policy or {
            "keep_waters": False,
            "keep_ions": False, 
            "keep_solvents": False,
            "keep_cofactors": False
        }
        
        # 3) Show what will happen
        if not interactive:
            print(f"\n🎯 CLEANING POLICY (Non-Interactive Mode):")
            print(f"   💧 Waters: {'KEEP' if policy['keep_waters'] else 'REMOVE'} ({summary['waters']} found)")
            print(f"   🧪 Solvents: {'KEEP' if policy['keep_solvents'] else 'REMOVE'} ({summary['solvents']} found)")
            print(f"   ⚛️  Ions: {'KEEP' if policy['keep_ions'] else 'REMOVE'} ({summary['ions']} found)")
            print(f"   🔋 Cofactors: {'KEEP' if policy['keep_cofactors'] else 'REMOVE'} ({summary['cofactors']} found)")
            print(f"   💡 Use --interactive flag to change these settings")
        else:
            print(f"\n🤔 Interactive mode: You'll be asked what to keep/remove...")
            policy = self._prompt_policy_interactive(summary, policy)
        
        # 4) Clean PDB
        cleaned = self.clean_pdb(self.pdb_path, policy)
        
        # 5) Try conversion on cleaned PDB
        for method in [self.prepare_with_mgltools, self.prepare_with_openbabel, self.prepare_with_pdb2pqr_fallback]:
            result = method(cleaned)
            if result:
                return result
        
        # 6) Fallback: try original PDB
        print("[WARN] Cleaned PDB conversion failed; retrying with original PDB...")
        for method in [self.prepare_with_mgltools, self.prepare_with_openbabel, self.prepare_with_pdb2pqr_fallback]:
            result = method(self.pdb_path)
            if result:
                return result
        
        print("[ERROR] ❌ All PDBQT preparation methods failed")
        print("[SOLUTION] 💡 Install ADFR Suite or ensure OpenBabel is properly configured")
        return None

    def _prompt_policy_interactive(self, counts: Dict[str,int], defaults: Dict[str,bool]) -> Dict[str,bool]:
        """Interactive prompts for cleaning policy."""
        def ask(key: str, display: str, default: bool) -> bool:
            count = counts.get(key, 0)
            default_char = "Y" if default else "N"
            other_char = "n" if default else "y"
            prompt = f"Keep {display}? ({count} found) [{default_char}/{other_char}]: "
            ans = input(prompt).strip().lower()
            if ans == 'y': 
                return True
            elif ans == 'n': 
                return False
            return default
        
        return {
            "keep_waters": ask("waters", "waters", defaults.get("keep_waters", False)),
            "keep_solvents": ask("solvents", "solvents", defaults.get("keep_solvents", False)), 
            "keep_ions": ask("ions", "ions/metals", defaults.get("keep_ions", False)),
            "keep_cofactors": ask("cofactors", "cofactors", defaults.get("keep_cofactors", False))
        }

def prepare_protein(pdb_file_path: str, interactive: bool = False, default_policy: Optional[Dict[str,bool]] = None) -> Optional[str]:
    """
    API function for programmatic use.
    
    Args:
        pdb_file_path: Path to input PDB file
        interactive: If True, prompt user about keeping removed items
        default_policy: Dict specifying what to keep (for non-interactive use)
    
    Returns:
        Path to prepared PDBQT file or None if failed
    """
    try:
        processor = ProteinPreprocessor(pdb_file_path)
        return processor.process(interactive=interactive, default_policy=default_policy)
    except Exception as e:
        print(f"[ERROR] ❌ Protein preparation error: {e}")
        return None

def main():
    """CLI entry point with intelligent file handling."""
    parser = argparse.ArgumentParser(
        description="Protein preprocessor for molecular docking",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Interactive mode with prompts
  python protein_preprocessor.py input.pdb --interactive
  
  # Non-interactive with specific policy
  python protein_preprocessor.py input.pdb --keep-cofactors
  
  # Keep multiple classes
  python protein_preprocessor.py input.pdb --keep-waters --keep-ions
        """
    )
    
    # Make PDB file optional initially
    parser.add_argument(
        "pdb_file", 
        nargs='?',  # Makes it optional
        help="Input PDB file path"
    )
    
    parser.add_argument(
        "--interactive", "-i",
        action="store_true",
        help="Enable interactive prompts for cleaning decisions"
    )
    
    parser.add_argument(
        "--keep-waters",
        action="store_true", 
        help="Keep crystallographic waters"
    )
    
    parser.add_argument(
        "--keep-ions",
        action="store_true",
        help="Keep metal ions"
    )
    
    parser.add_argument(
        "--keep-solvents", 
        action="store_true",
        help="Keep solvents and buffer molecules"
    )
    
    parser.add_argument(
        "--keep-cofactors",
        action="store_true",
        help="Keep cofactors (HEM, FAD, NAD, etc.)"
    )
    
    args = parser.parse_args()
    
    # If no PDB file provided, help user find one
    if not args.pdb_file:
        print("🔍 No PDB file specified. Looking for PDB files in current directory...")
        pdb_files = [f for f in os.listdir('.') if f.endswith(('.pdb', '.PDB'))]
        
        if not pdb_files:
            print("❌ No PDB files found in current directory.")
            print("💡 Please provide a PDB file:")
            print("   python protein_preprocessor.py your_file.pdb")
            sys.exit(1)
        
        print(f"📁 Found {len(pdb_files)} PDB file(s):")
        for i, pdb_file in enumerate(pdb_files, 1):
            print(f"   {i}. {pdb_file}")
        
        try:
            choice = input(f"\n🎯 Select file (1-{len(pdb_files)}): ")
            choice_idx = int(choice) - 1
            if 0 <= choice_idx < len(pdb_files):
                args.pdb_file = pdb_files[choice_idx]
                print(f"✅ Selected: {args.pdb_file}")
            else:
                print("❌ Invalid choice.")
                sys.exit(1)
        except (ValueError, KeyboardInterrupt):
            print("\n❌ Cancelled.")
            sys.exit(1)
    
    # Check if selected file exists
    if not Path(args.pdb_file).exists():
        print(f"❌ ERROR: PDB file not found: {args.pdb_file}")
        sys.exit(1)
    
    # Build policy from command line flags
    policy = {
        "keep_waters": args.keep_waters,
        "keep_ions": args.keep_ions, 
        "keep_solvents": args.keep_solvents,
        "keep_cofactors": args.keep_cofactors
    }
    
    try:
        result = prepare_protein(
            args.pdb_file, 
            interactive=args.interactive,
            default_policy=policy
        )
        
        if result:
            print(f"\n🎉 SUCCESS! Prepared PDBQT: {result}")
            sys.exit(0)
        else:
            print("\n❌ FAILED: Protein preparation unsuccessful")
            sys.exit(1)
            
    except Exception as e:
        print(f"\n❌ UNEXPECTED ERROR: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
