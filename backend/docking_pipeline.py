# backend/docking_pipeline.py
# Complete molecular docking pipeline with workspace-only file creation

import os
import sys
import argparse
import json
import time
import uuid
import shutil
from pathlib import Path
from datetime import datetime
from typing import Dict, Optional, List

from modules.protein_fetcher import ProteinFetcher
from modules.ligand_preparer import process_ligand, _convert_to_pdbqt
from modules.pocket_identifier import identify_binding_site
from modules.docking_engine import DockingEngine
from modules.result_parser import parse_vina_results

# Import the modified protein preprocessor
from modules.protein_preprocessor import ProteinPreprocessor


class CompleteFileTracker:
    """Tracks ALL files created during pipeline execution"""
    
    def __init__(self, workspace: Path):
        self.workspace = workspace
        self.external_files: List[Path] = []
        self.data_dirs_to_check = [
            Path("data/prepared_proteins"),
            Path("data/prepared_ligands"), 
            Path("data/proteins"),
            Path("data/ligands")
        ]
    
    def track_external_file(self, file_path: Path):
        """Track files created outside workspace"""
        if file_path.exists() and not self.is_in_workspace(file_path):
            self.external_files.append(file_path)
    
    def is_in_workspace(self, file_path: Path) -> bool:
        """Check if file is within pipeline workspace"""
        try:
            file_path.relative_to(self.workspace)
            return True
        except ValueError:
            return False
    
    def find_external_files_by_timestamp(self, timestamp: str):
        """Find files created during this run outside workspace"""
        import re
        pattern = re.compile(rf".*{timestamp}.*\.(pdb|pdbqt|json|sdf|mol2)$", re.IGNORECASE)
        
        for data_dir in self.data_dirs_to_check:
            if data_dir.exists():
                for file_path in data_dir.rglob("*"):
                    if file_path.is_file() and pattern.match(file_path.name):
                        self.track_external_file(file_path)
    
    def cleanup_all_tracked_files(self):
        """Remove all external files created during this run"""
        cleaned = 0
        for file_path in self.external_files:
            try:
                if file_path.exists():
                    file_path.unlink()
                    cleaned += 1
                    print(f"🗑️  Removed: {file_path}")
            except Exception as e:
                print(f"Warning: Could not remove {file_path}: {e}")
        
        # Remove workspace
        if self.workspace.exists():
            shutil.rmtree(self.workspace)
            print(f"🗑️  Removed workspace: {self.workspace}")
        
        if cleaned > 0:
            print(f"🧹 Total cleaned: {cleaned} external files + workspace")
    
    def get_total_size(self) -> float:
        """Get total size of all tracked files in MB"""
        total_size = 0
        
        # External files
        for file_path in self.external_files:
            try:
                if file_path.exists():
                    total_size += file_path.stat().st_size
            except:
                pass
        
        # Workspace
        if self.workspace.exists():
            total_size += sum(f.stat().st_size for f in self.workspace.rglob('*') if f.is_file())
        
        return total_size / (1024 * 1024)


def prepare_protein_in_workspace(protein_pdb_path: str, workspace_prepared_dir: Path, 
                                run_timestamp: str, cleaning_policy: Dict = None) -> str:
    """
    Prepare protein with output directed to pipeline workspace
    """
    workspace_prepared_dir.mkdir(parents=True, exist_ok=True)
    
    # Create a custom protein preprocessor that outputs to workspace
    processor = ProteinPreprocessor(Path(protein_pdb_path))
    
    # Override the output directory to use workspace
    original_output_dir = processor.output_dir
    processor.output_dir = workspace_prepared_dir
    
    # Update file paths to use workspace
    processor.base_stem = f"{Path(protein_pdb_path).stem}_{run_timestamp}"
    processor.tmp_cleaned_pdb = workspace_prepared_dir / f"{processor.base_stem}_clean.pdb"
    processor.clean_report = workspace_prepared_dir / f"{processor.base_stem}_clean_report.json"
    processor.output_pdbqt = workspace_prepared_dir / f"{processor.base_stem}_prepared.pdbqt"
    
    if not cleaning_policy:
        cleaning_policy = {
            "keep_waters": False,
            "keep_ions": True,
            "keep_solvents": False,
            "keep_cofactors": True
        }
    
    # Process with workspace output
    result = processor.process(interactive=False, default_policy=cleaning_policy)
    
    return result


class CleanDockingPipeline:
    """Production-ready docking pipeline with complete file control"""
    
    def __init__(self, pipeline_id: str = None):
        self.pipeline_id = pipeline_id or str(uuid.uuid4())[:8]
        self.run_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Create workspace
        self.workspace = Path(f"temp_runs/{self.run_timestamp}_{self.pipeline_id}")
        self.raw_dir = self.workspace / "raw"
        self.prepared_dir = self.workspace / "prepared"
        self.docking_dir = self.workspace / "docking"
        
        # Create all directories
        for directory in [self.workspace, self.raw_dir, self.prepared_dir, self.docking_dir]:
            directory.mkdir(parents=True, exist_ok=True)
        
        # Initialize file tracker
        self.file_tracker = CompleteFileTracker(self.workspace)
        
        print(f"🚀 Docking Pipeline [{self.pipeline_id}]")
        print(f"   Workspace: {self.workspace}")
        
        # Cleanup old temp files
        self._cleanup_old_temp_runs()
    
    def _cleanup_old_temp_runs(self, days_old: int = 7):
        """Cleanup old temporary runs"""
        temp_runs_dir = Path("temp_runs")
        if not temp_runs_dir.exists():
            return
        
        cutoff_time = time.time() - (days_old * 24 * 3600)
        cleaned = 0
        
        for item in temp_runs_dir.iterdir():
            if item.is_dir() and item.stat().st_mtime < cutoff_time:
                shutil.rmtree(item, ignore_errors=True)
                cleaned += 1
        
        if cleaned > 0:
            print(f"🧹 Cleaned {cleaned} old temp runs")
    
    def update_progress(self, percent: int, message: str):
        """Update progress with visual bar"""
        bar_length = 30
        filled_length = int(bar_length * percent // 100)
        bar = '█' * filled_length + '░' * (bar_length - filled_length)
        
        print(f"\rProgress: |{bar}| {percent}% - {message}", end='', flush=True)
        
        # Emit for frontend
        progress_data = {
            "pipeline_id": self.pipeline_id,
            "percent": percent,
            "message": message,
            "timestamp": datetime.now().isoformat()
        }
        
        progress_file = self.workspace / "progress.json"
        with open(progress_file, 'w') as f:
            json.dump(progress_data, f)
    
    def get_inputs(self, protein_input: str = None, ligand_input: str = None):
        """Get inputs with prompts"""
        if not protein_input:
            pdb_files = [f for f in os.listdir('.') if f.endswith(('.pdb', '.PDB'))]
            if pdb_files:
                print(f"\n📁 Available: {', '.join(pdb_files[:3])}{'...' if len(pdb_files) > 3 else ''}")
            protein_input = input("🧬 Enter protein: ").strip()
        
        if not ligand_input:
            ligand_files = [f for f in os.listdir('.') if f.endswith(('.sdf', '.mol', '.mol2', '.pdbqt'))]
            if ligand_files:
                print(f"📁 Available: {', '.join(ligand_files[:3])}{'...' if len(ligand_files) > 3 else ''}")
            ligand_input = input("💊 Enter ligand: ").strip()
        
        return protein_input, ligand_input
    
    def resolve_protein_input(self, protein_input: str) -> str:
        """Resolve protein and store in workspace"""
        p = Path(protein_input)
        if p.exists() and p.is_file():
            dst = self.raw_dir / p.name
            dst.write_bytes(p.read_bytes())
            return str(dst)
        
        # Fetch to workspace
        fetcher = ProteinFetcher(download_dir=str(self.raw_dir))
        
        if len(protein_input) == 4 and protein_input.isalnum():
            result = fetcher.fetch_from_pdb(protein_input)
        else:
            result = fetcher.fetch_from_alphafold(protein_input)
        
        return result
    
    def resolve_ligand_input(self, ligand_input: str) -> str:
        """Resolve ligand and store in workspace"""
        p = Path(ligand_input)
        if p.exists() and p.is_file():
            dst = self.raw_dir / p.name
            dst.write_bytes(p.read_bytes())
            
            if p.suffix.lower() == ".pdbqt":
                return str(dst)
            elif p.suffix.lower() in {".sdf", ".mol2", ".pdb", ".mol"}:
                result = _convert_to_pdbqt(str(dst), str(self.prepared_dir))
                return result
        
        # Fetch to workspace
        result = process_ligand(ligand_input, ligand_dir=str(self.raw_dir), prepared_dir=str(self.prepared_dir))
        return result
    
    def display_binding_sites(self, pocket_info: Dict):
        """Display top 5 binding sites"""
        print(f"\n\n🎯 TOP 5 BINDING SITES:")
        print("="*60)
        
        all_sites = []
        if "primary" in pocket_info:
            all_sites.append(pocket_info["primary"])
        if "modes" in pocket_info:
            all_sites.extend(pocket_info["modes"])
        
        for i, site in enumerate(all_sites[:5], 1):
            print(f"Site {i}: ({site.get('center_x', 0):.1f}, {site.get('center_y', 0):.1f}, {site.get('center_z', 0):.1f}) "
                  f"Box: {site.get('size_x', 20):.0f}×{site.get('size_y', 20):.0f}×{site.get('size_z', 20):.0f}Å")
        
        if len(all_sites) > 5:
            print(f"... and {len(all_sites) - 5} more sites")
        print("="*60)
    
    def handle_data_retention(self, summary_data: Dict) -> bool:
        """Complete data retention with external file tracking"""
        print(f"\n\n" + "="*80)
        print("💾 DATA MANAGEMENT - ALL FILES TRACKED")
        print("="*80)
        
        # Find any external files created during this run
        self.file_tracker.find_external_files_by_timestamp(self.run_timestamp)
        
        print(f"📊 Results Summary:")
        results = summary_data.get('results', {})
        print(f"   • Best Affinity: {results.get('best_affinity', 'N/A')} kcal/mol")
        print(f"   • Total Poses: {results.get('total_poses', 'N/A')}")
        print(f"   • Duration: {summary_data['pipeline_info']['duration']}")
        
        # Show what will be cleaned
        total_size = self.file_tracker.get_total_size()
        print(f"\n📁 ALL FILES TO BE MANAGED:")
        print(f"   • Workspace: {self.workspace}")
        print(f"   • External files found: {len(self.file_tracker.external_files)}")
        for ext_file in self.file_tracker.external_files:
            print(f"     - {ext_file}")
        print(f"   • Total size: {total_size:.2f} MB")
        
        while True:
            print(f"\nWhat would you like to do?")
            print(f"1. 💾 Save Permanently (move everything to results)")
            print(f"2. 🗑️  Delete Everything (remove ALL tracked files)")
            print(f"3. ⏳ Keep Temporarily (auto-delete in 7 days)")
            
            choice = input("\nChoose (1-3): ").strip()
            
            if choice == '1':
                return self._save_permanently(summary_data)
            elif choice == '2':
                return self._delete_everything()
            elif choice == '3':
                return self._keep_temporarily()
            else:
                print("❌ Invalid choice. Please enter 1, 2, or 3.")
    
    def _save_permanently(self, summary_data: Dict) -> bool:
        """Save everything permanently"""
        try:
            default_name = f"docking_{summary_data['inputs']['protein']}_{summary_data['inputs']['ligand']}_{datetime.now().strftime('%Y%m%d')}"
            project_name = input(f"\n📝 Project name [{default_name}]: ").strip()
            if not project_name:
                project_name = default_name
            
            # Create permanent location
            permanent_path = Path(f"data/results/{project_name}")
            permanent_path.mkdir(parents=True, exist_ok=True)
            
            # Copy workspace
            workspace_dest = permanent_path / "workspace"
            shutil.copytree(self.workspace, workspace_dest)
            
            # Copy external files
            external_dest = permanent_path / "external_files"
            external_dest.mkdir(exist_ok=True)
            for ext_file in self.file_tracker.external_files:
                if ext_file.exists():
                    dest_file = external_dest / ext_file.name
                    shutil.copy2(ext_file, dest_file)
            
            # Save summary
            summary_data['pipeline_info']['permanent_location'] = str(permanent_path)
            summary_file = permanent_path / "summary.json"
            with open(summary_file, 'w') as f:
                json.dump(summary_data, f, indent=2, default=str)
            
            print(f"✅ Everything saved to: {permanent_path}")
            
            # Now cleanup original files
            self.file_tracker.cleanup_all_tracked_files()
            
            return True
            
        except Exception as e:
            print(f"❌ Error saving: {e}")
            return False
    
    def _delete_everything(self) -> bool:
        """Delete everything"""
        try:
            total_size = self.file_tracker.get_total_size()
            print(f"\n⚠️  This will permanently delete {total_size:.2f} MB:")
            print(f"   • Workspace: {self.workspace}")
            print(f"   • External files: {len(self.file_tracker.external_files)} files")
            
            confirm = input(f"\nType 'DELETE' to confirm: ").strip()
            if confirm == 'DELETE':
                self.file_tracker.cleanup_all_tracked_files()
                print(f"🗑️  Everything deleted ({total_size:.2f} MB freed)")
                return True
            else:
                print("❌ Deletion cancelled")
                return self._keep_temporarily()
        except Exception as e:
            print(f"❌ Error deleting: {e}")
            return False
    
    def _keep_temporarily(self) -> bool:
        """Keep temporarily"""
        print(f"⏳ All files kept temporarily")
        print(f"🗑️  Will be auto-deleted in 7 days")
        return True
    
    def run_docking(self, protein_input: str = None, ligand_input: str = None, 
                   cleaning_policy: Optional[Dict[str,bool]] = None, detailed: bool = True):
        """Run complete docking with workspace-only files"""
        start_time = time.time()
        
        try:
            # Step 1: Inputs
            self.update_progress(5, "Collecting inputs")
            protein_input, ligand_input = self.get_inputs(protein_input, ligand_input)
            
            # Step 2: Resolve inputs
            self.update_progress(15, "Resolving protein")
            raw_protein_path = self.resolve_protein_input(protein_input)
            
            self.update_progress(25, "Resolving ligand")
            prepared_ligand_pdbqt = self.resolve_ligand_input(ligand_input)
            
            # Step 3: Prepare protein IN WORKSPACE ONLY
            self.update_progress(45, "Preparing protein (workspace-only)")
            prepared_protein_pdbqt = prepare_protein_in_workspace(
                raw_protein_path, 
                self.prepared_dir, 
                self.run_timestamp, 
                cleaning_policy
            )
            
            if not prepared_protein_pdbqt:
                raise Exception("Protein preparation failed")
            
            # Step 4: Find binding sites
            self.update_progress(55, "Finding binding sites")
            pocket_info = identify_binding_site(
                raw_protein_path, prepared_protein_pdbqt, prepared_ligand_pdbqt,
                use_validated=True, return_n=5, detailed=detailed
            )
            
            if not pocket_info or "primary" not in pocket_info:
                raise Exception("No binding sites found")
            
            self.display_binding_sites(pocket_info)
            center = pocket_info["primary"]
            
            # Step 5: Run docking
            self.update_progress(65, "Running molecular docking")
            engine = DockingEngine(prepared_protein_pdbqt, prepared_ligand_pdbqt)
            docked_file = engine.run_docking(
                center_x=center["center_x"], center_y=center["center_y"], center_z=center["center_z"],
                box_size_x=center["size_x"], box_size_y=center["size_y"], box_size_z=center["size_z"],
                output_dir=str(self.docking_dir)
            )
            
            if not docked_file:
                raise Exception("Docking failed")
            
            # Step 6: Parse results
            self.update_progress(85, "Analyzing results")
            results_csv = self.docking_dir / "binding_scores.csv"
            df = parse_vina_results(docked_file, output_csv_path=results_csv, pocket_info=pocket_info)
            
            if df is None or df.empty:
                raise Exception("Results parsing failed")
            
            best_affinity = float(df["Binding Affinity (kcal/mol)"].min())
            total_poses = len(df)
            total_time = time.time() - start_time
            
            # Step 7: Summary
            self.update_progress(95, "Generating summary")
            summary_data = {
                "pipeline_info": {
                    "id": self.pipeline_id,
                    "run_timestamp": self.run_timestamp,
                    "workspace": str(self.workspace),
                    "duration": f"{total_time:.1f}s",
                    "timestamp": datetime.now().isoformat()
                },
                "inputs": {"protein": protein_input, "ligand": ligand_input},
                "binding_sites": pocket_info,
                "results": {
                    "best_affinity": best_affinity,
                    "total_poses": total_poses,
                    "poses_file": docked_file,
                    "results_csv": str(results_csv)
                },
                "files": {
                    "raw_protein": raw_protein_path,
                    "prepared_protein": prepared_protein_pdbqt,
                    "prepared_ligand": prepared_ligand_pdbqt,
                    "docked_poses": docked_file,
                    "results_csv": str(results_csv)
                }
            }
            
            # Save summary
            temp_summary_file = self.workspace / "summary.json"
            with open(temp_summary_file, 'w') as f:
                json.dump(summary_data, f, indent=2, default=str)
            
            # Completion
            self.update_progress(100, "Pipeline completed")
            print(f"\n\n🎉 DOCKING COMPLETED!")
            print(f"   Best Affinity: {best_affinity:.2f} kcal/mol")
            print(f"   Total Poses: {total_poses}")
            print(f"   Duration: {total_time:.1f}s")
            
            # Handle data retention with complete tracking
            self.handle_data_retention(summary_data)
            
            return summary_data
            
        except Exception as e:
            total_time = time.time() - start_time
            print(f"\n❌ Pipeline failed after {total_time:.1f}s: {str(e)}")
            
            # Cleanup on failure
            choice = input("\n🗑️  Delete all files? [Y/n]: ").strip().lower()
            if choice != 'n':
                self.file_tracker.cleanup_all_tracked_files()
                print("🧹 All files cleaned up")
            
            return None


def main():
    """CLI entry point"""
    parser = argparse.ArgumentParser(description="Clean Molecular Docking Pipeline")
    parser.add_argument("--protein", "-p", help="Protein input")
    parser.add_argument("--ligand", "-l", help="Ligand input")
    parser.add_argument("--keep-waters", action="store_true", help="Keep waters")
    parser.add_argument("--keep-ions", action="store_true", help="Keep ions")
    parser.add_argument("--keep-solvents", action="store_true", help="Keep solvents")
    parser.add_argument("--keep-cofactors", action="store_true", help="Keep cofactors")
    parser.add_argument("--detailed", action="store_true", help="Detailed pocket analysis")
    
    args = parser.parse_args()
    
    cleaning_policy = None
    if any([args.keep_waters, args.keep_ions, args.keep_solvents, args.keep_cofactors]):
        cleaning_policy = {
            "keep_waters": args.keep_waters,
            "keep_ions": args.keep_ions,
            "keep_solvents": args.keep_solvents,
            "keep_cofactors": args.keep_cofactors
        }
    
    try:
        pipeline = CleanDockingPipeline()
        result = pipeline.run_docking(
            protein_input=args.protein,
            ligand_input=args.ligand,
            cleaning_policy=cleaning_policy,
            detailed=args.detailed
        )
        
        sys.exit(0 if result else 1)
        
    except KeyboardInterrupt:
        print("\n❌ Cancelled by user")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Fatal error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
