# backend/docking_pipeline.py
# Complete molecular docking pipeline with workspace-only file creation and robust pocket selection

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
from modules.pocket_identifier import identify_binding_site, get_pocket_analysis_config
from modules.docking_engine import DockingEngine
from modules.result_parser import parse_vina_results
from modules.protein_preprocessor import ProteinPreprocessor

class CompleteFileTracker:
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
        if file_path.exists() and not self.is_in_workspace(file_path):
            self.external_files.append(file_path)
    def is_in_workspace(self, file_path: Path) -> bool:
        try:
            file_path.relative_to(self.workspace)
            return True
        except ValueError:
            return False
    def find_external_files_by_timestamp(self, timestamp: str):
        import re
        pattern = re.compile(rf".*{timestamp}.*\.(pdb|pdbqt|json|sdf|mol2)$", re.IGNORECASE)
        for data_dir in self.data_dirs_to_check:
            if data_dir.exists():
                for fp in data_dir.rglob("*"):
                    if fp.is_file() and pattern.match(fp.name):
                        self.track_external_file(fp)
    def cleanup_all_tracked_files(self):
        cleaned = 0
        for file_path in self.external_files:
            try:
                if file_path.exists():
                    file_path.unlink()
                    cleaned += 1
                    print(f"🗑️  Removed: {file_path}")
            except Exception as e:
                print(f"Warning: Could not remove {file_path}: {e}")
        if self.workspace.exists():
            shutil.rmtree(self.workspace)
            print(f"🗑️  Removed workspace: {self.workspace}")
        if cleaned > 0:
            print(f"🧹 Total cleaned: {cleaned} external files + workspace")
    def get_total_size(self) -> float:
        total_size = 0
        for file_path in self.external_files:
            try:
                if file_path.exists():
                    total_size += file_path.stat().st_size
            except:
                pass
        if self.workspace.exists():
            total_size += sum(f.stat().st_size for f in self.workspace.rglob('*') if f.is_file())
        return total_size / (1024 * 1024)

def prepare_protein_in_workspace(protein_pdb_path: str, workspace_prepared_dir: Path,
                                 run_timestamp: str, cleaning_policy: Dict = None) -> str:
    workspace_prepared_dir.mkdir(parents=True, exist_ok=True)
    processor = ProteinPreprocessor(Path(protein_pdb_path))
    processor.output_dir = workspace_prepared_dir
    processor.base_stem = f"{Path(protein_pdb_path).stem}_{run_timestamp}"
    processor.tmp_cleaned_pdb = workspace_prepared_dir / f"{processor.base_stem}_clean.pdb"
    processor.clean_report = workspace_prepared_dir / f"{processor.base_stem}_clean_report.json"
    processor.output_pdbqt = workspace_prepared_dir / f"{processor.base_stem}_prepared.pdbqt"
    if not cleaning_policy:
        cleaning_policy = {"keep_waters": False,"keep_ions": True,"keep_solvents": False,"keep_cofactors": True}
    result = processor.process(interactive=False, default_policy=cleaning_policy)
    return result

class CleanDockingPipeline:
    def __init__(self, pipeline_id: str = None):
        self.pipeline_id = pipeline_id or str(uuid.uuid4())[:8]
        self.run_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.workspace = Path(f"temp_runs/{self.run_timestamp}_{self.pipeline_id}")
        self.raw_dir = self.workspace / "raw"
        self.prepared_dir = self.workspace / "prepared"
        self.docking_dir = self.workspace / "docking"
        for directory in [self.workspace, self.raw_dir, self.prepared_dir, self.docking_dir]:
            directory.mkdir(parents=True, exist_ok=True)
        self.file_tracker = CompleteFileTracker(self.workspace)
        print(f"🚀 Docking Pipeline [{self.pipeline_id}]")
        print(f"   Workspace: {self.workspace}")
        self._cleanup_old_temp_runs()

    def _cleanup_old_temp_runs(self, days_old: int = 7):
        temp_runs_dir = Path("temp_runs")
        if not temp_runs_dir.exists(): return
        cutoff = time.time() - (days_old * 24 * 3600)
        cleaned = 0
        for item in temp_runs_dir.iterdir():
            if item.is_dir() and item.stat().st_mtime < cutoff:
                shutil.rmtree(item, ignore_errors=True)
                cleaned += 1
        if cleaned: print(f"🧹 Cleaned {cleaned} old temp runs")

    def update_progress(self, percent: int, message: str):
        bar_length = 30
        filled_length = int(bar_length * percent // 100)
        bar = '█' * filled_length + '░' * (bar_length - filled_length)
        print(f"\rProgress: |{bar}| {percent}% - {message}", end='', flush=True)
        progress_data = {"pipeline_id": self.pipeline_id,"percent": percent,"message": message,"timestamp": datetime.now().isoformat()}
        (self.workspace / "progress.json").write_text(json.dumps(progress_data))

    def get_inputs(self, protein_input: str = None, ligand_input: str = None):
        if not protein_input:
            pdbs = [f for f in os.listdir('.') if f.lower().endswith('.pdb')]
            if pdbs: print(f"\n📁 Available: {', '.join(pdbs[:3])}{'...' if len(pdbs)>3 else ''}")
            protein_input = input("🧬 Enter protein: ").strip()
        if not ligand_input:
            ligs = [f for f in os.listdir('.') if f.lower().endswith(('.sdf','.mol','.mol2','.pdbqt'))]
            if ligs: print(f"📁 Available: {', '.join(ligs[:3])}{'...' if len(ligs)>3 else ''}")
            ligand_input = input("💊 Enter ligand: ").strip()
        return protein_input, ligand_input

    def resolve_protein_input(self, protein_input: str) -> str:
        p = Path(protein_input)
        if p.exists() and p.is_file():
            dst = self.raw_dir / p.name
            dst.write_bytes(p.read_bytes())
            return str(dst)
        fetcher = ProteinFetcher(download_dir=str(self.raw_dir))
        return fetcher.fetch_from_pdb(protein_input) if (len(protein_input)==4 and protein_input.isalnum()) else fetcher.fetch_from_alphafold(protein_input)

    def resolve_ligand_input(self, ligand_input: str) -> str:
        p = Path(ligand_input)
        if p.exists() and p.is_file():
            dst = self.raw_dir / p.name
            dst.write_bytes(p.read_bytes())
            if p.suffix.lower()==".pdbqt": return str(dst)
            if p.suffix.lower() in {".sdf",".mol2",".pdb",".mol"}:
                return _convert_to_pdbqt(str(dst), str(self.prepared_dir))
        return process_ligand(ligand_input, ligand_dir=str(self.raw_dir), prepared_dir=str(self.prepared_dir))

    def display_binding_sites(self, pocket_info: Dict):
        print(f"\n\n🎯 TOP BINDING SITES")
        print("="*60)
        all_sites=[]
        if "primary" in pocket_info: all_sites.append(pocket_info["primary"])
        if "modes" in pocket_info: all_sites.extend(pocket_info["modes"])
        for i, s in enumerate(all_sites[:5], 1):
            print(f"Site {i}: ({s.get('center_x',0):.1f}, {s.get('center_y',0):.1f}, {s.get('center_z',0):.1f}) "
                  f"Box: {s.get('size_x',20):.0f}×{s.get('size_y',20):.0f}×{s.get('size_z',20):.0f} Å  "
                  f"via {s.get('method','?')}  final={s.get('final_score','-')}")
        print("="*60)

    def _protein_bounds_ca(self, pdb_path: str):
        xs=[]; ys=[]; zs=[]
        with open(pdb_path) as f:
            for line in f:
                if line.startswith("ATOM") and line[12:16].strip()=="CA":
                    parts = line[30:54].split()
                    if len(parts)>=3:
                        x,y,z = map(float, parts[:3])
                        xs.append(x); ys.append(y); zs.append(z)
        if not xs: return None
        return (min(xs),max(xs)),(min(ys),max(ys)),(min(zs),max(zs))

    def _clamp_box_to_protein(self, raw_protein_path: str, center: Dict) -> Dict:
        bounds = self._protein_bounds_ca(raw_protein_path)
        if not bounds: return center
        (xmin,xmax),(ymin,ymax),(zmin,zmax) = bounds
        cx,cy,cz = center["center_x"], center["center_y"], center["center_z"]
        sx,sy,sz = center["size_x"], center["size_y"], center["size_z"]
        def clamp_center(c, s, lo, hi):
            half = s/2.0
            if c+half < lo: c = lo + half + 1.5
            if c-half > hi: c = hi - half - 1.5
            return c
        cx = clamp_center(cx, sx, xmin, xmax)
        cy = clamp_center(cy, sy, ymin, ymax)
        cz = clamp_center(cz, sz, zmin, zmax)
        sx = max(14.0, min(26.0, sx)); sy = max(14.0, min(26.0, sy)); sz = max(14.0, min(26.0, sz))
        center.update({"center_x":cx,"center_y":cy,"center_z":cz,"size_x":sx,"size_y":sy,"size_z":sz})
        return center

    def run_docking(self, protein_input: str = None, ligand_input: str = None,
                    cleaning_policy: Optional[Dict[str,bool]] = None, detailed: bool = True):
        start_time = time.time()
        try:
            # Inputs
            self.update_progress(5, "Collecting inputs")
            protein_input, ligand_input = self.get_inputs(protein_input, ligand_input)

            # Resolve inputs
            self.update_progress(15, "Resolving protein")
            raw_protein_path = self.resolve_protein_input(protein_input)

            self.update_progress(25, "Resolving ligand")
            prepared_ligand_pdbqt = self.resolve_ligand_input(ligand_input)

            # Prepare protein in workspace
            self.update_progress(45, "Preparing protein (workspace-only)")
            prepared_protein_pdbqt = prepare_protein_in_workspace(raw_protein_path, self.prepared_dir, self.run_timestamp, cleaning_policy)
            if not prepared_protein_pdbqt:
                raise Exception("Protein preparation failed")

            # Pocket identification with consensus/template/geometric
            self.update_progress(55, "Finding binding sites")
            pocket_info = identify_binding_site(
                raw_protein_path, prepared_protein_pdbqt, prepared_ligand_pdbqt,
                **get_pocket_analysis_config("production" if detailed else "fast")
            )
            if not pocket_info or "primary" not in pocket_info:
                raise Exception("No binding sites found")

            # Show and clamp box to protein bounds
            self.display_binding_sites(pocket_info)
            center = pocket_info["primary"]
            center = self._clamp_box_to_protein(raw_protein_path, center)

            # Optional: probe alternate
            modes = pocket_info.get("modes", [])
            if len(modes) >= 2:
                self.update_progress(63, "Probing alternate pocket")
                def quick_probe(c):
                    try:
                        eng = DockingEngine(prepared_protein_pdbqt, prepared_ligand_pdbqt)
                        return eng.run_docking(center_x=c["center_x"], center_y=c["center_y"], center_z=c["center_z"],
                                               box_size_x=c["size_x"], box_size_y=c["size_y"], box_size_z=c["size_z"],
                                               output_dir=str(self.docking_dir), quick=True)
                    except Exception:
                        return None
                q1 = quick_probe(center)
                q2 = quick_probe(modes[1])
                def best_aff(path):
                    if not path or not Path(path).exists(): return None
                    vals=[]
                    try:
                        with open(path) as f: 
                            for line in f:
                                if "REMARK VINA RESULT" in line:
                                    parts = line.split()
                                    vals.append(float(parts[3]))
                    except Exception:
                        return None
                    return min(vals) if vals else None
                a1 = best_aff(q1); a2 = best_aff(q2)
                if a2 is not None and (a1 is None or a2 < a1):
                    center = modes[1]
                    center = self._clamp_box_to_protein(raw_protein_path, center)
                    print(f"\n🔁 Switched to alternate pocket (quick better: {a2:.2f} < {a1 if a1 is not None else float('inf')})")

            # Docking
            self.update_progress(65, "Running molecular docking")
            engine = DockingEngine(prepared_protein_pdbqt, prepared_ligand_pdbqt)
            docked_file = engine.run_docking(
                center_x=center["center_x"], center_y=center["center_y"], center_z=center["center_z"],
                box_size_x=center["size_x"], box_size_y=center["size_y"], box_size_z=center["size_z"],
                output_dir=str(self.docking_dir)
            )
            if not docked_file:
                raise Exception("Docking failed")

            # Results
            self.update_progress(85, "Analyzing results")
            results_csv = self.docking_dir / "binding_scores.csv"
            df = parse_vina_results(docked_file, output_csv_path=results_csv, pocket_info=pocket_info)
            if df is None or df.empty:
                raise Exception("Results parsing failed")

            best_affinity = float(df["Binding Affinity (kcal/mol)"].min())
            total_poses = len(df)
            total_time = time.time() - start_time

            # Summary
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
            (self.workspace / "summary.json").write_text(json.dumps(summary_data, indent=2, default=str))

            self.update_progress(100, "Pipeline completed")
            print(f"\n\n🎉 DOCKING COMPLETED!")
            print(f"   Best Affinity: {best_affinity:.2f} kcal/mol")
            print(f"   Total Poses: {total_poses}")
            print(f"   Duration: {total_time:.1f}s")

            self.handle_data_retention(summary_data)
            return summary_data

        except Exception as e:
            total_time = time.time() - start_time
            print(f"\n❌ Pipeline failed after {total_time:.1f}s: {str(e)}")
            choice = input("\n🗑️  Delete all files? [Y/n]: ").strip().lower()
            if choice != 'n':
                self.file_tracker.cleanup_all_tracked_files()
                print("🧹 All files cleaned up")
            return None

    # Data retention methods (unchanged)
    def handle_data_retention(self, summary_data: Dict) -> bool:
        print(f"\n\n" + "="*80)
        print("💾 DATA MANAGEMENT - ALL FILES TRACKED")
        print("="*80)
        self.file_tracker.find_external_files_by_timestamp(self.run_timestamp)
        results = summary_data.get('results', {})
        print(f"   • Best Affinity: {results.get('best_affinity', 'N/A')} kcal/mol")
        print(f"   • Total Poses: {results.get('total_poses', 'N/A')}")
        print(f"   • Duration: {summary_data['pipeline_info']['duration']}")
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
            if choice == '1': return self._save_permanently(summary_data)
            elif choice == '2': return self._delete_everything()
            elif choice == '3': return self._keep_temporarily()
            else: print("❌ Invalid choice. Please enter 1, 2, or 3.")

    def _save_permanently(self, summary_data: Dict) -> bool:
        try:
            default_name = f"docking_{summary_data['inputs']['protein']}_{summary_data['inputs']['ligand']}_{datetime.now().strftime('%Y%m%d')}"
            project_name = input(f"\n📝 Project name [{default_name}]: ").strip() or default_name
            permanent_path = Path(f"data/results/{project_name}")
            permanent_path.mkdir(parents=True, exist_ok=True)
            shutil.copytree(self.workspace, permanent_path / "workspace")
            external_dest = permanent_path / "external_files"
            external_dest.mkdir(exist_ok=True)
            for ext_file in self.file_tracker.external_files:
                if ext_file.exists():
                    shutil.copy2(ext_file, external_dest / ext_file.name)
            summary_data['pipeline_info']['permanent_location'] = str(permanent_path)
            (permanent_path / "summary.json").write_text(json.dumps(summary_data, indent=2, default=str))
            self.file_tracker.cleanup_all_tracked_files()
            print(f"✅ Everything saved to: {permanent_path}")
            return True
        except Exception as e:
            print(f"❌ Error saving: {e}")
            return False

    def _delete_everything(self) -> bool:
        try:
            total_size = self.file_tracker.get_total_size()
            print(f"\n⚠️  This will permanently delete {total_size:.2f} MB")
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
        print(f"⏳ All files kept temporarily (auto-delete policy applies)")
        return True

def main():
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
