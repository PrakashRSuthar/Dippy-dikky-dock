# backend/docking_pipeline.py

import sys
import argparse
from api.docking_api import DockingPipeline, update_job_status
import uuid
from datetime import datetime

def console_status_callback(job_id: str, status: str, progress: int, message: str):
    """Print status updates to console"""
    bar_length = 30
    filled_length = int(bar_length * progress // 100)
    bar = '█' * filled_length + '-' * (bar_length - filled_length)
    print(f'\r[{bar}] {progress}% - {message}', end='', flush=True)
    if status in ["completed", "failed"]:
        print()  # New line when done

def main():
    parser = argparse.ArgumentParser(description="Dippy Docking Pipeline")
    parser.add_argument("--protein", required=True, help="Protein ID or file path")
    parser.add_argument("--ligand", required=True, help="Ligand name or file path")
    parser.add_argument("--api", action="store_true", help="Start API server instead")
    
    args = parser.parse_args()
    
    if args.api:
        print("🚀 Starting Dippy Docking API server...")
        import uvicorn
        from api.docking_api import app
        uvicorn.run(app, host="0.0.0.0", port=8000)
        return
    
    print("🦾 === ULTIMATE DOCKING BEAST (CLI Mode) ===")
    print(f"Protein: {args.protein}")
    print(f"Ligand: {args.ligand}")
    
    # Run pipeline
    job_id = str(uuid.uuid4())
    pipeline = DockingPipeline()
    pipeline.set_status_callback(console_status_callback)
    
    result = pipeline.run_pipeline(job_id, args.protein, args.ligand)
    
    if result.get("status") == "completed":
        print(f"\n🏆 SUCCESS!")
        print(f"🥇 Best Affinity: {result['best_affinity']:.2f} kcal/mol")
        print(f"📁 Results: {result['run_dir']}")
    else:
        print(f"\n💥 FAILED: {result.get('error_message', 'Unknown error')}")

if __name__ == "__main__":
    if len(sys.argv) == 1:
        # Interactive mode
        print("🦾 === ULTIMATE DOCKING BEAST ===")
        protein = input("Protein (ID or file path): ").strip()
        ligand = input("Ligand (name or file path): ").strip()
        
        if not protein or not ligand:
            print("❌ Both protein and ligand required!")
            sys.exit(1)
        
        sys.argv.extend(["--protein", protein, "--ligand", ligand])
    
    main()
