# backend/modules/result_parser.py

import pandas as pd
from pathlib import Path
from typing import Dict, Optional


def parse_vina_results(docked_pdbqt_path: str, output_csv_path: str = None, pocket_info: Dict = None):
    """
    Parses a Vina output PDBQT file to extract binding affinity scores and coordinates for each pose.
    
    Args:
        docked_pdbqt_path (str): The path to the docked PDBQT file from Vina.
        output_csv_path (str, optional): Path to save the results as a CSV file.
        pocket_info (Dict, optional): Information about binding sites and coordinates.
    
    Returns:
        pandas.DataFrame: A DataFrame with pose information including coordinates and box sizes.
    """
    pdbqt_file = Path(docked_pdbqt_path)
    if not pdbqt_file.exists():
        print(f"❌ ERROR: Docked file not found at {pdbqt_file}")
        return pd.DataFrame()

    results = []
    pose_coordinates = []
    
    try:
        print(f"[INFO] Parsing results from {pdbqt_file.name}...")
        with open(pdbqt_file, 'r') as f:
            lines = f.readlines()
        
        current_pose = 0
        pose_atoms = []
        
        for line in lines:
            # Vina results are in lines like: "REMARK VINA RESULT:  -7.5      0.000    0.000"
            if line.startswith("REMARK VINA RESULT:"):
                parts = line.split()
                if len(parts) >= 4:
                    affinity = float(parts[3])
                    results.append(affinity)
                    current_pose += 1
                    
                    # Calculate center of mass for this pose if we have atoms
                    if pose_atoms:
                        center_x = sum(atom[0] for atom in pose_atoms) / len(pose_atoms)
                        center_y = sum(atom[1] for atom in pose_atoms) / len(pose_atoms)
                        center_z = sum(atom[2] for atom in pose_atoms) / len(pose_atoms)
                        pose_coordinates.append((center_x, center_y, center_z))
                        pose_atoms = []
            
            # Extract atomic coordinates for center of mass calculation
            elif line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip()) 
                    z = float(line[46:54].strip())
                    pose_atoms.append((x, y, z))
                except (ValueError, IndexError):
                    continue
            
            # Reset atoms list for new model
            elif line.startswith("MODEL"):
                pose_atoms = []

        # Handle last pose if it exists
        if pose_atoms and len(pose_coordinates) < len(results):
            center_x = sum(atom[0] for atom in pose_atoms) / len(pose_atoms)
            center_y = sum(atom[1] for atom in pose_atoms) / len(pose_atoms)
            center_z = sum(atom[2] for atom in pose_atoms) / len(pose_atoms)
            pose_coordinates.append((center_x, center_y, center_z))

        if not results:
            print(f"[WARN] No binding affinity scores found in {pdbqt_file.name}")
            return pd.DataFrame()

        # Create enhanced DataFrame
        data = {
            "Pose": range(1, len(results) + 1),
            "Binding Affinity (kcal/mol)": results
        }
        
        # Add coordinate information if available
        if pose_coordinates:
            if len(pose_coordinates) >= len(results):
                data["Pose Center X"] = [coord[0] for coord in pose_coordinates[:len(results)]]
                data["Pose Center Y"] = [coord[1] for coord in pose_coordinates[:len(results)]]
                data["Pose Center Z"] = [coord[2] for coord in pose_coordinates[:len(results)]]
        
        # Add binding site information if provided
        if pocket_info and "primary" in pocket_info:
            primary = pocket_info["primary"]
            data["Binding Site X"] = [primary.get("center_x", 0)] * len(results)
            data["Binding Site Y"] = [primary.get("center_y", 0)] * len(results)
            data["Binding Site Z"] = [primary.get("center_z", 0)] * len(results)
            data["Box Size X"] = [primary.get("size_x", 20)] * len(results)
            data["Box Size Y"] = [primary.get("size_y", 20)] * len(results)
            data["Box Size Z"] = [primary.get("size_z", 20)] * len(results)
        
        df = pd.DataFrame(data)
        
        # Display results
        print("\n" + "="*100)
        print("📊 DOCKING RESULTS SUMMARY")
        print("="*100)
        print(df.to_string(index=False, float_format='%.2f'))
        
        # Display top 5 binding sites if available
        if pocket_info:
            print(f"\n🎯 TOP 5 BINDING SITES DETAILS:")
            
            # Combine primary and modes
            all_sites = []
            if "primary" in pocket_info:
                all_sites.append(pocket_info["primary"])
            if "modes" in pocket_info:
                all_sites.extend(pocket_info["modes"])
            
            # Show top 5 sites
            for i, site in enumerate(all_sites[:5], 1):
                print(f"   Site {i}: ({site.get('center_x', 0):.2f}, {site.get('center_y', 0):.2f}, {site.get('center_z', 0):.2f}) - Box: {site.get('size_x', 20):.1f}×{site.get('size_y', 20):.1f}×{site.get('size_z', 20):.1f} Å")
            
            if len(all_sites) > 5:
                print(f"   ... and {len(all_sites) - 5} more sites found")
            
            # Show which site was used for docking
            primary = pocket_info.get("primary", {})
            if primary:
                print(f"\n   🎯 DOCKING PERFORMED AT:")
                print(f"   Site 1: ({primary.get('center_x', 0):.2f}, {primary.get('center_y', 0):.2f}, {primary.get('center_z', 0):.2f}) - Box: {primary.get('size_x', 20):.1f}×{primary.get('size_y', 20):.1f}×{primary.get('size_z', 20):.1f} Å")
        
        print("="*100)

        # Save the results to a CSV file if a path is provided
        if output_csv_path is None:
            # Default to saving next to the input file
            output_csv_path = pdbqt_file.parent / f"{pdbqt_file.stem}_scores.csv"
        
        df.to_csv(output_csv_path, index=False)
        print(f"\n[INFO] ✅ Results table saved to: {output_csv_path}")

        return df

    except Exception as e:
        print(f"❌ ERROR: Failed to parse results file. Reason: {e}")
        return pd.DataFrame()


# This block allows you to run this file directly for testing
if __name__ == "__main__":
    print("--- Running result parser test ---")
    
    # This test assumes you have a docked result file from the docking_engine step
    docked_file = input("Enter path to a docked PDBQT file (e.g., data/docking_results/1A3N_aspirin_docked.pdbqt): ").strip()
    
    if not docked_file:
        print("No file provided. Exiting.")
    elif not Path(docked_file).exists():
        print(f"❌ ERROR: File not found: {docked_file}")
    else:
        parse_vina_results(docked_file)
        
    print("--- Test complete ---")
