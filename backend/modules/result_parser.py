# backend/modules/result_parser.py

import pandas as pd
from pathlib import Path

def parse_vina_results(docked_pdbqt_path: str, output_csv_path: str = None):
    """
    Parses a Vina output PDBQT file to extract binding affinity scores for each pose.

    Args:
        docked_pdbqt_path (str): The path to the docked PDBQT file from Vina.
        output_csv_path (str, optional): Path to save the results as a CSV file. 
                                         Defaults to the same directory as the input file.

    Returns:
        pandas.DataFrame: A DataFrame with the pose number and binding affinity,
                          or an empty DataFrame if parsing fails.
    """
    pdbqt_file = Path(docked_pdbqt_path)
    if not pdbqt_file.exists():
        print(f"❌ ERROR: Docked file not found at {pdbqt_file}")
        return pd.DataFrame()

    results = []
    try:
        print(f"[INFO] Parsing results from {pdbqt_file.name}...")
        with open(pdbqt_file, 'r') as f:
            lines = f.readlines()
        
        for line in lines:
            # Vina results are in lines like: "REMARK VINA RESULT:  -7.5      0.000    0.000"
            if line.startswith("REMARK VINA RESULT:"):
                parts = line.split()
                if len(parts) >= 4:
                    affinity = float(parts[3])
                    results.append(affinity)

        if not results:
            print(f"[WARN] No binding affinity scores found in {pdbqt_file.name}")
            return pd.DataFrame()

        # Create a DataFrame from the extracted scores
        df = pd.DataFrame({
            "Pose": range(1, len(results) + 1),
            "Binding Affinity (kcal/mol)": results
        })
        
        print("\n[INFO] ✅ Successfully parsed results:")
        print(df.to_string(index=False))

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