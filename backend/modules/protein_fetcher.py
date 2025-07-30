import os
import requests
import time
from pathlib import Path
from datetime import datetime
import platform
import subprocess
import webbrowser


class ProteinFetcher:
    def __init__(self, download_dir="downloads/proteins"):
        self.download_dir = Path(download_dir)
        self.download_dir.mkdir(parents=True, exist_ok=True)

    def _open_file(self, path):
        try:
            if platform.system() == "Windows":
                os.startfile(path)
            elif platform.system() == "Darwin":  # macOS
                subprocess.run(["open", path])
            else:  # Linux/Unix
                subprocess.run(["xdg-open", path])
        except Exception as e:
            print(f"[WARNING] Couldn't open the file automatically: {e}")

    def _open_folder(self, path):
        try:
            if platform.system() == "Windows":
                os.startfile(path)
            elif platform.system() == "Darwin":  # macOS
                subprocess.run(["open", path])
            else:  # Linux/Unix
                subprocess.run(["xdg-open", path])
            print(f"[INFO] Opened download directory: {path}")
        except Exception as e:
            print(f"[WARNING] Couldn't open the download directory automatically: {e}")

    def _generate_unique_filename(self, base_name):
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        return f"{base_name}_{timestamp}.pdb"

    def is_online(self):
        try:
            requests.get("https://www.google.com", timeout=5)
            return True
        except:
            return False

    def fetch_from_pdb(self, pdb_id):
        pdb_id = pdb_id.lower()
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        response = requests.get(url)

        if response.status_code == 200:
            filename = self._generate_unique_filename(pdb_id)
            file_path = self.download_dir / filename
            with open(file_path, "w") as f:
                f.write(response.text)
            print(f"[INFO] ✅ Downloaded {pdb_id}.pdb from RCSB to {file_path}")
            self._open_file(str(file_path))
            return str(file_path)
        else:
            raise Exception(f"[ERROR] ❌ RCSB fetch failed for {pdb_id}. Status code: {response.status_code}")

    def fetch_from_alphafold(self, uniprot_id):
        uniprot_id = uniprot_id.upper()
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
        response = requests.get(url)

        if response.status_code == 200:
            filename = self._generate_unique_filename(f"{uniprot_id}_AF")
            file_path = self.download_dir / filename
            with open(file_path, "w") as f:
                f.write(response.text)
            print(f"[INFO] ✅ Downloaded AlphaFold model for {uniprot_id} to {file_path}")
            self._open_file(str(file_path))
            return str(file_path)
        else:
            raise Exception(f"[ERROR] ❌ AlphaFold fetch failed for {uniprot_id}. Status code: {response.status_code}")

    def generate_from_sequence(self, sequence_placeholder=None): # Added placeholder to match call signature if needed elsewhere
        if not self.is_online():
            print("⚠️ You are offline. Please connect to the internet to use AlphaFold2 Colab.")
            return None

        colab_url = "https://colab.research.google.com/drive/1gAaNunQVQQvbzyBdgj3VpnH7MAa6bA9K"

        print("\n[INFO] To build a protein from a FASTA sequence:")
        print("1. Go to the AlphaFold2 Colab notebook.")
        print("2. Paste your FASTA sequence into the appropriate field.")
        print("3. Run the notebook to generate the protein structure.")
        print(f"\n🔗 Open Colab Notebook: {colab_url}")

        webbrowser.open(colab_url)

        filename_prefix = "generated_from_sequence"
        print(f"\n💾 After generating and downloading your PDB, please save it as '{filename_prefix}_<timestamp>.pdb'")
        print(f"📁 And place it manually in: {self.download_dir.resolve()}")

        input("\n🕐 Press Enter once you've saved the file manually in the download directory...")

        # You still need to check if the user actually placed a file,
        # but you don't know its exact name without the timestamp.
        # A more robust check here would be to list files in the directory
        # that match the pattern "generated_from_sequence_*.pdb" and pick the newest.
        # For simplicity, let's keep the existing check assuming the user saves with the suggested name.
        
        # This part assumes the user will manually create a file with the correct naming convention
        # You could also prompt the user for the actual filename they saved it as.
        dummy_filename = self._generate_unique_filename(filename_prefix)
        expected_path_example = self.download_dir / dummy_filename

        if any(f.name.startswith(f"{filename_prefix}_") and f.suffix == ".pdb" for f in self.download_dir.iterdir()):
             print(f"✅ It seems a generated structure (e.g., '{filename_prefix}_*.pdb') is in {self.download_dir}")
             # You might want to ask the user which file it is, or open the folder directly
             self._open_folder(str(self.download_dir.resolve()))
             return "Path to the manually saved file (needs user input or more advanced search)"
        else:
            print("❌ No new generated file found matching the expected pattern. Please ensure it was saved correctly.")
            self._open_folder(str(self.download_dir.resolve())) # Still open folder for user
            return None


if __name__ == "__main__":
    fetcher = ProteinFetcher()

    user_input = input("🔎 Enter PDB ID, UniProt ID, or FASTA sequence (or 'fasta' to generate from sequence): ").strip()

    if len(user_input) == 4 and user_input.isalnum():
        try:
            fetcher.fetch_from_pdb(user_input)
        except Exception as e:
            print(e)
            print("[INFO] Trying AlphaFold fallback...")
            try:
                fetcher.fetch_from_alphafold(user_input)
            except Exception as e2:
                print(e2)
    elif user_input.lower() == 'fasta': # Changed condition to 'fasta' keyword
        fetcher.generate_from_sequence()
    else:
        print("[ERROR] ❌ Invalid input. Please provide a valid 4-character PDB ID, UniProt ID, or type 'fasta' to generate from sequence.")