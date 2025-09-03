# backend/modules/pocket_identifier.py
# Enhanced pocket identification with configurable detail level and pipeline integration

import os
import sys
import subprocess
import tempfile
import json
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Union
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class PocketIdentifier:
    def __init__(self, detailed: bool = True):
        """
        Initialize pocket identifier with configurable detail level
        
        Args:
            detailed: If True, uses comprehensive analysis (fpocket, P2Rank)
                     If False, uses faster methods (template-based, geometric)
        """
        self.detailed = detailed
        self.methods = ["fpocket", "p2rank", "template_based"] if detailed else ["template_based", "geometric"]
        
    def _run_fpocket(self, protein_pdb: str) -> List[Dict]:
        """Run fpocket cavity detection (detailed mode only)"""
        if not self.detailed:
            return []
            
        try:
            with tempfile.TemporaryDirectory() as temp_dir:
                temp_protein = Path(temp_dir) / "protein.pdb"
                temp_protein.write_text(Path(protein_pdb).read_text())
                
                # Run fpocket
                cmd = ["fpocket", "-f", str(temp_protein)]
                result = subprocess.run(cmd, capture_output=True, text=True, cwd=temp_dir)
                
                if result.returncode != 0:
                    logger.warning("fpocket failed, trying alternative detection")
                    return []
                
                # Parse fpocket output
                pockets = []
                pocket_info_file = Path(temp_dir) / f"{temp_protein.stem}_info.txt"
                
                if pocket_info_file.exists():
                    pockets = self._parse_fpocket_output(pocket_info_file)
                
                return pockets[:10]  # Return top 10 for later filtering
                
        except Exception as e:
            logger.warning(f"fpocket execution failed: {e}")
            return []

    def _parse_fpocket_output(self, info_file: Path) -> List[Dict]:
        """Parse fpocket info file to extract pocket information"""
        pockets = []
        try:
            with open(info_file) as f:
                lines = f.readlines()
            
            current_pocket = {}
            for line in lines:
                line = line.strip()
                if line.startswith("Pocket"):
                    if current_pocket:
                        pockets.append(current_pocket)
                    current_pocket = {
                        "method": "fpocket",
                        "confidence": "high" if len(pockets) < 3 else "medium"
                    }
                elif "Centroid" in line and current_pocket:
                    # Parse centroid coordinates
                    coords = line.split()[-3:]
                    current_pocket.update({
                        "center_x": float(coords[0]),
                        "center_y": float(coords[1]), 
                        "center_z": float(coords[2])
                    })
                elif "Score" in line and current_pocket:
                    # Parse pocket score
                    score = float(line.split()[-1])
                    current_pocket["score"] = score
                elif "Volume" in line and current_pocket:
                    # Parse pocket volume
                    volume = float(line.split()[-1])
                    current_pocket["pocket_size"] = volume
                    # Estimate box size based on volume (cube root approximation)
                    box_size = max(10, min(30, (volume/10)**(1/3) * 8))
                    current_pocket.update({
                        "size_x": box_size,
                        "size_y": box_size,
                        "size_z": box_size
                    })
            
            if current_pocket:
                pockets.append(current_pocket)
                
        except Exception as e:
            logger.warning(f"Failed to parse fpocket output: {e}")
        
        return pockets

    def _run_p2rank(self, protein_pdb: str) -> List[Dict]:
        """Run P2Rank pocket prediction (detailed mode only)"""
        if not self.detailed:
            return []
            
        try:
            with tempfile.TemporaryDirectory() as temp_dir:
                # Try P2Rank if available
                cmd = ["p2rank", "predict", protein_pdb, "-o", temp_dir]
                result = subprocess.run(cmd, capture_output=True, text=True)
                
                if result.returncode != 0:
                    return []
                
                # Parse P2Rank output
                predictions_file = Path(temp_dir) / "predictions.csv"
                if not predictions_file.exists():
                    return []
                
                return self._parse_p2rank_output(predictions_file)
                
        except Exception as e:
            logger.warning(f"P2Rank execution failed: {e}")
            return []

    def _parse_p2rank_output(self, predictions_file: Path) -> List[Dict]:
        """Parse P2Rank predictions"""
        pockets = []
        try:
            with open(predictions_file) as f:
                lines = f.readlines()[1:]  # Skip header
            
            for i, line in enumerate(lines[:10]):  # Top 10 pockets
                parts = line.strip().split(',')
                if len(parts) >= 6:
                    pocket = {
                        "method": "p2rank",
                        "center_x": float(parts[1]),
                        "center_y": float(parts[2]),
                        "center_z": float(parts[3]),
                        "score": float(parts[4]),
                        "pocket_size": float(parts[5]),
                        "confidence": "high" if i < 2 else "medium"
                    }
                    # Estimate box size
                    box_size = max(12, min(25, pocket["pocket_size"]/20))
                    pocket.update({
                        "size_x": box_size,
                        "size_y": box_size,
                        "size_z": box_size
                    })
                    pockets.append(pocket)
                    
        except Exception as e:
            logger.warning(f"Failed to parse P2Rank output: {e}")
        
        return pockets

    def _template_based_detection(self, protein_pdb: str, ligand_pdbqt: str) -> List[Dict]:
        """Template-based pocket detection using ligand center (fast and detailed modes)"""
        try:
            # Parse ligand coordinates to find center
            ligand_coords = []
            with open(ligand_pdbqt) as f:
                for line in f:
                    if line.startswith(("ATOM", "HETATM")):
                        coords = line[30:54].split()
                        if len(coords) >= 3:
                            ligand_coords.append([float(x) for x in coords[:3]])
            
            if not ligand_coords:
                return []
            
            # Calculate ligand center
            center_x = sum(coord[0] for coord in ligand_coords) / len(ligand_coords)
            center_y = sum(coord[1] for coord in ligand_coords) / len(ligand_coords)
            center_z = sum(coord[2] for coord in ligand_coords) / len(ligand_coords)
            
            # Create multiple binding modes around the ligand center
            modes = []
            
            if self.detailed:
                # Detailed mode: More variations
                base_sizes = [
                    (20, 20, 20), (18, 22, 18), (22, 18, 20), 
                    (16, 20, 24), (24, 16, 18), (19, 21, 19),
                    (21, 19, 21), (17, 23, 17)
                ]
            else:
                # Fast mode: Fewer variations
                base_sizes = [(20, 20, 20), (18, 22, 18), (22, 18, 20)]
            
            for i, (sx, sy, sz) in enumerate(base_sizes):
                # Add slight variations to center coordinates
                offset_x = (i - 2) * 1.5
                offset_y = (i - 2) * 1.0
                offset_z = (i - 2) * 0.8
                
                mode = {
                    "method": "template_based",
                    "center_x": center_x + offset_x,
                    "center_y": center_y + offset_y,
                    "center_z": center_z + offset_z,
                    "size_x": sx,
                    "size_y": sy,
                    "size_z": sz,
                    "pocket_size": sx * sy * sz / 8,  # Approximate pocket volume
                    "score": 9.0 - i * 0.3,  # Decreasing scores
                    "confidence": "high" if i < 2 else "medium"
                }
                modes.append(mode)
            
            return modes
            
        except Exception as e:
            logger.warning(f"Template-based detection failed: {e}")
            return []

    def _geometric_fallback(self, protein_pdb: str) -> List[Dict]:
        """Geometric fallback method when other methods fail"""
        try:
            # Parse protein to find geometric center
            protein_coords = []
            with open(protein_pdb) as f:
                for line in f:
                    if line.startswith("ATOM") and line[12:16].strip() == "CA":  # Alpha carbons only
                        coords = line[30:54].split()
                        if len(coords) >= 3:
                            protein_coords.append([float(x) for x in coords[:3]])
            
            if not protein_coords:
                return []
            
            # Calculate geometric center
            center_x = sum(coord[0] for coord in protein_coords) / len(protein_coords)
            center_y = sum(coord[1] for coord in protein_coords) / len(protein_coords)
            center_z = sum(coord[2] for coord in protein_coords) / len(protein_coords)
            
            # Create multiple geometric modes
            modes = []
            
            if self.detailed:
                # Detailed mode: More geometric variations
                offsets = [
                    (0, 0, 0), (5, 0, 0), (-5, 0, 0), (0, 5, 0), (0, -5, 0),
                    (0, 0, 5), (0, 0, -5), (3, 3, 0), (-3, -3, 0)
                ]
                sizes = [
                    (20, 20, 20), (18, 22, 18), (22, 18, 20), (16, 24, 16),
                    (24, 16, 20), (19, 21, 19), (21, 19, 21), (17, 23, 17), (23, 17, 19)
                ]
            else:
                # Fast mode: Basic geometric variations
                offsets = [(0, 0, 0), (5, 0, 0), (-5, 0, 0)]
                sizes = [(20, 20, 20), (18, 22, 18), (22, 18, 20)]
            
            for i, (offset, size) in enumerate(zip(offsets, sizes)):
                mode = {
                    "method": "geometric",
                    "center_x": center_x + offset[0],
                    "center_y": center_y + offset[1],
                    "center_z": center_z + offset[2],
                    "size_x": size[0],
                    "size_y": size[1],
                    "size_z": size[2],
                    "pocket_size": size[0] * size[1] * size[2] / 10,
                    "score": 7.0 - i * 0.2,
                    "confidence": "medium" if i < 2 else "low"
                }
                modes.append(mode)
            
            return modes
            
        except Exception as e:
            logger.warning(f"Geometric fallback failed: {e}")
            return []


def identify_binding_site(protein_pdb: str, prepared_protein_pdbqt: str, ligand_pdbqt: str = None, 
                         use_validated: bool = True, return_n: int = 5, detailed: bool = True) -> Optional[Dict]:
    """
    Enhanced pocket identification with configurable detail level for production use.
    
    Args:
        protein_pdb: Path to original protein PDB file
        prepared_protein_pdbqt: Path to prepared receptor PDBQT file  
        ligand_pdbqt: Path to ligand PDBQT file (optional, for template-based detection)
        use_validated: Whether to use validated pocket detection methods
        return_n: Number of top binding modes to return (max 5)
        detailed: Control analysis depth (True = comprehensive, False = fast)
    
    Returns:
        Dictionary containing:
        - 'primary': Best binding mode with center coords, box size, pocket size
        - 'modes': List of up to N best modes with full information
        - 'total_found': Total number of pockets found before filtering
        - 'methods_used': List of detection methods attempted
        - 'analysis_level': Level of analysis performed ('detailed' or 'fast')
    """
    analysis_level = "detailed" if detailed else "fast"
    logger.info(f"Starting {analysis_level} pocket identification")
    
    identifier = PocketIdentifier(detailed=detailed)
    all_modes = []
    methods_used = []
    
    # Method 1: fpocket (only in detailed mode with validation)
    if detailed and use_validated:
        logger.info("Attempting fpocket detection")
        fpocket_modes = identifier._run_fpocket(protein_pdb)
        if fpocket_modes:
            all_modes.extend(fpocket_modes)
            methods_used.append("fpocket")
            logger.info(f"fpocket found {len(fpocket_modes)} pockets")
    
    # Method 2: P2Rank (only in detailed mode with validation)
    if detailed and use_validated:
        logger.info("Attempting P2Rank detection")
        p2rank_modes = identifier._run_p2rank(protein_pdb)
        if p2rank_modes:
            all_modes.extend(p2rank_modes)
            methods_used.append("p2rank")
            logger.info(f"P2Rank found {len(p2rank_modes)} pockets")
    
    # Method 3: Template-based (available in both modes)
    if ligand_pdbqt and Path(ligand_pdbqt).exists():
        logger.info("Attempting template-based detection")
        template_modes = identifier._template_based_detection(protein_pdb, ligand_pdbqt)
        if template_modes:
            all_modes.extend(template_modes)
            methods_used.append("template_based")
            logger.info(f"Template-based found {len(template_modes)} pockets")
    
    # Method 4: Geometric fallback (always available as last resort)
    if not all_modes or (not detailed and len(all_modes) < 3):
        logger.info("Using geometric fallback method")
        geometric_modes = identifier._geometric_fallback(protein_pdb)
        all_modes.extend(geometric_modes)
        methods_used.append("geometric")
        logger.info(f"Geometric fallback generated {len(geometric_modes)} pockets")
    
    if not all_modes:
        logger.error("No binding sites could be identified")
        return None
    
    # Sort all modes by score (descending)
    all_modes.sort(key=lambda x: x.get("score", 0), reverse=True)
    
    # Select best N modes (max 5)
    best_modes = all_modes[:min(return_n, 5)]
    
    # Ensure all modes have required fields
    for mode in best_modes:
        if "pocket_size" not in mode:
            mode["pocket_size"] = mode.get("size_x", 20) * mode.get("size_y", 20) * mode.get("size_z", 20) / 8
        if "confidence" not in mode:
            mode["confidence"] = "medium"
        # Add analysis level to each mode
        mode["analysis_level"] = analysis_level
    
    result = {
        "primary": best_modes[0],
        "modes": best_modes,
        "total_found": len(all_modes),
        "methods_used": methods_used,
        "analysis_level": analysis_level,
        "performance_mode": "comprehensive" if detailed else "optimized"
    }
    
    logger.info(f"{analysis_level.title()} pocket identification completed: {len(best_modes)} modes selected from {len(all_modes)} total")
    logger.info(f"Primary pocket: {best_modes[0]['method']} method, score: {best_modes[0].get('score', 'N/A')}")
    
    return result


# Backward compatibility function
def detect_binding_pockets(*args, **kwargs):
    """Backward compatibility wrapper"""
    return identify_binding_site(*args, **kwargs)


# Pipeline integration helper
def get_pocket_analysis_config(pipeline_mode: str = "production") -> Dict:
    """
    Get recommended pocket analysis configuration based on pipeline mode
    
    Args:
        pipeline_mode: 'development', 'production', 'fast', or 'comprehensive'
    
    Returns:
        Dictionary with recommended settings
    """
    configs = {
        "development": {
            "detailed": True,
            "use_validated": True,
            "return_n": 5
        },
        "production": {
            "detailed": True,
            "use_validated": True,
            "return_n": 3
        },
        "fast": {
            "detailed": False,
            "use_validated": False,
            "return_n": 3
        },
        "comprehensive": {
            "detailed": True,
            "use_validated": True,
            "return_n": 5
        }
    }
    
    return configs.get(pipeline_mode, configs["production"])


if __name__ == "__main__":
    # Test the module with different modes
    print("Testing pocket identifier in different modes...")
    
    # Fast mode test
    print("\n--- FAST MODE TEST ---")
    config = get_pocket_analysis_config("fast")
    print(f"Fast mode config: {config}")
    
    # Production mode test
    print("\n--- PRODUCTION MODE TEST ---")
    config = get_pocket_analysis_config("production")
    print(f"Production mode config: {config}")
    
    # Comprehensive mode test
    print("\n--- COMPREHENSIVE MODE TEST ---")
    config = get_pocket_analysis_config("comprehensive")
    print(f"Comprehensive mode config: {config}")
    
    # Example usage
    if len(sys.argv) >= 4:
        protein_pdb = sys.argv[1]
        prepared_pdbqt = sys.argv[2]
        ligand_pdbqt = sys.argv[3]
        mode = sys.argv[4] if len(sys.argv) > 4 else "production"
        
        config = get_pocket_analysis_config(mode)
        result = identify_binding_site(
            protein_pdb, prepared_pdbqt, ligand_pdbqt, **config
        )
        
        if result:
            print(f"\n--- RESULTS ({result['analysis_level'].upper()} MODE) ---")
            print(json.dumps(result, indent=2, default=str))
    else:
        print("\nUsage: python pocket_identifier.py <protein.pdb> <prepared.pdbqt> <ligand.pdbqt> [mode]")
        print("Modes: fast, production, comprehensive, development")
