# backend/modules/template_pocket.py

from dataclasses import dataclass
from typing import List, Optional, Tuple
from pathlib import Path
import numpy as np

@dataclass
class TemplateHit:
    pdb_id: str
    ligand_id: str
    seq_identity: float
    lig_sim: float
    transform: np.ndarray  # 4x4 homogenous
    pocket_center: Tuple[float, float, float]
    pocket_size: Tuple[float, float, float]
    pocket_rmsd: Optional[float] = None

class TemplatePocketAugmentor:
    def __init__(self, index_dir: str = "data/template_index",
                 min_seqid: float = 0.35, min_ligsim: float = 0.4):
        self.index_dir = Path(index_dir)
        self.min_seqid = min_seqid
        self.min_ligsim = min_ligsim
        # lazy load small indices (ligand fingerprints, seq index, pocket centers)

    def search_templates(self, protein_fasta: str, ligand_smiles: str) -> List[TemplateHit]:
        # 1) ligand FP search (FP2/Morgan) => shortlist
        # 2) protein seq search (MMseqs2/BLAST DB-lite) => intersect
        # 3) optional pocket RMSD/pocket residue overlap filter if cached
        # Return top-k TemplateHit objects
        return []

    def transfer_pocket(self, hit: TemplateHit) -> Tuple[Tuple[float,float,float], Tuple[float,float,float], float]:
        # Apply hit.transform to template pocket center/box to map into query coords
        cx, cy, cz = hit.pocket_center
        size = hit.pocket_size
        # Confidence from seqid + ligsim (+ pocket_rmsd if present)
        conf = min(0.95, 0.55 + 0.25*hit.seq_identity + 0.25*hit.lig_sim)
        return (cx, cy, cz), size, conf

    def augment(self, protein_fasta: str, ligand_smiles: str, max_boxes=2):
        hits = self.search_templates(protein_fasta, ligand_smiles)
        boxes = []
        for h in hits[:max_boxes]:
            (cx,cy,cz), (sx,sy,sz), conf = self.transfer_pocket(h)
            boxes.append((cx,cy,cz,(sx,sy,sz), conf, "template"))
        return boxes
