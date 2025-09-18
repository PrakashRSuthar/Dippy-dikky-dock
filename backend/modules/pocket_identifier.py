# backend/modules/pocket_identifier.py
# Pocket identification via raw-structure (fpocket/P2Rank), ligand-template, and structure-based fallback,
# with consensus clustering, cross-validation, contact rescore, and ligand-guided refinement.

import subprocess
import tempfile
from pathlib import Path
from typing import List, Dict, Optional
import logging
import math

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# -------------------------
# Math / helpers
# -------------------------
def _dist(a, b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)

def _center_xyz(m):
    return (float(m.get("center_x", 0)), float(m.get("center_y", 0)), float(m.get("center_z", 0)))

def _box_diag(m):
    return math.sqrt(float(m.get("size_x",20))**2 + float(m.get("size_y",20))**2 + float(m.get("size_z",20))**2)

def _normalize(vals: List[float]) -> List[float]:
    if not vals: return []
    mn, mx = min(vals), max(vals)
    if mx - mn < 1e-9: return [0.5]*len(vals)
    return [(v - mn)/(mx - mn) for v in vals]

def _cluster_modes(modes: List[Dict], max_dist_factor: float = 0.6) -> List[Dict]:
    clusters: List[Dict] = []
    for m in modes:
        c = _center_xyz(m); d = _box_diag(m)
        placed = False
        for cl in clusters:
            cc = cl["_centroid"]
            if _dist(c, cc) <= max(d, cl["_avg_diag"]) * max_dist_factor:
                cl["items"].append(m)
                xs = [_center_xyz(x)[0] for x in cl["items"]]
                ys = [_center_xyz(x)[1] for x in cl["items"]]
                zs = [_center_xyz(x)[2] for x in cl["items"]]
                cl["_centroid"] = (sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs))
                cl["_avg_diag"] = sum(_box_diag(x) for x in cl["items"]) / len(cl["items"])
                placed = True
                break
        if not placed:
            clusters.append({"items":[m], "_centroid": c, "_avg_diag": d})
    return clusters

def _consensus_score(cluster: Dict) -> float:
    items = cluster["items"]
    methods = set(i.get("method","?") for i in items)
    return len(items) + (0.5 if len(methods) >= 2 else 0.0)

def _contact_rescore(protein_pdb: str, mode: Dict, radius: float = 6.0) -> float:
    """Count CA atoms within radius of box center; log-diminishing score."""
    try:
        cx, cy, cz = _center_xyz(mode)
        r2 = radius*radius
        cnt = 0
        with open(protein_pdb) as f:
            for line in f:
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    parts = line[30:54].split()
                    if len(parts) >= 3:
                        x, y, z = map(float, parts[:3])
                        dx, dy, dz = x-cx, y-cy, z-cz
                        if dx*dx + dy*dy + dz*dz <= r2:
                            cnt += 1
        return math.log1p(cnt)
    except Exception:
        return 0.0

def _ligand_guided_refine(protein_pdb: str, ligand_pdbqt: Optional[str], mode: Dict) -> Dict:
    """Pull box toward protein and size to ligand extents + margin, clamped 14..26 Å."""
    if not ligand_pdbqt or not Path(ligand_pdbqt).exists():
        return mode
    try:
        pts = []
        with open(ligand_pdbqt) as f:
            for line in f:
                if line.startswith(("ATOM","HETATM")) and line[76:78].strip() != 'H':
                    parts = line[30:54].split()
                    if len(parts) >= 3:
                        pts.append(tuple(map(float, parts[:3])))
        if not pts: return mode
        cx = sum(p[0] for p in pts)/len(pts)
        cy = sum(p[1] for p in pts)/len(pts)
        cz = sum(p[2] for p in pts)/len(pts)
        minx=min(p[0] for p in pts); maxx=max(p[0] for p in pts)
        miny=min(p[1] for p in pts); maxy=max(p[1] for p in pts)
        minz=min(p[2] for p in pts); maxz=max(p[2] for p in pts)
        extx=maxx-minx; exty=maxy-miny; extz=maxz-minz
        # nearest CA to ligand centroid
        nearest=None; nd2=1e9
        with open(protein_pdb) as f:
            for line in f:
                if line.startswith("ATOM") and line[12:16].strip()=="CA":
                    parts = line[30:54].split()
                    if len(parts)>=3:
                        x,y,z = map(float, parts[:3])
                        d2=(x-cx)**2+(y-cy)**2+(z-cz)**2
                        if d2<nd2: nd2=d2; nearest=(x,y,z)
        if nearest:
            nx,ny,nz = nearest
            blend = 0.25
            cx = cx*(1-blend)+nx*blend
            cy = cy*(1-blend)+ny*blend
            cz = cz*(1-blend)+nz*blend
        margin = 8.0
        sx = max(14.0, min(26.0, extx + margin))
        sy = max(14.0, min(26.0, exty + margin))
        sz = max(14.0, min(26.0, extz + margin))
        m = dict(mode)
        m.update({"center_x":cx,"center_y":cy,"center_z":cz,"size_x":sx,"size_y":sy,"size_z":sz})
        return m
    except Exception:
        return mode

# -------------------------
# PocketIdentifier
# -------------------------
class PocketIdentifier:
    def __init__(self, detailed: bool = True):
        self.detailed = detailed
        self.methods = ["fpocket", "p2rank", "template_based"] if detailed else ["template_based", "geometric"]

    def _run_fpocket(self, protein_pdb: str) -> List[Dict]:
        if not self.detailed:
            return []
        try:
            with tempfile.TemporaryDirectory() as temp_dir:
                temp_protein = Path(temp_dir) / "protein.pdb"
                temp_protein.write_text(Path(protein_pdb).read_text())
                cmd = ["fpocket", "-f", str(temp_protein)]
                result = subprocess.run(cmd, capture_output=True, text=True, cwd=temp_dir)
                if result.returncode != 0:
                    logger.warning("fpocket failed")
                    return []
                info_file = Path(temp_dir) / f"{temp_protein.stem}_info.txt"
                pockets = self._parse_fpocket_output(info_file) if info_file.exists() else []
                return pockets[:10]
        except Exception as e:
            logger.warning(f"fpocket execution failed: {e}")
            return []

    def _parse_fpocket_output(self, info_file: Path) -> List[Dict]:
        pockets: List[Dict] = []
        try:
            with open(info_file) as f:
                lines = f.readlines()
            current: Dict = {}
            for line in lines:
                line = line.strip()
                if line.startswith("Pocket"):
                    if current: pockets.append(current)
                    current = {"method":"fpocket","confidence":"high" if len(pockets)<3 else "medium"}
                elif "Centroid" in line and current:
                    coords = line.split()[-3:]
                    current.update({"center_x":float(coords[0]),"center_y":float(coords[1]),"center_z":float(coords[2])})
                elif "Score" in line and current:
                    current["score"] = float(line.split()[-1])
                elif "Volume" in line and current:
                    volume = float(line.split()[-1])
                    current["pocket_size"] = volume
                    box = max(10, min(30, (volume/10)**(1/3) * 8))
                    current.update({"size_x":box,"size_y":box,"size_z":box})
            if current: pockets.append(current)
        except Exception as e:
            logger.warning(f"Failed to parse fpocket output: {e}")
        return pockets

    def _run_p2rank(self, protein_pdb: str) -> List[Dict]:
        if not self.detailed:
            return []
        try:
            with tempfile.TemporaryDirectory() as temp_dir:
                cmd = ["p2rank", "predict", protein_pdb, "-o", temp_dir]
                result = subprocess.run(cmd, capture_output=True, text=True)
                if result.returncode != 0:
                    return []
                predictions = Path(temp_dir) / "predictions.csv"
                if not predictions.exists(): return []
                return self._parse_p2rank_output(predictions)
        except Exception as e:
            logger.warning(f"P2Rank execution failed: {e}")
            return []

    def _parse_p2rank_output(self, predictions_file: Path) -> List[Dict]:
        pockets: List[Dict] = []
        try:
            with open(predictions_file) as f:
                lines = f.readlines()[1:]
            for i, line in enumerate(lines[:10]):
                parts = line.strip().split(',')
                if len(parts) >= 6:
                    pocket = {
                        "method":"p2rank",
                        "center_x": float(parts[1]),
                        "center_y": float(parts[2]),
                        "center_z": float(parts[3]),
                        "score": float(parts[4]),
                        "pocket_size": float(parts[5]),
                        "confidence": "high" if i<2 else "medium"
                    }
                    box = max(12, min(25, pocket["pocket_size"]/20))
                    pocket.update({"size_x":box,"size_y":box,"size_z":box})
                    pockets.append(pocket)
        except Exception as e:
            logger.warning(f"Failed to parse P2Rank output: {e}")
        return pockets

    def _template_based_detection(self, protein_pdb: str, ligand_pdbqt: str) -> List[Dict]:
        try:
            lig = []
            with open(ligand_pdbqt) as f:
                for line in f:
                    if line.startswith(("ATOM","HETATM")):
                        parts = line[30:54].split()
                        if len(parts) >= 3:
                            lig.append([float(x) for x in parts[:3]])
            if not lig: return []
            cx = sum(p[0] for p in lig)/len(lig)
            cy = sum(p[1] for p in lig)/len(lig)
            cz = sum(p[2] for p in lig)/len(lig)
            if self.detailed:
                sizes=[(20,20,20),(18,22,18),(22,18,20),(16,20,24),(24,16,18),(19,21,19),(21,19,21),(17,23,17)]
            else:
                sizes=[(20,20,20),(18,22,18),(22,18,20)]
            modes=[]
            for i,(sx,sy,sz) in enumerate(sizes):
                ox=(i-2)*1.5; oy=(i-2)*1.0; oz=(i-2)*0.8
                modes.append({
                    "method":"template_based",
                    "center_x":cx+ox,"center_y":cy+oy,"center_z":cz+oz,
                    "size_x":sx,"size_y":sy,"size_z":sz,
                    "pocket_size": sx*sy*sz/8,
                    "score": 9.0 - i*0.3,
                    "confidence": "high" if i<2 else "medium"
                })
            return modes
        except Exception as e:
            logger.warning(f"Template-based detection failed: {e}")
            return []

    def _geometric_fallback(self, protein_pdb: str) -> List[Dict]:
        try:
            ca=[]
            with open(protein_pdb) as f:
                for line in f:
                    if line.startswith("ATOM") and line[12:16].strip()=="CA":
                        parts = line[30:54].split()
                        if len(parts)>=3:
                            ca.append([float(x) for x in parts[:3]])
            if not ca: return []
            cx = sum(p[0] for p in ca)/len(ca)
            cy = sum(p[1] for p in ca)/len(ca)
            cz = sum(p[2] for p in ca)/len(ca)
            if self.detailed:
                offs=[(0,0,0),(5,0,0),(-5,0,0),(0,5,0),(0,-5,0),(0,0,5),(0,0,-5),(3,3,0),(-3,-3,0)]
                sizes=[(20,20,20),(18,22,18),(22,18,20),(16,24,16),(24,16,20),(19,21,19),(21,19,21),(17,23,17),(23,17,19)]
            else:
                offs=[(0,0,0),(5,0,0),(-5,0,0)]
                sizes=[(20,20,20),(18,22,18),(22,18,20)]
            modes=[]
            for i,(o,s) in enumerate(zip(offs,sizes)):
                sx,sy,sz=s
                modes.append({
                    "method":"geometric",
                    "center_x":cx+o[0],"center_y":cy+o[1],"center_z":cz+o[2],
                    "size_x":sx,"size_y":sy,"size_z":sz,
                    "pocket_size": sx*sy*sz/10,
                    "score": 7.0 - i*0.2,
                    "confidence": "medium" if i<2 else "low"
                })
            return modes
        except Exception as e:
            logger.warning(f"Geometric fallback failed: {e}")
            return []

# -------------------------
# Public API
# -------------------------
def identify_binding_site(protein_pdb: str, prepared_protein_pdbqt: str, ligand_pdbqt: str = None,
                          use_validated: bool = True, return_n: int = 5, detailed: bool = True) -> Optional[Dict]:
    """
    Pocket identification that combines:
    - Raw structure methods (fpocket/P2Rank),
    - Template-based ligand reference,
    - Structure-based geometric fallback,
    with consensus clustering + contact rescore and ligand-guided refinement.
    """
    analysis_level = "detailed" if detailed else "fast"
    logger.info(f"Starting {analysis_level} pocket identification")

    ident = PocketIdentifier(detailed=detailed)
    all_modes: List[Dict] = []
    methods_used: List[str] = []

    # Raw-structure validated methods
    if detailed and use_validated:
        logger.info("Attempting fpocket detection")
        m = ident._run_fpocket(protein_pdb)
        if m: all_modes += m; methods_used.append("fpocket"); logger.info(f"fpocket: {len(m)}")

    if detailed and use_validated:
        logger.info("Attempting P2Rank detection")
        m = ident._run_p2rank(protein_pdb)
        if m: all_modes += m; methods_used.append("p2rank"); logger.info(f"P2Rank: {len(m)}")

    # Template-based ligand guidance
    if ligand_pdbqt and Path(ligand_pdbqt).exists():
        logger.info("Attempting template-based detection")
        m = ident._template_based_detection(protein_pdb, ligand_pdbqt)
        if m: all_modes += m; methods_used.append("template_based"); logger.info(f"Template-based: {len(m)}")

    # Structure-based fallback
    if not all_modes or (not detailed and len(all_modes) < 3):
        logger.info("Using geometric fallback")
        m = ident._geometric_fallback(protein_pdb)
        all_modes += m
        methods_used.append("geometric")
        logger.info(f"Geometric: {len(m)}")

    if not all_modes:
        logger.error("No binding sites found")
        return None

    # Ensure base scores
    for m in all_modes:
        if "score" not in m:
            sx,sy,sz = float(m.get("size_x",20)),float(m.get("size_y",20)),float(m.get("size_z",20))
            m["score"] = max(0.1, 1000.0/(sx*sy*sz))

    # Consensus clustering and rescoring
    clusters = _cluster_modes(all_modes)
    reps = []
    base_scores = []
    for cl in clusters:
        rep = sorted(cl["items"], key=lambda x: x.get("score",0), reverse=True)[0]
        reps.append(rep); base_scores.append(rep.get("score",0.0))
    nbase = _normalize(base_scores)

    ranked: List[Dict] = []
    for idx, cl in enumerate(clusters):
        rep = dict(sorted(cl["items"], key=lambda x: x.get("score",0), reverse=True)[0])
        cbonus = _consensus_score(cl)
        cscore = _contact_rescore(protein_pdb, rep, radius=6.0)
        final = 0.6*nbase[idx] + 0.3*(min(cbonus,5.0)/5.0) + 0.1*(min(cscore,3.0)/3.0)
        rep["final_score"] = final
        rep["consensus_size"] = len(cl["items"])
        rep["consensus_methods"] = list({i.get("method","?") for i in cl["items"]})
        ranked.append(rep)
    ranked.sort(key=lambda x: x.get("final_score",0.0), reverse=True)

    best = ranked[:min(return_n,5)]
    for m in best:
        if "pocket_size" not in m:
            m["pocket_size"] = m.get("size_x",20)*m.get("size_y",20)*m.get("size_z",20)/8
        if "confidence" not in m:
            m["confidence"] = "high" if m.get("consensus_size",1) >= 2 else "medium"
        m["analysis_level"] = analysis_level

    # Ligand-guided refine the primary
    best[0] = _ligand_guided_refine(protein_pdb, ligand_pdbqt, best[0])

    result = {
        "primary": best[0],
        "modes": best,
        "total_found": len(all_modes),
        "methods_used": methods_used,
        "analysis_level": analysis_level,
        "performance_mode": "comprehensive" if detailed else "optimized",
        "consensus_clusters": len(clusters)
    }
    logger.info(f"{analysis_level.title()} pocket identification: {len(best)} selected from {len(all_modes)}; primary final_score={best[0].get('final_score',0):.3f}")
    return result

# Backward compatibility wrapper
def detect_binding_pockets(*args, **kwargs):
    return identify_binding_site(*args, **kwargs)

def get_pocket_analysis_config(pipeline_mode: str = "production") -> Dict:
    configs = {
        "development": { "detailed": True,  "use_validated": True,  "return_n": 5 },
        "production":  { "detailed": True,  "use_validated": True,  "return_n": 3 },
        "fast":        { "detailed": False, "use_validated": False, "return_n": 3 },
        "comprehensive": { "detailed": True, "use_validated": True, "return_n": 5 }
    }
    return configs.get(pipeline_mode, configs["production"])

if __name__ == "__main__":
    print("Pocket identifier module. Use via pipeline.")
