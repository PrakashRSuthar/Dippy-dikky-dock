# pocket_identifier.py
"""
Multi-Method Consensus Pocket Identification Engine
Paper Section III-D/E — Eq. 1: final_score = 0.60*norm_base + 0.30*consensus + 0.10*contact
"""
import subprocess, tempfile, math, logging
from pathlib import Path
from typing import List, Dict, Optional, Tuple

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

W_BASE, W_CONSENSUS, W_CONTACT = 0.60, 0.30, 0.10
BOX_MIN_AXIS, BOX_MAX_AXIS, CLAMP_BUFFER, CONTACT_RADIUS = 14.0, 26.0, 1.5, 6.0
CLUSTER_MERGE_FACTOR = 0.60

# Ligand-aware box sizing: per-Axis maximum scales with ligand extent.
# For small ligands (≤ 5 heavy atoms) we keep the legacy 26 Å cap.
# For medium ligands (6-30 heavy atoms) we scale up to 36 Å.
# For large ligands (> 30 heavy atoms) we allow up to 44 Å.
# This is the key change that lets peptides / macrocycles dock correctly.
LIGAND_BOX_MIN = 14.0
LIGAND_BOX_SMALL_MAX = 26.0   # ≤5 heavy atoms — legacy behaviour
LIGAND_BOX_MEDIUM_MAX = 36.0  # 6-30 heavy atoms
LIGAND_BOX_LARGE_MAX = 44.0   # >30 heavy atoms
LIGAND_ATOM_THRESHOLDS = (5, 30)  # (small, medium) boundary


def _ligand_heavy_atom_count(ligand_path: str) -> int:
    """Count heavy (non-H) atoms in a PDBQT ligand file."""
    count = 0
    try:
        with open(ligand_path) as fh:
            for line in fh:
                if line.startswith(("ATOM", "HETATM")):
                    # PDBQT element column is 77-78, but not always reliable.
                    # Fallback: any non-H atom type.
                    atype = line[76:78].strip() if len(line) > 77 else ""
                    if atype and atype not in ("H", "HD", "HS"):
                        count += 1
    except Exception:
        pass
    return count


def _ligand_axis_max(ligand_path: str) -> float:
    """Return the per-axis box maximum appropriate for the ligand size."""
    n_heavy = _ligand_heavy_atom_count(ligand_path)
    if n_heavy <= LIGAND_ATOM_THRESHOLDS[0]:
        return LIGAND_BOX_SMALL_MAX
    elif n_heavy <= LIGAND_ATOM_THRESHOLDS[1]:
        return LIGAND_BOX_MEDIUM_MAX
    else:
        return LIGAND_BOX_LARGE_MAX

def _dist(a, b): return math.sqrt(sum((x-y)**2 for x,y in zip(a,b)))
def _center_xyz(m): return (float(m.get("center_x",0)), float(m.get("center_y",0)), float(m.get("center_z",0)))
def _box_diag(m): return math.sqrt(float(m.get("size_x",20))**2+float(m.get("size_y",20))**2+float(m.get("size_z",20))**2)

def _normalize(vals):
    if not vals: return []
    mn, mx = min(vals), max(vals)
    span = mx - mn
    return [0.5]*len(vals) if span < 1e-9 else [(v-mn)/span for v in vals]

def _cluster_modes(modes, merge_factor=CLUSTER_MERGE_FACTOR):
    clusters = []
    for m in modes:
        c = _center_xyz(m); d = _box_diag(m); placed = False
        for cl in clusters:
            if _dist(c, cl["_centroid"]) <= max(d, cl["_avg_diag"]) * merge_factor:
                cl["items"].append(m)
                n = len(cl["items"])
                cl["_centroid"] = tuple(sum(_center_xyz(x)[i] for x in cl["items"])/n for i in range(3))
                cl["_avg_diag"] = sum(_box_diag(x) for x in cl["items"])/n
                placed = True; break
        if not placed:
            clusters.append({"items":[m], "_centroid":c, "_avg_diag":d})
    return clusters

def _consensus_support(cl):
    methods = {i.get("method","?") for i in cl["items"]}
    return len(cl["items"]) + (0.5 if len(methods) >= 2 else 0.0)

def _contact_rescore(protein_pdb, mode, radius=CONTACT_RADIUS):
    try:
        cx, cy, cz = _center_xyz(mode); r2 = radius*radius; count = 0
        with open(protein_pdb) as fh:
            for line in fh:
                if line.startswith("ATOM") and line[12:16].strip()=="CA":
                    parts = line[30:54].split()
                    if len(parts)>=3:
                        x,y,z = map(float,parts[:3])
                        if (x-cx)**2+(y-cy)**2+(z-cz)**2 <= r2: count += 1
        return math.log1p(count)
    except: return 0.0

def clamp_box_to_protein(protein_pdb, mode, ligand_pdbqt=None):
    """Algorithm 1 — Geometry-Aware Bounding Box Clamping (paper Section III-E)
    
    Updated to use ligand-aware maximum box size: the per-axis upper bound
    is no longer fixed at 26 Å but scales with the ligand's heavy-atom count.
    This prevents truncating large/extended ligands into an insufficient search
    volume while still keeping grids compact for small drug-like molecules.
    """
    xs,ys,zs = [],[],[]
    try:
        with open(protein_pdb) as fh:
            for line in fh:
                if line.startswith("ATOM") and line[12:16].strip()=="CA":
                    parts = line[30:54].split()
                    if len(parts)>=3:
                        x,y,z = map(float,parts[:3])
                        xs.append(x); ys.append(y); zs.append(z)
    except: return mode
    if not xs: return mode
    bounds = [(min(xs),max(xs)),(min(ys),max(ys)),(min(zs),max(zs))]
    keys_c = ["center_x","center_y","center_z"]
    keys_s = ["size_x","size_y","size_z"]
    box_max = _ligand_axis_max(ligand_pdbqt) if ligand_pdbqt else BOX_MAX_AXIS
    m = dict(mode); drifted = False
    for i,(lo,hi) in enumerate(bounds):
        c = float(m.get(keys_c[i],0)); s = float(m.get(keys_s[i],20)); half = s/2.0; orig = c
        if c+half < lo: c = lo+half+CLAMP_BUFFER
        if c-half > hi: c = hi-half-CLAMP_BUFFER
        s = max(BOX_MIN_AXIS, min(box_max, s))
        if abs(c-orig) > 0.01: drifted = True
        m[keys_c[i]] = c; m[keys_s[i]] = s
    m["_clamped"] = drifted
    m["_box_max_axis"] = box_max  # record for downstream logging
    return m

def _ligand_guided_refine(protein_pdb, ligand_pdbqt, mode):
    if not ligand_pdbqt or not Path(ligand_pdbqt).exists(): return mode
    try:
        pts = []
        with open(ligand_pdbqt) as fh:
            for line in fh:
                if line.startswith(("ATOM","HETATM")) and line[76:78].strip() not in ("H","HD","HS"):
                    parts = line[30:54].split()
                    if len(parts)>=3: pts.append(tuple(map(float,parts[:3])))
        if not pts: return mode
        cx = sum(p[0] for p in pts)/len(pts)
        cy = sum(p[1] for p in pts)/len(pts)
        cz = sum(p[2] for p in pts)/len(pts)
        extx = max(p[0] for p in pts)-min(p[0] for p in pts)
        exty = max(p[1] for p in pts)-min(p[1] for p in pts)
        extz = max(p[2] for p in pts)-min(p[2] for p in pts)
        nearest=None; nd2=1e9
        with open(protein_pdb) as fh:
            for line in fh:
                if line.startswith("ATOM") and line[12:16].strip()=="CA":
                    parts=line[30:54].split()
                    if len(parts)>=3:
                        x,y,z=map(float,parts[:3]); d2=(x-cx)**2+(y-cy)**2+(z-cz)**2
                        if d2<nd2: nd2=d2; nearest=(x,y,z)
        BLEND=0.25
        if nearest:
            cx=cx*(1-BLEND)+nearest[0]*BLEND
            cy=cy*(1-BLEND)+nearest[1]*BLEND
            cz=cz*(1-BLEND)+nearest[2]*BLEND
        margin=8.0
        box_max = _ligand_axis_max(ligand_pdbqt)
        sx=max(BOX_MIN_AXIS,min(box_max,extx+margin))
        sy=max(BOX_MIN_AXIS,min(box_max,exty+margin))
        sz=max(BOX_MIN_AXIS,min(box_max,extz+margin))
        refined=dict(mode); refined.update({"center_x":cx,"center_y":cy,"center_z":cz,"size_x":sx,"size_y":sy,"size_z":sz})
        return refined
    except: return mode

class PocketIdentifier:
    def __init__(self, detailed=True): self.detailed = detailed

    def _run_fpocket(self, protein_pdb):
        if not self.detailed: return []
        try:
            with tempfile.TemporaryDirectory() as tmp:
                tp = Path(tmp)/"protein.pdb"; tp.write_text(Path(protein_pdb).read_text())
                result = subprocess.run(["fpocket","-f",str(tp)], capture_output=True, text=True, cwd=tmp)
                if result.returncode!=0: return []
                info_file = Path(tmp)/f"{tp.stem}_info.txt"
                return self._parse_fpocket(info_file)[:10] if info_file.exists() else []
        except FileNotFoundError: logger.info("fpocket not installed"); return []
        except Exception as e: logger.warning("fpocket: %s",e); return []

    def _parse_fpocket(self, info_file):
        pockets=[]; current={}
        try:
            with open(info_file) as fh:
                for line in fh:
                    line=line.strip()
                    if line.startswith("Pocket"):
                        if current and "center_x" in current: pockets.append(current)
                        current={"method":"fpocket","confidence":"high" if len(pockets)<3 else "medium"}
                    elif "Centroid" in line and current:
                        parts=line.split()[-3:]
                        try: current.update({"center_x":float(parts[0]),"center_y":float(parts[1]),"center_z":float(parts[2])})
                        except: pass
                    elif "Score" in line and current:
                        try: current["score"]=float(line.split()[-1])
                        except: pass
                    elif "Volume" in line and current:
                        try:
                            v=float(line.split()[-1]); current["pocket_size"]=v
                            box=max(BOX_MIN_AXIS,min(BOX_MAX_AXIS,(v/10)**(1/3)*8))
                            current.update({"size_x":box,"size_y":box,"size_z":box})
                        except: pass
            if current and "center_x" in current: pockets.append(current)
        except Exception as e: logger.warning("fpocket parse: %s",e)
        return pockets

    def _run_p2rank(self, protein_pdb):
        if not self.detailed: return []
        try:
            sd=Path(__file__).parent; p2d=sd.parent/"tools"/"p2rank_2.4.2"
            import platform
            exe = str(p2d/("prank.bat" if platform.system()=="Windows" else "prank")) if p2d.exists() else "p2rank"
            with tempfile.TemporaryDirectory() as tmp:
                result=subprocess.run([exe,"predict",protein_pdb,"-o",tmp],capture_output=True,text=True)
                if result.returncode!=0: return []
                pred=Path(tmp)/"predictions.csv"
                return self._parse_p2rank(pred) if pred.exists() else []
        except FileNotFoundError: logger.info("P2Rank not installed"); return []
        except Exception as e: logger.warning("P2Rank: %s",e); return []

    def _parse_p2rank(self, pred_file):
        pockets=[]
        try:
            with open(pred_file) as fh: lines=fh.readlines()[1:]
            for i,line in enumerate(lines[:10]):
                parts=line.strip().split(",")
                if len(parts)<6: continue
                try:
                    v=float(parts[5]); box=max(BOX_MIN_AXIS,min(BOX_MAX_AXIS,v/20))
                    pockets.append({"method":"p2rank","center_x":float(parts[1]),"center_y":float(parts[2]),"center_z":float(parts[3]),
                                    "score":float(parts[4]),"pocket_size":v,"size_x":box,"size_y":box,"size_z":box,
                                    "confidence":"high" if i<2 else "medium"})
                except: continue
        except Exception as e: logger.warning("P2Rank parse: %s",e)
        return pockets

    def _template_based(self, protein_pdb, ligand_pdbqt):
        try:
            pts=[]
            with open(ligand_pdbqt) as fh:
                for line in fh:
                    if line.startswith(("ATOM","HETATM")):
                        parts=line[30:54].split()
                        if len(parts)>=3: pts.append([float(x) for x in parts[:3]])
            if not pts: return []
            cx=sum(p[0] for p in pts)/len(pts); cy=sum(p[1] for p in pts)/len(pts); cz=sum(p[2] for p in pts)/len(pts)
            offsets=[(0,0,0),(2,0,0),(-2,0,0),(0,2,0),(0,-2,0),(0,0,2),(0,0,-2),(1.5,1.5,0),(-1.5,-1.5,0)] if self.detailed else [(0,0,0),(2,0,0),(-2,0,0)]
            sizes=[(20,20,20),(18,22,18),(22,18,20),(16,20,24),(24,16,18),(19,21,19),(21,19,21),(17,23,17),(23,17,19)] if self.detailed else [(20,20,20),(18,22,18),(22,18,20)]
            modes=[]
            for i,(o,s) in enumerate(zip(offsets,sizes)):
                sx,sy,sz=s
                modes.append({"method":"template_based","center_x":cx+o[0],"center_y":cy+o[1],"center_z":cz+o[2],
                               "size_x":sx,"size_y":sy,"size_z":sz,"pocket_size":sx*sy*sz/8.0,
                               "score":max(0.1,9.0-i*0.3),"confidence":"high" if i<2 else "medium"})
            return modes
        except Exception as e: logger.warning("Template: %s",e); return []

    def _geometric_fallback(self, protein_pdb):
        try:
            ca=[]
            with open(protein_pdb) as fh:
                for line in fh:
                    if line.startswith("ATOM") and line[12:16].strip()=="CA":
                        parts=line[30:54].split()
                        if len(parts)>=3: ca.append([float(x) for x in parts[:3]])
            if not ca: return []
            cx=sum(p[0] for p in ca)/len(ca); cy=sum(p[1] for p in ca)/len(ca); cz=sum(p[2] for p in ca)/len(ca)
            offsets=[(0,0,0),(5,0,0),(-5,0,0),(0,5,0),(0,-5,0),(0,0,5),(0,0,-5),(3,3,0),(-3,-3,0)]
            sizes=[(20,20,20),(18,22,18),(22,18,20),(16,24,16),(24,16,20),(19,21,19),(21,19,21),(17,23,17),(23,17,19)]
            if not self.detailed: offsets=offsets[:3]; sizes=sizes[:3]
            modes=[]
            for i,(o,s) in enumerate(zip(offsets,sizes)):
                sx,sy,sz=s
                modes.append({"method":"geometric","center_x":cx+o[0],"center_y":cy+o[1],"center_z":cz+o[2],
                               "size_x":sx,"size_y":sy,"size_z":sz,"pocket_size":sx*sy*sz/10.0,
                               "score":max(0.1,7.0-i*0.2),"confidence":"medium" if i<2 else "low"})
            return modes
        except Exception as e: logger.warning("Geometric: %s",e); return []

def identify_binding_site(protein_pdb, prepared_protein_pdbqt, ligand_pdbqt=None,
                           use_validated=True, return_n=5, detailed=True):
    level = "detailed" if detailed else "fast"
    logger.info("Pocket identification [%s]", level)
    ident = PocketIdentifier(detailed=detailed)
    all_modes=[]; methods_used=[]

    if detailed and use_validated:
        m=ident._run_fpocket(protein_pdb)
        if m: all_modes+=m; methods_used.append("fpocket"); logger.info("fpocket: %d",len(m))

    if detailed and use_validated:
        m=ident._run_p2rank(protein_pdb)
        if m: all_modes+=m; methods_used.append("p2rank"); logger.info("p2rank: %d",len(m))

    if ligand_pdbqt and Path(ligand_pdbqt).exists():
        m=ident._template_based(protein_pdb, ligand_pdbqt)
        if m: all_modes+=m; methods_used.append("template_based"); logger.info("template: %d",len(m))

    if not all_modes:
        m=ident._geometric_fallback(protein_pdb)
        all_modes+=m; methods_used.append("geometric")

    if not all_modes: return None

    for m in all_modes:
        if "score" not in m:
            sx,sy,sz = float(m.get("size_x",20)),float(m.get("size_y",20)),float(m.get("size_z",20))
            m["score"] = max(0.1,1000.0/(sx*sy*sz))

    clusters = _cluster_modes(all_modes)
    base_scores = [max(cl["items"],key=lambda x:x.get("score",0)).get("score",0.0) for cl in clusters]
    norm_bases = _normalize(base_scores)

    ranked=[]
    for idx,cl in enumerate(clusters):
        rep = dict(max(cl["items"],key=lambda x:x.get("score",0)))
        # Eq. 1 (paper Section III-D)
        final = (W_BASE*norm_bases[idx] +
                 W_CONSENSUS*min(_consensus_support(cl),5.0)/5.0 +
                 W_CONTACT*min(_contact_rescore(protein_pdb,rep),3.0)/3.0)
        rep["final_score"]=round(final,4)
        rep["consensus_size"]=len(cl["items"])
        rep["consensus_methods"]=sorted({i.get("method","?") for i in cl["items"]})
        if "pocket_size" not in rep:
            rep["pocket_size"]=float(rep.get("size_x",20))*float(rep.get("size_y",20))*float(rep.get("size_z",20))/8.0
        if "confidence" not in rep:
            rep["confidence"]="high" if len(cl["items"])>=2 else "medium"
        rep["analysis_level"]=level
        ranked.append(rep)

    ranked.sort(key=lambda x:x.get("final_score",0.0),reverse=True)
    best = ranked[:min(return_n,5)]
    best[0] = _ligand_guided_refine(protein_pdb, ligand_pdbqt, best[0])

    drift_events=0
    for i,pocket in enumerate(best):
        clamped=clamp_box_to_protein(protein_pdb, pocket, ligand_pdbqt=ligand_pdbqt)
        if clamped.get("_clamped",False): drift_events+=1
        best[i]=clamped

    logger.info("Done: %d pockets, %d clusters, drift_events=%d, primary_score=%.3f",
                len(best),len(clusters),drift_events,best[0].get("final_score",0))

    return {"primary":best[0],"modes":best,"total_found":len(all_modes),
            "methods_used":methods_used,"analysis_level":level,
            "performance_mode":"comprehensive" if detailed else "optimized",
            "consensus_clusters":len(clusters),"drift_events":drift_events}

def detect_binding_pockets(*args, **kwargs): return identify_binding_site(*args, **kwargs)
def get_pocket_analysis_config(pipeline_mode="production"):
    configs = {"development":{"detailed":True,"use_validated":True,"return_n":5},
               "production":{"detailed":True,"use_validated":True,"return_n":3},
               "fast":{"detailed":False,"use_validated":False,"return_n":3},
               "comprehensive":{"detailed":True,"use_validated":True,"return_n":5}}
    return configs.get(pipeline_mode, configs["production"])