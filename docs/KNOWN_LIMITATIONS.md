# Known Limitations of the Dippy Dock Pipeline

## 1. Search-Volume Ceiling for Large Ligands

**Status:** Partially mitigated (ligand-aware box sizing now implemented).

The per-axis box maximum was originally hard-capped at 26 Å regardless of
ligand size.  We have replaced this with a tiered ligand-aware maximum:

| Heavy atoms | Per-axis max |
|-------------|-------------|
| ≤ 5         | 26 Å (legacy) |
| 6–30        | 36 Å |
| > 30        | 44 Å |

**Remaining caveat:** Even at 44 Å the search volume is finite.  Very large
peptides (> 50 heavy atoms) or ligands that span multiple sub-pockets may
still require a bigger box than the pipeline allows.  Dedicated ensemble
docking tools (e.g., DINC) exist for this regime and are outside the scope
of this pipeline.

**Paper wording suggestion:** "We impose a per-axis upper bound on the
docking search volume that scales with ligand heavy-atom count (26–44 Å).
This ceiling is deliberately conservative to keep Vina grid dimensions
tractable, but it means that very large or extended ligands — such as
macrocyclic peptides — may require a bespoke protocol."

---

## 2. Exhaustiveness Scaling

**Status:** Now adaptive.

Exhaustiveness is no longer fixed at 8; it is computed as:

    exh = base × (V_box / 8000)^{1/3} × (1 + 0.15 × max(0, n_rot − 6))

with a hard cap of 32.  This keeps exhaustiveness at 8 for the reference
case (20³ Å box, 6 rotatable bonds) and increases it for larger boxes or
more flexible ligands.

**Remaining caveat:** The relationship between exhaustiveness and search
quality is not linear, and the Vina developers have not published a
general scaling law.  Our formula is a reasonable heuristic, not a
theoretically grounded scaling.  For publication-critical results on large
ligands, we recommend manual tuning with the Vina `--cpu` flag and
multiple exhaustiveness values.

---

## 3. Rotatable-Bond Accuracy Wall (Vina Limitation)

**Status:** Acknowledged — not fixable within Vina.

Docking ligands with ≤ 6 rotatable bonds is generally fast and accurate.
As the dimensionality of the ligand's conformational space increases, Vina's
Lamarckian genetic algorithm struggles to sample adequately, even when
exhaustiveness is increased.  This is a known, hard accuracy wall for the
AutoDock family of docking programs.

**Production caps observed in hosted Vina services:**
- Max 32 rotatable bonds
- Max 150 heavy atoms
- Max 300 total atoms

**Paper wording suggestion:** "AutoDock Vina's search quality degrades as
conformational complexity increases.  We adopt a tiered exhaustiveness
strategy to partially compensate, but acknowledge that dedicated
large-ligand protocols (e.g., DINC) may be required for ligands exceeding
~30 heavy atoms or ~15 rotatable bonds."

---

## 4. RMSD Validation: Atom-Order Sensitivity

**Status:** Mitigated (pose extraction + Kabsch alignment now implemented).

The previous implementation had two bugs:

1. **Multi-model merge bug:** `Chem.MolFromPDBFile` was called on the
   entire multi-model PDBQT output, causing RDKit to merge all poses into
   one molecule.  The `pose_index` parameter was ignored.  **Fixed:** We
   now extract only MODEL `pose_index + 1` into a temporary PDB file before
   RMSD computation.

2. **Naive atom-order fallback:** The coordinate RMSD fallback did a
   1:1 `zip()` of crystal and docked atom lists without reordering or
   alignment.  For large/flexible ligands where atom ordering diverges,
   this could return a silently incorrect (too-low) RMSD.  **Fixed:** We
   now use Kabsch alignment (optimal rotation) before computing RMSD.

**Remaining caveat:** RDKit's `GetBestRMS` uses substructure matching for
atom reordering, which works well for drug-like molecules but can fail for
highly modified or unusual scaffolds.  In such cases the Kabsch fallback
is used, which does not reorder atoms — it only finds the optimal rotation.
If atom ordering is significantly different, even the Kabsch RMSD may be
misleading.

---

## 5. Adaptive Docking Timeout

**Status:** Implemented.

The previous fixed 600-second timeout was replaced with an adaptive
computation based on box volume and ligand flexibility:

    timeout = 600 × (V_box / 8000)^{1/2} × (1 + 0.08 × max(0, n_rot − 6))

clamped to [300, 1800] seconds.  This gives:
- Small boxes + rigid ligands: ~5 min
- Medium boxes + average ligands: ~10 min
- Large boxes + flexible ligands: up to 30 min

**Remaining caveat:** The Vina subprocess itself has no internal timeout.
The Python-level timeout uses `concurrent.futures.ThreadPoolExecutor`,
which raises `TimeoutError` in the calling thread but may leave the Vina
process running.  A future improvement would be to use `subprocess.run(timeout=...)`
or `psutil` to ensure the process is killed.

---

## 6. Failure Categorisation

**Status:** Implemented.

The `failure_category` field in per-target `result.json` files now classifies
failures into:

| Category | Description |
|----------|-------------|
| `ligand_prep` | RDKit/Meeko/OpenBabel embedding or conversion failed |
| `protein_prep` | Protein preprocessor (cleaning/PDBQT) failed |
| `docking_timeout` | Vina subprocess timed out or was killed |
| `docking_failed` | Vina ran but returned no valid poses |
| `pocket_not_found` | No binding site could be identified |
| `rmsd_error` | Docking succeeded but RMSD computation failed |
| `download_error` | PDB / UniProt fetch failed |
| `unknown` | Could not determine category |

These categories are aggregated in `statistics.json` under the
`failure_categories` key, enabling the per-failure-type breakdown
recommended by the reviewer.

---

## 7. Benchmark Dataset Limitations

The Astex Diverse Set contains 85 crystallographic complexes with mostly
drug-like ligands (median ~25 heavy atoms, median ~5 rotatable bonds).
Performance on this set does not generalise to:

- **Peptide ligands** (> 50 heavy atoms, > 20 rotatable bonds)
- **Macrocycles** (conformational rigidity helps Vina, but large size hurts)
- **Covalent inhibitors** (Vina does not model covalent bond formation)
- **Metalloenzymes** (metal coordination is poorly handled by Vina's scoring)

The `failure_categories` + ligand-property table in `statistics.json`
provides the evidence base for making these scope claims precisely.
