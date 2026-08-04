# Generate single comprehensive PDF report of Astex benchmark results
import os, sys, json, math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.table as mtable

sys.path.insert(0, str(Path(__file__).parent))
from visualize_results import load_results, compute_statistics, classify_rmsd

RCSB_COLORS = {"very_good":"#2d6a4f","good":"#52b788","acceptable":"#ffd166","poor":"#e63946"}

OUTPUT_DIR = Path("D:\\Dippy-dikky-dock\\temp_runs\\astex_full_85\\pdf_report")

def make_title_page(stats, results):
    fig, ax = plt.subplots(figsize=(8.5, 11))
    ax.axis("off")
    lines = []
    lines.append(("Dippy-Dikky-Dock: Astex Diverse Set Benchmark", 20, True))
    lines.append(("", 12, False))
    lines.append((f"Date: {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M')}", 11, False))
    lines.append(("", 6, False))
    lines.append(("=" * 60, 10, False))
    lines.append(("", 6, False))
    lines.append(("Summary Statistics", 14, True))
    lines.append(("", 6, False))

    completed = stats["completed"]
    n_success = stats["success_2a"]
    lines.append((f"Total targets:       {stats['total']}", 11, False))
    lines.append((f"Completed:           {completed}", 11, False))
    lines.append((f"Failed:              {stats['failed']}", 11, False))
    lines.append((f"Success (RMSD<=2A):  {n_success}  ({n_success/max(completed,1)*100:.1f}%)", 11, False))
    lines.append((f"Success (RMSD<=1A):  {stats['success_1a']}  ({stats['success_1a']/max(completed,1)*100:.1f}%)", 11, False))
    if stats.get("mean_rmsd"):
        lines.append((f"Mean RMSD:           {stats['mean_rmsd']:.3f} A", 11, False))
        lines.append((f"Median RMSD:         {stats['median_rmsd']:.3f} A", 11, False))
        lines.append((f"Std Dev RMSD:        {stats['std_rmsd']:.3f} A", 11, False))
        lines.append((f"Min RMSD:            {stats['min_rmsd']:.3f} A", 11, False))
        lines.append((f"Max RMSD:            {stats['max_rmsd']:.3f} A", 11, False))
    if stats.get("mean_aff"):
        lines.append((f"Mean Binding Affinity: {stats['mean_aff']:.2f} kcal/mol", 11, False))
    lines.append(("", 6, False))
    lines.append(("=" * 60, 10, False))
    lines.append(("", 6, False))
    lines.append(("Methodology", 14, True))
    lines.append(("", 6, False))
    lines.append(("  - Docking tool: AutoDock Vina (exhaustiveness=4)", 10, False))
    lines.append(("  - Ligand preparation: RDKit ETKDG + Meeko", 10, False))
    lines.append(("  - Protein preparation: OpenBabel CLI", 10, False))
    lines.append(("  - Binding site: Crystal ligand center (self-docking)", 10, False))
    lines.append(("  - RMSD: RDKit GetBestRMS / Kabsch-aligned coordinate RMSD", 10, False))
    lines.append(("  - Success threshold: RMSD <= 2.0 A", 10, False))

    y = 0.95
    for text, size, bold in lines:
        w = 700 if bold else 400
        ax.text(0.05, y, text, fontsize=size, fontweight="bold" if bold else "normal",
                fontfamily="monospace" if not bold else "sans-serif",
                verticalalignment="top", transform=ax.transAxes)
        y -= size / 800 * 1.8

    fig.tight_layout()
    return fig

def make_per_target_table(results, stats):
    completed = [r for r in results if r["status"] == "completed"]
    completed.sort(key=lambda r: r.get("rmsd") or 999)

    # Split into pages (max 40 rows per page)
    rows_per_page = 40
    n_pages = max(1, math.ceil(len(completed) / rows_per_page))
    figs = []
    for page in range(n_pages):
        start = page * rows_per_page
        end = min(start + rows_per_page, len(completed))
        chunk = completed[start:end]
        n_rows = len(chunk) + 2
        fig, ax = plt.subplots(figsize=(8.5, max(4, n_rows * 0.35)))
        ax.axis("off")

        cell_text = [["PDB", "Lig", "RMSD(A)", "Affinity", "Duration", "Status"]]
        cell_colors = [["#1b4332"] * 6]
        for r in chunk:
            cls = classify_rmsd(r.get("rmsd"))
            c = RCSB_COLORS.get(cls, "#999999")
            aff = f'{r.get("best_affinity","?"):.1f}' if r.get("best_affinity") else "?"
            rmsd = f'{r.get("rmsd","?"):.3f}' if r.get("rmsd") else "N/A"
            dur = f'{r.get("duration_s","?"):.0f}' if r.get("duration_s") else "?"
            status = "OK" if r.get("success") else ("FAIL" if r["status"]=="failed" else "MISS")
            cell_text.append([r["pdb"], r["ligand_code"], rmsd, aff, dur, status])
            cell_colors.append([c]*6)

        table = ax.table(cellText=cell_text, cellColours=cell_colors,
                         loc="center", cellLoc="center", colWidths=[0.12]*6)
        table.auto_set_font_size(False)
        table.set_fontsize(7)
        table.scale(1, 1.3)
        for (row, col), cell in table.get_celld().items():
            if row == 0:
                cell.set_text_props(color="white", fontweight="bold", fontsize=8)
            cell.set_edgecolor("#dee2e6")
        ax.set_title(f"Per-Target Results (page {page+1}/{n_pages})", fontsize=12, fontweight="bold", pad=10)
        fig.tight_layout()
        figs.append(fig)
    return figs

def embed_existing_figures(fig_dir, pdf):
    """Embed existing generated PDF figures into the report"""
    order = ["rmsd_distribution", "cumulative_success", "success_rate",
             "rmsd_vs_affinity", "per_target_rmsd", "summary_statistics"]
    for name in order:
        path = fig_dir / f"{name}.pdf"
        if path.exists():
            try:
                im = plt.imread(str(fig_dir / f"{name}.png"))
                fig, ax = plt.subplots(figsize=(8.5, 6))
                ax.imshow(im)
                ax.axis("off")
                fig.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)
            except:
                pass

def generate_report():
    results_dir = Path("D:\\Dippy-dikky-dock\\temp_runs\\astex_full_85\\results")
    fig_dir = Path("D:\\Dippy-dikky-dock\\temp_runs\\astex_full_85\\figures")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    pdf_path = OUTPUT_DIR / "Astex_Benchmark_Report.pdf"

    results = load_results(str(results_dir))
    if not results:
        print("No results found!")
        return

    stats, rmsd_vals, _, _ = compute_statistics(results)
    print(f"Loaded {len(results)} results ({stats['completed']} completed)")

    with PdfPages(pdf_path) as pdf:
        # Title page
        print("Adding title page...")
        pdf.savefig(make_title_page(stats, results))
        plt.close()

        # Per-target tables
        print("Adding per-target tables...")
        for fig in make_per_target_table(results, stats):
            pdf.savefig(fig)
            plt.close()

        # Figures
        print("Adding figures...")
        embed_existing_figures(fig_dir, pdf)

        # Metadata
        d = pdf.infodict()
        d["Title"] = "Astex Diverse Set Benchmark - Dippy-Dikky-Dock"
        d["Author"] = "Dippy-Dikky-Dock"
        d["Subject"] = f"{stats['completed']} targets, {stats['success_2a']} successful ({stats['success_2a']/max(stats['completed'],1)*100:.1f}%)"

    print(f"\nPDF generated: {pdf_path}")
    print(f"Size: {os.path.getsize(pdf_path) / 1024 / 1024:.1f} MB")
    print(f"Pages embedded: figures from {fig_dir}")

if __name__ == "__main__":
    generate_report()
