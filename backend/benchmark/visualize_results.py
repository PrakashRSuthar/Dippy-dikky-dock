# backend/benchmark/visualize_results.py
# Generate publication-quality figures from Astex benchmark results

import os, sys, json, math
from pathlib import Path
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 11,
    "axes.labelsize": 13,
    "axes.titlesize": 14,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.linewidth": 1.2,
    "xtick.major.width": 1.0,
    "ytick.major.width": 1.0,
})

RCSB_COLORS = {
    "very_good": "#2d6a4f",
    "good": "#52b788",
    "acceptable": "#ffd166",
    "poor": "#e63946",
}

def load_results(results_dir):
    results = []
    for entry_dir in sorted(Path(results_dir).iterdir()):
        if not entry_dir.is_dir():
            continue
        result_file = entry_dir / "result.json"
        if result_file.exists():
            try:
                results.append(json.loads(result_file.read_text()))
            except:
                pass
    return results

def classify_rmsd(rmsd):
    if rmsd is None:
        return "failed"
    if rmsd <= 1.0:
        return "very_good"
    elif rmsd <= 2.0:
        return "good"
    elif rmsd <= 3.0:
        return "acceptable"
    else:
        return "poor"

def compute_statistics(results):
    completed = [r for r in results if r["status"] == "completed"]
    failed = [r for r in results if r["status"] == "failed"]
    rmsd_vals = [r["rmsd"] for r in completed if r["rmsd"] is not None]
    aff_vals = [r["best_affinity"] for r in completed if r["best_affinity"] is not None]

    stats = {
        "total": len(results),
        "completed": len(completed),
        "failed": len(failed),
        "success_2a": sum(1 for r in completed if r.get("success")),
        "success_1a": sum(1 for r in completed if r["rmsd"] is not None and r["rmsd"] <= 1.0),
        "n_rmsd": len(rmsd_vals),
        "n_aff": len(aff_vals),
    }
    if rmsd_vals:
        stats.update({
            "mean_rmsd": float(np.mean(rmsd_vals)),
            "median_rmsd": float(np.median(rmsd_vals)),
            "std_rmsd": float(np.std(rmsd_vals)),
            "min_rmsd": float(np.min(rmsd_vals)),
            "max_rmsd": float(np.max(rmsd_vals)),
            "q1_rmsd": float(np.percentile(rmsd_vals, 25)),
            "q3_rmsd": float(np.percentile(rmsd_vals, 75)),
        })
    if aff_vals:
        stats.update({
            "mean_aff": float(np.mean(aff_vals)),
            "min_aff": float(np.min(aff_vals)),
            "max_aff": float(np.max(aff_vals)),
        })
    # Classification breakdown (completed targets)
    classes = Counter()
    for r in completed:
        cls = classify_rmsd(r.get("rmsd"))
        classes[cls] += 1
    stats["classes"] = dict(classes)

    # Failure breakdown by category
    failure_cats = Counter()
    for r in failed:
        cat = r.get("failure_category") or "unknown"
        failure_cats[cat] += 1
    stats["failure_categories"] = dict(failure_cats)

    # Ligand property statistics (for completed + failed)
    all_rot = [r.get("ligand_rotatable_bonds") for r in results if r.get("ligand_rotatable_bonds") is not None]
    all_ha  = [r.get("ligand_heavy_atoms") for r in results if r.get("ligand_heavy_atoms") is not None]
    if all_rot:
        stats["ligand_rotatable_bonds"] = {
            "mean": float(np.mean(all_rot)),
            "median": float(np.median(all_rot)),
            "max": int(np.max(all_rot)),
        }
    if all_ha:
        stats["ligand_heavy_atoms"] = {
            "mean": float(np.mean(all_ha)),
            "median": float(np.median(all_ha)),
            "max": int(np.max(all_ha)),
        }

    return stats, rmsd_vals, aff_vals, completed

def fig_rmsd_histogram(rmsd_vals, stats, output_dir):
    fig, ax = plt.subplots(figsize=(7, 5))
    bins = np.linspace(0, max(rmsd_vals + [10.0]), 40)
    colors = ["#2d6a4f" if b <= 2.0 else "#ffd166" if b <= 3.0 else "#e63946" for b in bins[:-1]]
    n, bins_p, patches = ax.hist(rmsd_vals, bins=bins, edgecolor="white", linewidth=0.5)
    for patch, color in zip(patches, colors):
        patch.set_facecolor(color)

    ax.axvline(2.0, color="#d00000", linestyle="--", linewidth=1.5, label="Success threshold (2.0 Å)")
    ax.axvline(stats["median_rmsd"], color="#023e8a", linestyle=":", linewidth=1.5,
               label=f'Median = {stats["median_rmsd"]:.2f} Å')

    ax.set_xlabel("RMSD (Å)")
    ax.set_ylabel("Number of targets")
    ax.set_title("RMSD Distribution — Astex Diverse Set")
    ax.legend(frameon=True, fancybox=False)
    ax.set_xlim(0, min(max(rmsd_vals) + 1, 12))
    fig.tight_layout()
    fig.savefig(output_dir / "rmsd_distribution.png")
    fig.savefig(output_dir / "rmsd_distribution.pdf")
    plt.close(fig)

def fig_success_rate(stats, output_dir):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

    categories = ["RMSD ≤ 1.0 Å", "RMSD ≤ 2.0 Å", "RMSD ≤ 3.0 Å", "Failed"]
    n_completed = stats["completed"]
    values = [
        stats["success_1a"],
        stats["success_2a"],
        n_completed - stats["classes"].get("poor", 0),
        stats["failed"],
    ]
    colors_bar = ["#2d6a4f", "#52b788", "#ffd166", "#e63946"]

    bars = ax1.bar(categories, values, color=colors_bar, edgecolor="black", linewidth=0.8, width=0.6)
    for bar, val in zip(bars, values):
        pct = val / stats["total"] * 100
        ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                 f"{val} ({pct:.1f}%)", ha="center", va="bottom", fontsize=10)
    ax1.set_ylabel("Number of targets")
    ax1.set_title("Success Rate at Different Thresholds")
    ax1.set_ylim(0, max(values) * 1.2)

    # Pie chart: classification breakdown
    classes = stats["classes"]
    labels = {"very_good": "Very Good ($\\leq$1.0 Å)", "good": "Good (1–2 Å)",
              "acceptable": "Acceptable (2–3 Å)", "poor": "Poor ($>$3.0 Å)", "failed": "Failed"}
    pie_labels = []
    pie_values = []
    pie_colors = []
    for key in ["very_good", "good", "acceptable", "poor", "failed"]:
        val = classes.get(key, 0)
        if val > 0:
            pie_labels.append(labels[key])
            pie_values.append(val)
            pie_colors.append(RCSB_COLORS.get(key, "#999999"))

    wedges, texts, autotexts = ax2.pie(
        pie_values, labels=pie_labels, autopct="%1.1f%%",
        colors=pie_colors, startangle=90,
        textprops={"fontsize": 9},
        wedgeprops={"edgecolor": "white", "linewidth": 0.8},
    )
    for t in autotexts:
        t.set_fontsize(9)
    ax2.set_title("RMSD Classification Breakdown")

    fig.tight_layout()
    fig.savefig(output_dir / "success_rate.png")
    fig.savefig(output_dir / "success_rate.pdf")
    plt.close(fig)

def fig_rmsd_vs_affinity(results, output_dir):
    completed = [r for r in results if r["status"] == "completed"
                 and r["rmsd"] is not None and r["best_affinity"] is not None]
    rmsds = [r["rmsd"] for r in completed]
    affs = [r["best_affinity"] for r in completed]
    labels = [f'{r["pdb"]}_{r["ligand_code"]}' for r in completed]

    fig, ax = plt.subplots(figsize=(8, 6))
    colors = [RCSB_COLORS[classify_rmsd(r)] for r in rmsds]
    sc = ax.scatter(affs, rmsds, c=colors, s=40, edgecolors="black", linewidth=0.5, alpha=0.85)

    # Correlation
    if len(rmsds) > 5:
        coeffs = np.polyfit(affs, rmsds, 1)
        line_x = np.linspace(min(affs), max(affs), 100)
        line_y = np.polyval(coeffs, line_x)
        ax.plot(line_x, line_y, color="#023e8a", linestyle="--", linewidth=1.0,
                label=f"Linear fit (r={np.corrcoef(affs, rmsds)[0,1]:.3f})")
        ax.legend(frameon=True, fancybox=False)

    ax.axhline(2.0, color="#d00000", linestyle=":", linewidth=1.2, alpha=0.7)
    ax.set_xlabel("Binding Affinity (kcal/mol)")
    ax.set_ylabel("RMSD (Å)")
    ax.set_title("RMSD vs Binding Affinity")

    # Legend for color classes
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#2d6a4f", markersize=8, label="Very Good ($\\leq$1.0 Å)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#52b788", markersize=8, label="Good (1–2 Å)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#ffd166", markersize=8, label="Acceptable (2–3 Å)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#e63946", markersize=8, label="Poor ($>$3.0 Å)"),
    ]
    ax.legend(handles=legend_elements, frameon=True, fancybox=False, fontsize=9)

    fig.tight_layout()
    fig.savefig(output_dir / "rmsd_vs_affinity.png")
    fig.savefig(output_dir / "rmsd_vs_affinity.pdf")
    plt.close(fig)

def fig_cumulative_success(rmsd_vals, output_dir):
    fig, ax = plt.subplots(figsize=(7, 5))

    sorted_rmsd = np.sort(rmsd_vals)
    n = len(sorted_rmsd)
    cumulative = np.arange(1, n + 1) / n

    ax.plot(sorted_rmsd, cumulative, color="#1b4332", linewidth=2.0)
    ax.fill_between(sorted_rmsd, cumulative, alpha=0.15, color="#52b788")

    for threshold, color in [(1.0, "#2d6a4f"), (2.0, "#d00000"), (3.0, "#e63946")]:
        idx = np.searchsorted(sorted_rmsd, threshold, side="right")
        pct = idx / n * 100
        ax.axvline(threshold, color=color, linestyle="--", linewidth=1.0, alpha=0.6)
        ax.text(threshold + 0.05, max(cumulative) * 0.15, f"{threshold} Å ({pct:.0f}%)",
                rotation=90, fontsize=9, color=color, va="bottom")

    ax.set_xlabel("RMSD Threshold (Å)")
    ax.set_ylabel("Cumulative Fraction of Targets")
    ax.set_title("Cumulative Success Rate")
    ax.set_xlim(0, min(max(sorted_rmsd) + 0.5, 12))
    ax.set_ylim(0, 1.05)
    ax.yaxis.set_major_formatter(ticker.PercentFormatter(1.0))
    ax.grid(True, alpha=0.3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(output_dir / "cumulative_success.png")
    fig.savefig(output_dir / "cumulative_success.pdf")
    plt.close(fig)

def fig_per_target_rmsd(results, output_dir, max_targets=85):
    completed = [r for r in results if r["status"] == "completed" and r["rmsd"] is not None]
    completed.sort(key=lambda r: r["rmsd"])
    completed = completed[:max_targets]
    labels = [f'{r["pdb"]}_{r["ligand_code"]}' for r in completed]
    rmsds = [r["rmsd"] for r in completed]
    colors_rmsd = [RCSB_COLORS[classify_rmsd(r)] for r in rmsds]

    fig, ax = plt.subplots(figsize=(14, max(5, len(completed) * 0.25)))
    y_pos = np.arange(len(completed))
    bars = ax.barh(y_pos, rmsds, color=colors_rmsd, edgecolor="black", linewidth=0.3, height=0.7)

    ax.axvline(2.0, color="#d00000", linestyle="--", linewidth=1.2, label="Success threshold (2.0 Å)")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("RMSD (Å)")
    ax.set_title("Per-Target RMSD (sorted)")
    ax.set_xlim(0, max(rmsds + [5.0]) * 1.1)
    ax.invert_yaxis()
    ax.legend(frameon=True, fancybox=False, fontsize=9)

    # Add text labels
    for i, (bar, rmsd) in enumerate(zip(bars, rmsds)):
        ax.text(bar.get_width() + 0.05, bar.get_y() + bar.get_height() / 2,
                f"{rmsd:.2f}", va="center", fontsize=6)

    fig.tight_layout()
    fig.savefig(output_dir / "per_target_rmsd.png")
    fig.savefig(output_dir / "per_target_rmsd.pdf")
    plt.close(fig)

def fig_stats_table(stats, output_dir):
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.axis("off")

    rows = [
        ["Metric", "Value"],
        ["Total targets", str(stats["total"])],
        ["Completed", str(stats["completed"])],
        ["Failed", str(stats["failed"])],
        ["Success rate (RMSD ≤ 2.0 Å)", f'{stats["success_2a"] / max(stats["completed"], 1) * 100:.1f}% ({stats["success_2a"]}/{stats["completed"]})'],
        ["Success rate (RMSD ≤ 1.0 Å)", f'{stats["success_1a"] / max(stats["completed"], 1) * 100:.1f}% ({stats["success_1a"]}/{stats["completed"]})'],
        ["Mean RMSD", f'{stats.get("mean_rmsd", "N/A"):.3f} Å'],
        ["Median RMSD", f'{stats.get("median_rmsd", "N/A"):.3f} Å'],
        ["Std Dev RMSD", f'{stats.get("std_rmsd", "N/A"):.3f} Å'],
        ["Q1 / Q3 RMSD", f'{stats.get("q1_rmsd", "N/A"):.3f} / {stats.get("q3_rmsd", "N/A"):.3f} Å'],
        ["Min / Max RMSD", f'{stats.get("min_rmsd", "N/A"):.3f} / {stats.get("max_rmsd", "N/A"):.3f} Å'],
        ["Mean Binding Affinity", f'{stats.get("mean_aff", "N/A"):.2f} kcal/mol'],
        ["Affinity range", f'{stats.get("min_aff", "N/A"):.2f} to {stats.get("max_aff", "N/A"):.2f} kcal/mol'],
    ]

    table = ax.table(cellText=rows, loc="center", cellLoc="left", colWidths=[0.45, 0.45])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.6)

    for (row, col), cell in table.get_celld().items():
        if row == 0:
            cell.set_text_props(fontweight="bold", fontsize=11)
            cell.set_facecolor("#1b4332")
            cell.set_text_props(color="white")
        elif row % 2 == 0:
            cell.set_facecolor("#e9ecef")
        else:
            cell.set_facecolor("white")
        cell.set_edgecolor("#dee2e6")

    ax.set_title("Astex Diverse Set — Benchmark Summary Statistics", fontweight="bold", pad=20)
    fig.tight_layout()
    fig.savefig(output_dir / "summary_statistics.png")
    fig.savefig(output_dir / "summary_statistics.pdf")
    plt.close(fig)

def generate_all(results_dir, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results = load_results(results_dir)
    if not results:
        print(f"No results found in {results_dir}")
        return

    stats, rmsd_vals, aff_vals, completed = compute_statistics(results)
    print(f"Loaded {len(results)} results ({len(completed)} completed)")

    print("Generating figures...")
    if rmsd_vals:
        fig_rmsd_histogram(rmsd_vals, stats, output_dir)
        print("  rmsd_distribution.png/pdf")
        fig_success_rate(stats, output_dir)
        print("  success_rate.png/pdf")
        fig_rmsd_vs_affinity(results, output_dir)
        print("  rmsd_vs_affinity.png/pdf")
        fig_cumulative_success(rmsd_vals, output_dir)
        print("  cumulative_success.png/pdf")
        fig_per_target_rmsd(results, output_dir)
        print("  per_target_rmsd.png/pdf")
        fig_stats_table(stats, output_dir)
        print("  summary_statistics.png/pdf")

    # Save statistics as JSON alongside figures
    stats_path = output_dir / "statistics.json"
    stats_path.write_text(json.dumps(stats, indent=2))
    print(f"\nStatistics saved to {stats_path}")
    print("All figures generated successfully.")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Generate benchmark visualization figures")
    parser.add_argument("results_dir", help="Path to benchmark results directory (containing per-target subdirs)")
    parser.add_argument("--output", "-o", default="benchmark_figures", help="Output directory for figures")
    args = parser.parse_args()
    generate_all(args.results_dir, args.output)
