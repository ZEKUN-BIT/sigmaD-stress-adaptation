#!/usr/bin/env python3
"""Plot three stacked growth-model heatmaps; all rates use min^-1.

Run: python plot_growth_heatmap.py --input source_data.csv
Dependencies: numpy, pandas, matplotlib
"""

from pathlib import Path
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from calculate_growth_model import CONDITIONS as MODEL_CONDITIONS
from calculate_growth_model import STRAINS, calculate, load_rates

CONDITIONS = MODEL_CONDITIONS[1:]
TITLES = CONDITIONS
SIZE_MM = (183, 210)
DISPLAY_SCALE = 1e-3  # A displayed value of 1 means 0.001 min^-1.


def prepare_data(source):
    results, curves = calculate(load_rates(source))
    p = np.sort(curves["p_stress"].unique())
    matrices = {}
    for condition in CONDITIONS:
        selected = curves[curves["stress"] == condition]
        matrices[condition] = selected.pivot(
            index="strain", columns="p_stress", values="s_min_inv"
        ).loc[STRAINS, p].to_numpy()
    summary = results[[
        "strain", "stress", "delta_normal_min_inv", "delta_stress_min_inv",
        "classification", "p_star", "advantage",
    ]].copy()
    return p, matrices, summary, curves


def create_figure(p, matrices, summary):
    arial_path = font_manager.findfont("Arial", fallback_to_default=False)
    plt.rcParams.update({
        "font.family": "Arial",
        "font.sans-serif": ["Arial"],
        "font.size": 9.5, "axes.labelsize": 8, "axes.titlesize": 12,
        "xtick.labelsize": 9.5, "ytick.labelsize": 9.5,
        "mathtext.fontset": "custom", "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
        "mathtext.sf": "Arial", "mathtext.default": "regular",
        "mathtext.fallback": None,
        "axes.linewidth": 0.5, "xtick.major.width": 0.5,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "axes.unicode_minus": True, "text.usetex": False,
        "figure.constrained_layout.use": False,
    })
    fig = plt.figure(figsize=np.array(SIZE_MM) / 25.4, dpi=150, facecolor="white")
    grid = fig.add_gridspec(3, 1, left=0.105, right=0.98, bottom=0.33, top=0.92, hspace=0.34)
    max_abs = max(float(np.max(np.abs(matrix))) for matrix in matrices.values())
    limit = float(np.ceil(max_abs / DISPLAY_SCALE))
    norm = TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)
    cmap = plt.get_cmap("RdBu_r")
    p_edges = np.r_[0, (p[1:] + p[:-1])/2, 1] * 100
    strain_edges = np.arange(len(STRAINS) + 1) - 0.5
    metadata = {"size_mm": list(SIZE_MM), "rate_unit": "min^-1", "display_scale": DISPLAY_SCALE,
                "color_limits_display_units": [-limit, limit], "p_points": len(p),
                "layout": "stacked_transposed", "column_order": STRAINS,
                "panel_order": CONDITIONS, "stress_fraction_increases": "upward",
                "strain_label_rotation_degrees": 45,
                "threshold_markers": int(summary.p_star.notna().sum()),
                "all_p_symbols": int(summary.p_star.isna().sum()), "row_standardization": False,
                "clipping": False, "interpolation_between_strains": False,
                "font_family": "Arial", "font_path": arial_path,
                "tick_font_pt": 9.5, "axis_label_font_pt": 8,
                "panel_title_font_pt": 12, "overall_title_font_pt": 12,
                "colour_map": "RdBu_r", "panel_letters": False}

    for k, (condition, title) in enumerate(zip(CONDITIONS, TITLES)):
        inner = grid[k].subgridspec(2, 1, height_ratios=[0.10, 1], hspace=0.13)
        ax = fig.add_subplot(inner[1])
        status_ax = fig.add_subplot(inner[0], sharex=ax)
        ax.set_gid(f"heatmap-{condition}")
        status_ax.set_gid(f"all-p-{condition}")
        mesh = ax.pcolormesh(strain_edges, p_edges, matrices[condition].T/DISPLAY_SCALE,
                             shading="flat", cmap=cmap, norm=norm, rasterized=False)
        ax.vlines(strain_edges[1:-1], 0, 100, colors="white", linewidth=0.22, alpha=0.65)
        ax.set(xlim=(-0.5, len(STRAINS)-0.5), ylim=(0, 100), yticks=[0, 25, 50, 75, 100])
        ax.set_xticks(np.arange(len(STRAINS)))
        ax.set_xticklabels(STRAINS, rotation=45, ha="right", va="top", rotation_mode="anchor")
        ax.tick_params(axis="x", labelbottom=k == len(CONDITIONS)-1,
                       bottom=k == len(CONDITIONS)-1, pad=3)
        ax.tick_params(axis="y", pad=4)
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_edgecolor("#777777")
        status_ax.set(ylim=(0, 1))
        status_ax.axis("off")
        status_ax.set_title(title, pad=6, fontweight="bold")
        status_ax.text(-0.01, 0.5, "All p", ha="right", va="center", fontsize=6.5,
                       transform=status_ax.transAxes)
        selected = summary[summary.stress == condition].set_index("strain").loc[STRAINS]
        for i, row in enumerate(selected.itertuples()):
            if pd.notna(row.p_star):
                marker = "^" if row.advantage == "p > p_star" else "v"
                ax.plot(i, row.p_star * 100, marker=marker, color="#111111", markersize=3.2,
                        markeredgewidth=0.3, linestyle="none", clip_on=False, zorder=4)
            else:
                sign = "+" if row.classification == "beneficial_both" else "\u2212"
                status_ax.text(i, 0.5, sign, ha="center", va="center", fontsize=8, fontweight="bold")

    fig.text(0.5, 0.985, "Model-predicted growth advantage", ha="center", va="top", fontsize=12)
    fig.text(0.025, 0.625, "Stress time fraction, p (%)", ha="center", va="center",
             rotation=90, fontsize=8)
    fig.text(0.5425, 0.248, "Mutant", ha="center", va="center", fontsize=8)
    color_ax = fig.add_axes([0.235, 0.18, 0.55, 0.016])
    colorbar = fig.colorbar(mesh, cax=color_ax, orientation="horizontal", ticks=[-limit, -limit/2, 0, limit/2, limit])
    colorbar.solids.set_rasterized(False)
    colorbar.set_label(r"Relative log-growth rate, $s(p)$ ($10^{-3}$ min$^{-1}$)", labelpad=3, fontsize=8)
    colorbar.outline.set_linewidth(0.5)
    colorbar.ax.tick_params(labelsize=9.5, length=2)
    fig.text(0.235, 0.209, "WT advantage", ha="left", fontsize=7, color="#2166AC")
    fig.text(0.785, 0.209, "Mutant advantage", ha="right", fontsize=7, color="#B2182B")
    legend_handles = [
        Line2D([], [], marker="^", color="black", linestyle="none", markersize=4, label=r"$p^*$; advantage at higher $p$"),
        Line2D([], [], marker="v", color="black", linestyle="none", markersize=4, label=r"$p^*$; advantage at lower $p$"),
        Line2D([], [], marker="$+$", color="black", linestyle="none", markersize=6, label="All p: mutant advantage"),
        Line2D([], [], marker="$\u2212$", color="black", linestyle="none", markersize=6, label="All p: WT advantage"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", bbox_to_anchor=(0.52, 0.055),
               ncol=2, frameon=False, fontsize=7, columnspacing=2.5, handletextpad=0.6,
               labelspacing=0.7, handlelength=1.3)
    fig.text(0.52, 0.015, "Point estimates from supplied \u03bcmax values; uncertainty not estimated.",
             ha="center", va="center", fontsize=7, color="#444444")
    return fig, metadata


def export_figure(fig, output):
    output.mkdir(parents=True, exist_ok=True)
    # Keep the specified physical canvas size; do not crop with bbox_inches='tight'.
    for suffix in ["pdf", "svg", "png"]:
        fig.savefig(output / f"growth_advantage_heatmap_vertical.{suffix}", dpi=600 if suffix == "png" else 300, facecolor="white")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Path to a user-supplied input CSV.")
    parser.add_argument("--output", type=Path, default=Path("figures"))
    args = parser.parse_args()
    if not args.input.is_file():
        parser.error(f"Input CSV not found: {args.input}")
    p, matrices, summary, curves = prepare_data(args.input)
    fig, _ = create_figure(p, matrices, summary)
    export_figure(fig, args.output)
    summary.to_csv(args.output / "threshold_summary.csv", index=False, encoding="utf-8-sig")
    curves.to_csv(args.output / "heatmap_values.csv", index=False, encoding="utf-8-sig")
    plt.close(fig)
    print(f"Saved 42 \u00d7 101 \u00d7 3 heatmap and source values to {args.output}")


if __name__ == "__main__":
    main()
