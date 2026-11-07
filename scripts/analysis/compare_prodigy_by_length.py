#!/usr/bin/env python3
"""
Raincloud-style comparison of initial vs optimized ΔG grouped by B-chain AA count.
Data sources:
  - Initial ΔG: xiao_haddock3_cyclic_batch/global_score_with_prodigy.tsv
  - Optimized ΔG and AA count: xiao_cyclic_2_MPNN/final_summary_with_AA.tsv

Filter: include only samples with |ΔG_initial| ≥ 5.

Outputs:
  deltaG_by_aa_initial_vs_optimized_long.tsv
  deltaG_by_aa_raincloud_initial_vs_optimized.pdf
"""


# Shared path settings; scripts can be invoked from any working directory.
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "lib"))
from lsp_paths import WORKFLOW, DATA, ANALYSIS, RF_ROOT, env_path
_Path(ANALYSIS).mkdir(parents=True, exist_ok=True)

from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import to_rgba


BASE_DIR = Path(ANALYSIS)
INITIAL_TSV = Path(env_path("LSP_INITIAL_TSV", f"{DATA}/xiao_haddock3_cyclic_batch/global_score_with_prodigy.tsv"))
OPTIM_TSV = Path(env_path("LSP_OPTIMIZED_TSV", f"{DATA}/xiao_cyclic_2_MPNN/final_summary_with_AA.tsv"))

OUT_LONG_TSV = BASE_DIR / "deltaG_by_aa_initial_vs_optimized_long.tsv"
OUT_PDF = BASE_DIR / "deltaG_by_aa_raincloud_initial_vs_optimized.pdf"


def load_initial() -> pd.DataFrame:
    df = pd.read_csv(INITIAL_TSV, sep="\t")
    if "Model" not in df.columns or "ΔG(kcal/mol)" not in df.columns:
        raise ValueError("Initial TSV must have 'Model' and 'ΔG(kcal/mol)'")
    df = df[["Model", "ΔG(kcal/mol)"]].copy()
    df.rename(columns={"ΔG(kcal/mol)": "deltaG_initial"}, inplace=True)
    return df


def load_optimized() -> pd.DataFrame:
    df = pd.read_csv(OPTIM_TSV, sep="\t")
    required = {"Backbone", "ΔG(kcal/mol)", "B_chain_AA_count"}
    if not required.issubset(df.columns):
        raise ValueError("Optimized TSV must have columns: Backbone, ΔG(kcal/mol), B_chain_AA_count")
    df = df[["Backbone", "ΔG(kcal/mol)", "B_chain_AA_count"]].copy()
    df.rename(columns={"ΔG(kcal/mol)": "deltaG_optimized", "B_chain_AA_count": "aa_count"}, inplace=True)
    return df


def backbone_to_initial_model(backbone: str) -> str:
    idx = backbone.find("_emref")
    if idx != -1:
        return backbone[:idx]
    return backbone


def build_long() -> pd.DataFrame:
    df_init = load_initial()
    df_opt = load_optimized()

    df_opt["Initial_Model"] = df_opt["Backbone"].astype(str).map(backbone_to_initial_model)
    merged = pd.merge(df_opt, df_init, left_on="Initial_Model", right_on="Model", how="left")
    missing = merged[merged["deltaG_initial"].isna()]["Backbone"].tolist()
    if missing:
        print(f"[warn] Missing initial ΔG for {len(missing)} backbone(s). Examples: {missing[:3]}")
    merged = merged.dropna(subset=["deltaG_initial"]).copy()

    merged["initial_abs"] = merged["deltaG_initial"].abs()
    subset = merged[merged["initial_abs"] >= 5].copy()

    # Long format: two rows per sample (Initial, Optimized)
    df_initial = subset[["aa_count", "Backbone", "deltaG_initial"]].copy()
    df_initial.rename(columns={"deltaG_initial": "deltaG"}, inplace=True)
    df_initial["state"] = "Initial"

    df_optimized = subset[["aa_count", "Backbone", "deltaG_optimized"]].copy()
    df_optimized.rename(columns={"deltaG_optimized": "deltaG"}, inplace=True)
    df_optimized["state"] = "Optimized"

    long_df = pd.concat([df_initial, df_optimized], axis=0, ignore_index=True)
    long_df["aa_count"] = long_df["aa_count"].astype(int)
    long_df.to_csv(OUT_LONG_TSV, sep="\t", index=False)
    return long_df


def plot_raincloud(long_df: pd.DataFrame):
    # Half violins with jittered points: initial on the left, optimized on the right.
    sns.set_theme(style="whitegrid", context="paper")
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["font.family"] = "DejaVu Sans"

    order = sorted(long_df["aa_count"].unique())
    palette = {"Initial": "#1f77b4", "Optimized": "#ff7f0e"}

    fig, ax = plt.subplots(figsize=(max(8, 1.25 * len(order)), 5.8), constrained_layout=True)

    # Draw centered violins separately without dodge, then clip each category into halves.
    n_before = len(ax.collections)
    sns.violinplot(
        data=long_df[long_df["state"] == "Initial"],
        x="aa_count", y="deltaG",
        order=order,
        color=palette["Initial"],
        inner=None,
        cut=0,
        density_norm="width",
        linewidth=0.8,
        ax=ax,
    )
    initial_patches = ax.collections[n_before:]

    n_before_opt = len(ax.collections)
    sns.violinplot(
        data=long_df[long_df["state"] == "Optimized"],
        x="aa_count", y="deltaG",
        order=order,
        color=palette["Optimized"],
        inner=None,
        cut=0,
        density_norm="width",
        linewidth=0.8,
        ax=ax,
    )
    optimized_patches = ax.collections[n_before_opt:]

    # Build clipping rectangles in data coordinates.
    y_min, y_max = long_df["deltaG"].min(), long_df["deltaG"].max()
    y_pad = 0.05 * (y_max - y_min + 1e-9)
    total_height = (y_max - y_min) + 2 * y_pad

    # Use transparent fills and colored outlines: initial left, optimized right.
    initial_rgba_transparent = to_rgba(palette["Initial"], alpha=0.0)
    optimized_rgba_transparent = to_rgba(palette["Optimized"], alpha=0.0)

    # Clip each category at x=i: initial left, optimized right.
    for i, _aa in enumerate(order):
        # Left rectangle: [i-0.5, i].
        left_rect = plt.Rectangle((i - 0.5, y_min - y_pad), 0.5, total_height, transform=ax.transData)
        # Right rectangle: [i, i+0.5].
        right_rect = plt.Rectangle((i, y_min - y_pad), 0.5, total_height, transform=ax.transData)

        # Each violinplot call creates one PolyCollection per category.
        if i < len(initial_patches):
            p = initial_patches[i]
            p.set_clip_path(left_rect)
            p.set_facecolor(initial_rgba_transparent)
            p.set_edgecolor(palette["Initial"])
            p.set_linewidth(1.0)
            p.set_zorder(3)
        if i < len(optimized_patches):
            p = optimized_patches[i]
            p.set_clip_path(right_rect)
            p.set_facecolor(optimized_rgba_transparent)
            p.set_edgecolor(palette["Optimized"])
            p.set_linewidth(1.0)
            p.set_zorder(3)

    # Overlay jittered points within each half: initial left, optimized right.
    rng = np.random.default_rng(42)
    jitter_max = 0.18
    size_pts = 16

    for i, aa in enumerate(order):
        ys_initial = long_df[(long_df["state"] == "Initial") & (long_df["aa_count"] == aa)]["deltaG"].to_numpy()
        ys_optimized = long_df[(long_df["state"] == "Optimized") & (long_df["aa_count"] == aa)]["deltaG"].to_numpy()

        xs_initial = i + rng.uniform(-jitter_max, 0.0, size=ys_initial.size)
        xs_optimized = i + rng.uniform(0.0, jitter_max, size=ys_optimized.size)

        ax.scatter(xs_initial, ys_initial, s=size_pts, alpha=0.65, color=palette["Initial"], edgecolor="none", zorder=10, label=None)
        ax.scatter(xs_optimized, ys_optimized, s=size_pts, alpha=0.65, color=palette["Optimized"], edgecolor="none", zorder=10, label=None)

    # Create a custom legend for the two states.
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=palette['Initial'], markersize=6, label='Initial'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=palette['Optimized'], markersize=6, label='Optimized'),
    ]
    ax.legend(handles=handles, title="State", frameon=False, loc="upper right")

    ax.set_xlabel("B-chain AA Count")
    ax.set_ylabel("ΔG (kcal/mol)")
    ax.set_title("Initial vs Optimized ΔG by AA Count (Half Violins + Jitter, |ΔG_initial| ≥ 5)")
    ax.grid(True, linestyle="--", linewidth=0.6, alpha=0.35)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    fig.savefig(OUT_PDF, dpi=300)
    plt.close(fig)


def main():
    long_df = build_long()
    print(f"[info] AA counts plotted: {sorted(long_df['aa_count'].unique())}")
    plot_raincloud(long_df)
    print(f"[done] Saved: {OUT_LONG_TSV.name}, {OUT_PDF.name}")


if __name__ == "__main__":
    main()
