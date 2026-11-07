
# Shared path settings; scripts can be invoked from any working directory.
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "lib"))
from lsp_paths import WORKFLOW, DATA, ANALYSIS, RF_ROOT, env_path
_Path(ANALYSIS).mkdir(parents=True, exist_ok=True)

import os
import re
import math
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

BEST_PDB_DIR = env_path('LSP_INPUT_DIR', f'{DATA}/xiao_haddock3_cyclic_batch/best_pdb_all')
SCORES_TSV = env_path('LSP_SCORES_TSV', f'{DATA}/xiao_haddock3_cyclic_batch/global_score_with_prodigy.tsv')
OUTPUT_DIR = ANALYSIS

POCKET1_RES = {118, 136}
POCKET2_RES = {129, 130, 132, 140}
CUTOFF = 4.0


def parse_pdb_atoms(pdb_path, target_residues):
    a_atoms_by_res = {r: [] for r in target_residues}
    b_atoms = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            if len(line) < 54:
                continue
            chain = line[21].strip()
            if chain not in ('A', 'B'):
                continue
            try:
                resi = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except Exception:
                continue
            element = line[76:78].strip().upper() if len(line) >= 78 else ''
            if element == 'H':
                continue
            if chain == 'A' and resi in a_atoms_by_res:
                a_atoms_by_res[resi].append((x, y, z))
            elif chain == 'B':
                b_atoms.append((x, y, z))
    return a_atoms_by_res, b_atoms


def min_distance(a_atoms_by_res, b_atoms, residues):
    m = float('inf')
    for r in residues:
        atoms_a = a_atoms_by_res.get(r, [])
        if not atoms_a:
            continue
        for ax, ay, az in atoms_a:
            for bx, by, bz in b_atoms:
                dx = ax - bx
                dy = ay - by
                dz = az - bz
                d = math.sqrt(dx * dx + dy * dy + dz * dz)
                if d < m:
                    m = d
    return m


def load_delta_g_map(tsv_path):
    df = pd.read_csv(tsv_path, sep='\t')
    model_col = next((c for c in df.columns if 'Model' in c), None)
    dg_col = next((c for c in df.columns if ('ΔG' in c) or ('dG' in c) or ('G(' in c)), None)
    if not model_col or not dg_col:
        raise ValueError('Could not locate Model or ΔG column in TSV')
    return dict(zip(df[model_col].astype(str), pd.to_numeric(df[dg_col], errors='coerce')))


def derive_model_id_from_filename(filename):
    return os.path.splitext(filename)[0].removesuffix('_emref_1')


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    dg_map = load_delta_g_map(SCORES_TSV)

    target_residues = POCKET1_RES.union(POCKET2_RES)
    results = []

    for fname in sorted(os.listdir(BEST_PDB_DIR)):
        if not fname.endswith('.pdb'):
            continue
        fpath = os.path.join(BEST_PDB_DIR, fname)
        a_atoms_by_res, b_atoms = parse_pdb_atoms(fpath, target_residues)
        if not b_atoms:
            # If no B chain atoms found, classify as other with infinite distances
            min1 = float('inf')
            min2 = float('inf')
        else:
            min1 = min_distance(a_atoms_by_res, b_atoms, POCKET1_RES)
            min2 = min_distance(a_atoms_by_res, b_atoms, POCKET2_RES)

        if min1 < CUTOFF:
            pocket = 'pocket1'
        elif min2 < CUTOFF:
            pocket = 'pocket2'
        else:
            pocket = 'other'

        model_id = derive_model_id_from_filename(fname)
        delta_g = dg_map.get(model_id, float('nan'))

        results.append({
            'pdb_file': fname,
            'model_id': model_id,
            'pocket': pocket,
            'min_dist_pocket1': min1,
            'min_dist_pocket2': min2,
            'ΔG(kcal/mol)': delta_g,
        })

    df_res = pd.DataFrame(results)
    all_tsv = os.path.join(OUTPUT_DIR, 'pocket_classification.tsv')
    df_res.to_csv(all_tsv, sep='\t', index=False)

    avg_df = df_res.groupby('pocket', dropna=False)['ΔG(kcal/mol)'].mean().reset_index()
    avg_tsv = os.path.join(OUTPUT_DIR, 'pocket_deltaG_avg.tsv')
    avg_df.to_csv(avg_tsv, sep='\t', index=False)

    # Set scientific plotting style
    sns.set_theme(style='whitegrid', context='paper')
    plt.rcParams.update({
        'font.family': 'DejaVu Sans',
        'font.size': 10,
        'axes.labelsize': 12,
        'axes.titlesize': 12,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'pdf.fonttype': 42,  # Editable text in PDF
    })

    # Bar chart with fixed order and tighter horizontal spacing
    order = ['pocket1', 'pocket2', 'other']
    display_labels = ['Pocket 1', 'Pocket 2', 'Other']
    count_df = df_res.groupby('pocket', dropna=False).size().reset_index(name='count')
    counts_map = dict(zip(count_df['pocket'], count_df['count']))
    categories = order
    counts = [counts_map.get(cat, 0) for cat in categories]
    colors = sns.color_palette('Set2', n_colors=len(categories))

    fig, ax = plt.subplots(figsize=(4.8, 4))
    bar_width = 0.55  # balanced width; reduce gap but keep bars slender
    ax.bar(display_labels, counts, width=bar_width, color=colors, edgecolor='black', linewidth=1)
    ax.set_xlabel('Pocket')
    ax.set_ylabel('Count')
    ax.set_title('Counts per pocket')
    sns.despine()
    # Only major gridlines; remove minor gridlines
    ax.grid(False)
    ax.yaxis.grid(True, which='major', linestyle='--', alpha=0.5)
    ax.minorticks_off()
    # Remove extra horizontal margins to tighten spacing
    ax.margins(x=0.02)
    plt.tight_layout()
    bar_pdf = os.path.join(OUTPUT_DIR, 'pocket_counts_bar.pdf')
    plt.savefig(bar_pdf)
    plt.close()

    plt.figure(figsize=(4.8, 4))
    ax = sns.violinplot(
        data=df_res.dropna(subset=['ΔG(kcal/mol)']),
        x='pocket', y='ΔG(kcal/mol)', order=order, palette='Set2', cut=0, scale='width', inner='quartile', linewidth=1
    )
    plt.xlabel('Pocket')
    plt.ylabel('ΔG (kcal/mol)')
    plt.title('ΔG distribution per pocket')
    # Replace category tick labels with display names
    plt.xticks(ticks=range(len(display_labels)), labels=display_labels)
    sns.despine()
    # Only major gridlines; remove minor gridlines
    ax.grid(False)
    ax.yaxis.grid(True, which='major', linestyle='--', alpha=0.4)
    ax.minorticks_off()
    plt.tight_layout()
    violin_pdf = os.path.join(OUTPUT_DIR, 'pocket_deltaG_violin.pdf')
    plt.savefig(violin_pdf)
    plt.close()

    print('Wrote files:')
    print(all_tsv)
    print(avg_tsv)
    print(bar_pdf)
    print(violin_pdf)


if __name__ == '__main__':
    main()
