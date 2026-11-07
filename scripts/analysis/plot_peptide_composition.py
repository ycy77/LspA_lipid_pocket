
# Shared path settings; scripts can be invoked from any working directory.
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "lib"))
from lsp_paths import WORKFLOW, DATA, ANALYSIS, RF_ROOT, env_path
_Path(ANALYSIS).mkdir(parents=True, exist_ok=True)

import os
import glob
import pickle
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import numpy as np
from collections import Counter, defaultdict
from Bio.PDB import PDBParser

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.unicode_minus'] = False

AMINO_ACIDS_SORTED = [
    'A', 'V', 'L', 'I', 'M', 'F', 'W', 'P',
    'G', 'S', 'T', 'C', 'Y', 'N', 'Q',
    'D', 'E',
    'K', 'R', 'H'
]

# 氨基酸三字母到单字母的映射
aa_3to1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

def main():
    pdb_dir = env_path('LSP_INPUT_DIR', f'{DATA}/cyclic_youhua_1')
    output_dir = f'{ANALYSIS}'
    pickle_file = os.path.join(output_dir, 'b_chain_peptide_data.pkl')

    if os.path.exists(pickle_file):
        print('Loading existing data from pickle file...')
        with open(pickle_file, 'rb') as f:
            data = pickle.load(f)
            all_sequences = data['all_sequences']
            lengths = data['lengths']
            aa_counts = data['aa_counts']
            total_aa = data['total_aa']
            aa_freq = data['aa_freq']
    else:
        print('Extracting sequences from PDB files...')

        # 查找所有PDB文件
        pdb_files = glob.glob(os.path.join(pdb_dir, '**', '*.pdb'), recursive=True)
        print(f'Found {len(pdb_files)} PDB files')

        # 初始化PDB解析器
        parser = PDBParser()

        # 存储所有B链序列
        all_sequences = []

        # 提取B链序列
        for pdb_file in pdb_files:
            try:
                # 解析PDB文件
                structure = parser.get_structure(os.path.basename(pdb_file), pdb_file)

                # 获取B链
                for model in structure:
                    if 'B' in model:
                        chain = model['B']
                        sequence = []
                        for residue in chain:
                            if residue.id[0] == ' ':  # 只处理标准氨基酸
                                resname = residue.get_resname()
                                if resname in aa_3to1:
                                    sequence.append(aa_3to1[resname])

                        if sequence:
                            seq_str = ''.join(sequence)
                            all_sequences.append(seq_str)
                        break
            except Exception as e:
                continue

        print(f'Total sequences collected: {len(all_sequences)}')

        lengths = [len(seq) for seq in all_sequences]

        aa_counts = defaultdict(int)
        total_aa = 0
        for seq in all_sequences:
            for aa in seq:
                if aa in AMINO_ACIDS_SORTED:
                    aa_counts[aa] += 1
                    total_aa += 1

        aa_freq = {aa: aa_counts.get(aa, 0) / total_aa if total_aa > 0 else 0 for aa in AMINO_ACIDS_SORTED}

        data = {
            'all_sequences': all_sequences,
            'lengths': lengths,
            'aa_counts': aa_counts,
            'total_aa': total_aa,
            'aa_freq': aa_freq
        }

        with open(pickle_file, 'wb') as f:
            pickle.dump(data, f)
        print(f'Data saved to: {pickle_file}')

    # 绘制序列长度分布柱状图
    length_counter = Counter(lengths)
    sorted_lengths = sorted(length_counter.items())
    x_lengths = [item[0] for item in sorted_lengths]
    y_lengths = [item[1] for item in sorted_lengths]

    # 黄色渐变配色
    yellow_colors = ['#FFF59D', '#FFEE58', '#FFD54F', '#FFCA28', '#FFC107', '#FFB300', '#FFA000']
    colors_for_lengths = [yellow_colors[i % len(yellow_colors)] for i in range(len(x_lengths))]

    fig1, ax1 = plt.subplots(figsize=(12, 7))
    bars = ax1.bar(x_lengths, y_lengths, color=colors_for_lengths, edgecolor='#333333', linewidth=1.2)
    ax1.set_xlabel('Sequence Length', fontsize=14, fontweight='bold')
    ax1.set_ylabel('Number of Sequences', fontsize=14, fontweight='bold')
    ax1.set_title('Cyclic Peptide B-Chain Sequence Length Distribution', fontsize=16, fontweight='bold', pad=20)
    ax1.grid(axis='y', alpha=0.3, linestyle='--')
    ax1.set_axisbelow(True)

    for bar in bars:
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height)}',
                ha='center', va='bottom', fontsize=11, fontweight='medium')

    plt.tight_layout()
    output_pdf1 = os.path.join(output_dir, 'b_chain_length_distribution.pdf')
    plt.savefig(output_pdf1, dpi=300, bbox_inches='tight')
    print(f'Length distribution plot saved to: {output_pdf1}')
    plt.close()

    # 绘制氨基酸组成热图
    unique_lengths = sorted(list(set(lengths)), reverse=True)
    length_aa_freq = np.zeros((len(unique_lengths), len(AMINO_ACIDS_SORTED)))

    for seq in all_sequences:
        seq_len = len(seq)
        if seq_len in unique_lengths:
            len_idx = unique_lengths.index(seq_len)
            for aa in seq:
                if aa in AMINO_ACIDS_SORTED:
                    aa_idx = AMINO_ACIDS_SORTED.index(aa)
                    length_aa_freq[len_idx, aa_idx] += 1

    for i in range(len(unique_lengths)):
        row_sum = np.sum(length_aa_freq[i, :])
        if row_sum > 0:
            length_aa_freq[i, :] = length_aa_freq[i, :] / row_sum

    fig2, ax2 = plt.subplots(figsize=(16, 10))
    max_freq = np.max(length_aa_freq)

    # 黄色渐变色图
    custom_cmap = mcolors.LinearSegmentedColormap.from_list(
        'yellow_cmap',
        ['#FFF9C4', '#FFF59D', '#FFEE58', '#FFD54F', '#FFCA28', '#FFC107', '#FFB300', '#FFA000', '#FF8F00'],
        N=256
    )

    sns.heatmap(length_aa_freq, annot=False, cmap=custom_cmap,
                ax=ax2, cbar_kws={'label': 'Frequency'},
                xticklabels=AMINO_ACIDS_SORTED, yticklabels=unique_lengths,
                linewidths=0.5, vmin=0, vmax=max_freq * 0.7)
    ax2.set_xlabel('Amino Acid', fontsize=14, fontweight='bold')
    ax2.set_ylabel('Sequence Length', fontsize=14, fontweight='bold')
    ax2.set_title('B-Chain Amino Acid Composition by Sequence Length', fontsize=16, fontweight='bold', pad=20)
    plt.setp(ax2.get_xticklabels(), rotation=0, ha='center', fontsize=11)
    plt.setp(ax2.get_yticklabels(), rotation=0, fontsize=11)

    plt.tight_layout()
    output_pdf2 = os.path.join(output_dir, 'b_chain_aa_composition.pdf')
    plt.savefig(output_pdf2, dpi=300, bbox_inches='tight')
    print(f'Amino acid composition plot saved to: {output_pdf2}')
    plt.close()

    print('All plots generated successfully!')

if __name__ == '__main__':
    main()
