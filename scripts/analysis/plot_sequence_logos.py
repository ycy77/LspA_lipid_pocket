
# Shared path settings; scripts can be invoked from any working directory.
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "lib"))
from lsp_paths import WORKFLOW, DATA, ANALYSIS, RF_ROOT, env_path
_Path(ANALYSIS).mkdir(parents=True, exist_ok=True)

import os
import glob
from collections import defaultdict
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import logomaker
from Bio.PDB import PDBParser

# 氨基酸三字母到单字母的映射
aa_3to1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

# 用户提供的配色方案
simple_light_colors = {
    # 疏水（非极性）
    'A': '#cb5148', 'V': '#cb5148', 'L': '#cb5148', 'I': '#cb5148','Y': '#cb5148',
    'M': '#cb5148', 'F': '#cb5148', 'W': '#cb5148', 'P': '#cb5148','G': '#cb5148',
    # 极性（不带电）
    'S': "#0cd26f", 'T': '#0cd26f', 'N': '#0cd26f', 'Q': '#0cd26f', 'C': '#0cd26f',
    # 带电荷（负电）
    'D': "#00b7ef", 'E': '#00b7ef',
    # 带电荷（正电）
    'K': "#2d22c7", 'R': '#2d22c7', 'H': '#2d22c7'
}

# PDB文件目录
pdb_dir = env_path('LSP_INPUT_DIR', f'{DATA}/cyclic_youhua_1')

# 查找所有PDB文件
pdb_files = glob.glob(os.path.join(pdb_dir, '**', '*.pdb'), recursive=True)
print(f"找到 {len(pdb_files)} 个PDB文件")

# 初始化PDB解析器
parser = PDBParser()

# 存储所有B链序列
b_chain_sequences = defaultdict(list)

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
                    b_chain_sequences[len(seq_str)].append(seq_str)
                break
    except Exception as e:
        continue

# 按长度分组统计
print("\n按长度分组统计:")
length_counts = {}
for length, sequences in b_chain_sequences.items():
    length_counts[length] = len(sequences)
    print(f"长度 {length}: {len(sequences)} 个序列")

# 保存序列信息到CSV
all_sequences = []
for length, sequences in b_chain_sequences.items():
    for seq in sequences:
        all_sequences.append({'length': length, 'sequence': seq})

df = pd.DataFrame(all_sequences)
df.to_csv(f'{ANALYSIS}/b_chain_sequences.csv', index=False)
print(f"\n已保存序列信息到: {ANALYSIS}/b_chain_sequences.csv")

# 为每个长度组绘制LOGO图
for length in sorted(b_chain_sequences.keys()):
    sequences = b_chain_sequences[length]
    print(f"\n绘制长度 {length} 的LOGO图，共 {len(sequences)} 个序列")

    if len(sequences) < 2:
        print(f"长度 {length} 的序列数量不足，跳过")
        continue

    # 计算每个位置的氨基酸频率
    max_len = max(len(seq) for seq in sequences)

    # 氨基酸顺序
    aa_order = sorted(set(''.join(sequences)))

    # 统计每个位置的氨基酸出现次数
    counts = defaultdict(lambda: [0] * max_len)
    for seq in sequences:
        for i, aa in enumerate(seq):
            counts[aa][i] += 1

    # 计算频率矩阵
    freq_data = {}
    for aa in aa_order:
        freq_data[aa] = [count / len(sequences) for count in counts[aa]]

    # 创建DataFrame（logomaker要求列是位置，行是氨基酸）
    df_logo = pd.DataFrame(freq_data)

    # 绘制LOGO图
    fig, ax = plt.subplots(figsize=(max_len * 0.5, 3))

    try:
        # 使用logomaker绘制，应用用户提供的配色方案
        logo = logomaker.Logo(df_logo, ax=ax, color_scheme=simple_light_colors)

        # 设置标题和标签
        ax.set_title(f'Amino Acid Logo - Length {length} (n={len(sequences)})', fontsize=14, fontweight='bold')
        ax.set_xlabel('Position', fontsize=12)
        ax.set_ylabel('Frequency', fontsize=12)

        # 保存图片
        output_path = f'{ANALYSIS}/logo_length_{length}.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"已保存: {output_path}")

        output_path_pdf = f'{ANALYSIS}/logo_length_{length}.pdf'
        plt.savefig(output_path_pdf, dpi=300, bbox_inches='tight', format='pdf')
        print(f"已保存: {output_path_pdf}")
    except Exception as e:
        print(f"绘制LOGO图时出错: {e}")
    finally:
        plt.close()

print("\n所有LOGO图绘制完成！")
