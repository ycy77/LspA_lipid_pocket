
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

# Map three-letter amino acid codes to one-letter codes.
aa_3to1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

# User-specified color palette.
simple_light_colors = {
    # Hydrophobic (nonpolar).
    'A': '#cb5148', 'V': '#cb5148', 'L': '#cb5148', 'I': '#cb5148','Y': '#cb5148',
    'M': '#cb5148', 'F': '#cb5148', 'W': '#cb5148', 'P': '#cb5148','G': '#cb5148',
    # Polar (uncharged).
    'S': "#0cd26f", 'T': '#0cd26f', 'N': '#0cd26f', 'Q': '#0cd26f', 'C': '#0cd26f',
    # Negatively charged.
    'D': "#00b7ef", 'E': '#00b7ef',
    # Positively charged.
    'K': "#2d22c7", 'R': '#2d22c7', 'H': '#2d22c7'
}

# PDB input directory.
pdb_dir = env_path('LSP_INPUT_DIR', f'{DATA}/cyclic_youhua_1')

# Find all PDB files.
pdb_files = glob.glob(os.path.join(pdb_dir, '**', '*.pdb'), recursive=True)
print(f"Found {len(pdb_files)} PDB files")

# Initialize the PDB parser.
parser = PDBParser()

# Store all chain B sequences.
b_chain_sequences = defaultdict(list)

# Extract chain B sequences.
for pdb_file in pdb_files:
    try:
        # Parse the PDB file.
        structure = parser.get_structure(os.path.basename(pdb_file), pdb_file)

        # Get chain B.
        for model in structure:
            if 'B' in model:
                chain = model['B']
                sequence = []
                for residue in chain:
                    if residue.id[0] == ' ':  # Process standard amino acids only.
                        resname = residue.get_resname()
                        if resname in aa_3to1:
                            sequence.append(aa_3to1[resname])

                if sequence:
                    seq_str = ''.join(sequence)
                    b_chain_sequences[len(seq_str)].append(seq_str)
                break
    except Exception as e:
        continue

# Summarize sequences by length.
print("\nSequence counts by length:")
length_counts = {}
for length, sequences in b_chain_sequences.items():
    length_counts[length] = len(sequences)
    print(f"Length {length}: {len(sequences)} sequences")

# Save sequence information to CSV.
all_sequences = []
for length, sequences in b_chain_sequences.items():
    for seq in sequences:
        all_sequences.append({'length': length, 'sequence': seq})

df = pd.DataFrame(all_sequences)
df.to_csv(f'{ANALYSIS}/b_chain_sequences.csv', index=False)
print(f"\nSequence information saved to: {ANALYSIS}/b_chain_sequences.csv")

# Plot a sequence logo for each length group.
for length in sorted(b_chain_sequences.keys()):
    sequences = b_chain_sequences[length]
    print(f"\nPlotting logo for length {length}, using {len(sequences)} sequences")

    if len(sequences) < 2:
        print(f"Insufficient sequences for length {length}; skipping")
        continue

    # Calculate amino acid frequencies at each position.
    max_len = max(len(seq) for seq in sequences)

    # Amino acid order.
    aa_order = sorted(set(''.join(sequences)))

    # Count occurrences of each amino acid at each position.
    counts = defaultdict(lambda: [0] * max_len)
    for seq in sequences:
        for i, aa in enumerate(seq):
            counts[aa][i] += 1

    # Calculate the frequency matrix.
    freq_data = {}
    for aa in aa_order:
        freq_data[aa] = [count / len(sequences) for count in counts[aa]]

    # Create the DataFrame for logomaker (positions as rows, amino acids as columns).
    df_logo = pd.DataFrame(freq_data)

    # Plot the sequence logo.
    fig, ax = plt.subplots(figsize=(max_len * 0.5, 3))

    try:
        # Draw with logomaker using the specified color palette.
        logo = logomaker.Logo(df_logo, ax=ax, color_scheme=simple_light_colors)

        # Set titles and labels.
        ax.set_title(f'Amino Acid Logo - Length {length} (n={len(sequences)})', fontsize=14, fontweight='bold')
        ax.set_xlabel('Position', fontsize=12)
        ax.set_ylabel('Frequency', fontsize=12)

        # Save the figure.
        output_path = f'{ANALYSIS}/logo_length_{length}.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {output_path}")

        output_path_pdf = f'{ANALYSIS}/logo_length_{length}.pdf'
        plt.savefig(output_path_pdf, dpi=300, bbox_inches='tight', format='pdf')
        print(f"Saved: {output_path_pdf}")
    except Exception as e:
        print(f"Error plotting sequence logo: {e}")
    finally:
        plt.close()

print("\nAll sequence logos completed!")
