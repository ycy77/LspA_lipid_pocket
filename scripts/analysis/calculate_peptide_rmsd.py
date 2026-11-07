
# Shared path settings; scripts can be invoked from any working directory.
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "lib"))
from lsp_paths import WORKFLOW, DATA, ANALYSIS, RF_ROOT, env_path
_Path(ANALYSIS).mkdir(parents=True, exist_ok=True)

import os
import glob
import pickle
import copy
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.PDB import PDBParser, Superimposer
from collections import defaultdict

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.unicode_minus'] = False

def get_ca_atoms(structure, chain_id):
    atoms = []
    for model in structure:
        if chain_id in model:
            chain = model[chain_id]
            for residue in chain:
                if residue.id[0] == ' ' and 'CA' in residue:
                    atoms.append(residue['CA'])
    return atoms

def get_all_atoms_from_struct(structure):
    atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atoms.append(atom)
    return atoms

def calculate_rmsd_b_chain(initial_file, docked_file):
    parser = PDBParser(QUIET=True)
    try:
        initial_struct = parser.get_structure('initial', initial_file)
        docked_struct = parser.get_structure('docked', docked_file)

        initial_b_ca = get_ca_atoms(initial_struct, 'B')
        docked_b_ca = get_ca_atoms(docked_struct, 'B')

        if len(initial_b_ca) == 0 or len(docked_b_ca) == 0:
            return None, None

        b_length = len(initial_b_ca)

        min_len = min(len(initial_b_ca), len(docked_b_ca))
        initial_b_trim = initial_b_ca[:min_len]
        docked_b_trim = docked_b_ca[:min_len]

        docked_struct_copy = copy.deepcopy(docked_struct)
        docked_all_atoms_copy = get_all_atoms_from_struct(docked_struct_copy)

        sup = Superimposer()
        sup.set_atoms(initial_b_trim, docked_b_trim)
        sup.apply(docked_all_atoms_copy)

        initial_b_ca_final = get_ca_atoms(initial_struct, 'B')
        docked_b_ca_final = get_ca_atoms(docked_struct_copy, 'B')

        min_b_len = min(len(initial_b_ca_final), len(docked_b_ca_final))
        initial_b_trim_final = initial_b_ca_final[:min_b_len]
        docked_b_trim_final = docked_b_ca_final[:min_b_len]

        diff_sum = 0.0
        for a1, a2 in zip(initial_b_trim_final, docked_b_trim_final):
            diff = a1.get_coord() - a2.get_coord()
            diff_sum += np.sum(diff**2)

        rmsd = np.sqrt(diff_sum / len(initial_b_trim_final))

        return rmsd, b_length
    except Exception as e:
        return None, None

def check_cyclic_connected(struct):
    for model in struct:
        if 'B' in model:
            chain = model['B']
            residues = [r for r in chain if r.id[0] == ' ']
            if len(residues) < 2:
                return False

            # Check all adjacent peptide bonds
            for i in range(len(residues) - 1):
                current_res = residues[i]
                next_res = residues[i + 1]
                if 'C' in current_res and 'N' in next_res:
                    c_coord = current_res['C'].get_coord()
                    n_coord = next_res['N'].get_coord()
                    dist = np.linalg.norm(c_coord - n_coord)
                    # Normal peptide bond distance is ~1.32Å, allow up to 2.0Å
                    if dist > 2.0:
                        return False

            # Check N-terminal to C-terminal cyclic bond
            first_res = residues[0]
            last_res = residues[-1]
            if 'N' in first_res and 'C' in last_res:
                n_coord = first_res['N'].get_coord()
                c_coord = last_res['C'].get_coord()
                dist = np.linalg.norm(n_coord - c_coord)
                # Cyclic bond distance should be similar to peptide bond
                if dist > 2.0:
                    return False
            else:
                return False

            return True
    return False

def find_initial_file(initial_dir, base_name):
    subdir = os.path.join(initial_dir, base_name)
    if not os.path.exists(subdir):
        return None

    packed_dir = os.path.join(subdir, 'packed')
    if not os.path.exists(packed_dir):
        return None

    pdb_files = glob.glob(os.path.join(packed_dir, f'{base_name}_packed_*.pdb'))
    if not pdb_files:
        return None

    if len(pdb_files) != 1:
        raise ValueError(f'Ambiguous initial match ({len(pdb_files)}): {packed_dir}. Use the selected one-model-per-backbone directory.')
    return pdb_files[0]

def main():
    initial_dir = env_path('LSP_INITIAL_DIR', f'{DATA}/xiao_cyclic_1600_noRosettaFold/cyclic')
    docked_dir = env_path('LSP_DOCKED_DIR', f'{DATA}/xiao_cyclic_2_MPNN/best_structures')
    output_dir = f'{ANALYSIS}'
    pickle_file = os.path.join(output_dir, 'rmsd_data_b_chain_only_new.pkl')
    csv_file = os.path.join(output_dir, 'rmsd_results_b_chain_only_new.csv')

    if os.path.exists(pickle_file):
        os.remove(pickle_file)

    print('Calculating RMSD values (B-chain only alignment, cyclic filtered, threshold=1.5A)...')

    docked_files = glob.glob(os.path.join(docked_dir, '*_best.pdb'))
    print(f'Found {len(docked_files)} docked files')

    # Create a directory for intact structures.
    intact_dir = os.path.join(output_dir, 'C-N_cylic_intact')
    os.makedirs(intact_dir, exist_ok=True)

    # Create subdirectories for initial and docked structures.
    initial_intact_dir = os.path.join(intact_dir, 'initial')
    docked_intact_dir = os.path.join(intact_dir, 'docked')
    os.makedirs(initial_intact_dir, exist_ok=True)
    os.makedirs(docked_intact_dir, exist_ok=True)

    skipped_disconnected = 0
    rmsd_by_length = defaultdict(list)
    all_rmsds = []
    all_lengths = []
    all_results = []

    parser = PDBParser(QUIET=True)

    for idx, docked_file in enumerate(docked_files, 1):
        if idx % 10 == 0:
            print(f'Processing {idx}/{len(docked_files)}...')

        docked_basename = os.path.basename(docked_file)
        base_name = docked_basename.replace('_B_chain_emref_1_best.pdb', '')

        initial_file = find_initial_file(initial_dir, base_name)

        if initial_file is None:
            continue

        docked_struct = parser.get_structure('docked_check', docked_file)
        initial_struct = parser.get_structure('initial_check', initial_file)

        # Check that both initial and docked structures are intact.
        if check_cyclic_connected(docked_struct) and check_cyclic_connected(initial_struct):
            # Save the intact initial and docked structures.
            import shutil
            initial_basename = os.path.basename(initial_file)
            shutil.copy2(initial_file, os.path.join(initial_intact_dir, initial_basename))
            shutil.copy2(docked_file, os.path.join(docked_intact_dir, docked_basename))
        else:
            skipped_disconnected += 1
            continue

        rmsd, b_length = calculate_rmsd_b_chain(initial_file, docked_file)
        if rmsd is not None and b_length is not None:
            rmsd_by_length[b_length].append(rmsd)
            all_rmsds.append(rmsd)
            all_lengths.append(b_length)
            all_results.append({
                'Filename': docked_basename,
                'Base_Name': base_name,
                'B_Chain_Length': b_length,
                'RMSD': rmsd
            })

    print(f'Skipped {skipped_disconnected} disconnected cyclic peptides')
    print(f'Successfully calculated {len(all_rmsds)} RMSD values')
    print(f'Intact structures saved to: {intact_dir}')

    data = {
        'rmsd_by_length': dict(rmsd_by_length),
        'all_rmsds': all_rmsds,
        'all_lengths': all_lengths
    }

    with open(pickle_file, 'wb') as f:
        pickle.dump(data, f)
    print(f'Data saved to: {pickle_file}')

    df = pd.DataFrame(all_results)
    df = df.sort_values(['B_Chain_Length', 'RMSD'])
    df.to_csv(csv_file, index=False)
    print(f'CSV results saved to: {csv_file}')

    print(f'\nOverall RMSD Statistics (B-chain only alignment):')
    print(f'Min: {np.min(all_rmsds):.3f} Å')
    print(f'Max: {np.max(all_rmsds):.3f} Å')
    print(f'Mean: {np.mean(all_rmsds):.3f} Å')
    print(f'Median: {np.median(all_rmsds):.3f} Å')
    print(f'Std: {np.std(all_rmsds):.3f} Å')

    print(f'\nRMSD by B-chain length:')
    sorted_lengths = sorted(rmsd_by_length.keys())
    for length in sorted_lengths:
        rmsds = rmsd_by_length[length]
        print(f'  Length {length:2d}: {len(rmsds):3d} structures, Mean RMSD: {np.mean(rmsds):.3f} Å')

    fig, ax = plt.subplots(figsize=(14, 8))

    lengths_list = []
    mean_rmsds_list = []
    std_rmsds_list = []
    count_list = []

    for length in sorted_lengths:
        rmsds = rmsd_by_length[length]
        lengths_list.append(length)
        mean_rmsds_list.append(np.mean(rmsds))
        std_rmsds_list.append(np.std(rmsds))
        count_list.append(len(rmsds))

    morning_sea_colors = [
        '#05357d', '#166cae', '#3eaad0', '#9ccee7',
        '#f4f9ff', '#f6d353', '#ffaa2c', '#f78b15'
    ]
    color_dict = {length: morning_sea_colors[i % len(morning_sea_colors)]
                  for i, length in enumerate(sorted_lengths)}

    ax.fill_between(lengths_list,
                    np.array(mean_rmsds_list) - np.array(std_rmsds_list),
                    np.array(mean_rmsds_list) + np.array(std_rmsds_list),
                    alpha=0.2, color='#4A90D9', zorder=8)

    ax.plot(lengths_list, mean_rmsds_list, marker='o', markersize=10,
            linewidth=3, color='#1E5AA8', label='Mean RMSD', zorder=10)

    for length in sorted_lengths:
        rmsds = rmsd_by_length[length]
        x_jitter = np.random.normal(length, 0.12, size=len(rmsds))
        ax.scatter(x_jitter, rmsds, alpha=0.6, s=120,
                   color=color_dict[length], label=f'Length {length} (n={len(rmsds)})',
                   edgecolor='#333333', linewidth=1.5, zorder=5)

    for i, (length, count) in enumerate(zip(lengths_list, count_list)):
        pass

    ax.set_xlabel('B-Chain Length', fontsize=16, fontweight='bold')
    ax.set_ylabel('RMSD (Å)', fontsize=16, fontweight='bold')
    ax.set_title('B-Chain Conformational Change RMSD (Direct B-Chain Alignment)',
                 fontsize=18, fontweight='bold', pad=20)
    ax.grid(alpha=0.3, linestyle='--', linewidth=1.2)
    ax.set_axisbelow(True)
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=12)
    ax.set_xticks(lengths_list)
    ax.tick_params(axis='both', which='major', labelsize=12)

    plt.tight_layout()
    output_pdf = os.path.join(output_dir, 'b_chain_rmsd_conformational_change_new.pdf')
    plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
    print(f'\nPlot saved to: {output_pdf}')
    plt.close()

if __name__ == '__main__':
    main()
