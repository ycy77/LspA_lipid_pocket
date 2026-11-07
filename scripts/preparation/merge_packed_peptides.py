#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose:
Batch alignment and replacement of chain B using side-chain-packed PDB files:
1. Align chain B in each packed PDB to chain B in the original PDB.
2. Replace the original chain B with the aligned packed chain B.
3. Save results in each sample's merge_pdb directory.
"""


# Shared path settings; scripts can be invoked from any working directory.
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "lib"))
from lsp_paths import WORKFLOW, DATA, ANALYSIS, RF_ROOT, env_path

import os
import glob
import pymol
from pymol import cmd
import logging

# ============ 1. Path configuration ============
BASE_DIR = f"{WORKFLOW}"
# Original PDB directory containing reference chain B structures.
ORIGINAL_PDB_DIR = env_path("LSP_REFERENCE_DIR", f"{WORKFLOW}/02_rf_diffusion_gai/cyclic_refine")
# Root directory of side-chain-packed PDB files.
PACKED_ROOT = env_path("LSP_INPUT_DIR", f"{WORKFLOW}/03_MPNN/cyclic_refine")
# Root output directory.
DEST_ROOT = env_path("LSP_OUTPUT_DIR", f"{WORKFLOW}/04_select_MPNN/cyclic_refine")
# Log file path.
LOG_FILE = f"{DEST_ROOT}/pymol_align.log"

# Create the output root before initializing file logging.
os.makedirs(DEST_ROOT, exist_ok=True)

# Log successes and failures.
logging.basicConfig(
    filename=LOG_FILE,
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

# ============ 2. Initialize PyMOL in quiet, headless mode ============
pymol.finish_launching(['pymol', '-cq'])  # Run without a GUI and suppress startup output.

# ============ 3. Batch processing function ============
def process_sample(sample_name, sample_path):
    """Align and replace chain B for one sample."""
    # Sample-specific paths.
    # The original PDB basename should match the sample prefix, e.g. cyclic_T100_018_0_g.
    base_sample_name = sample_name.split('_B_')[0]
    original_pdb = os.path.join(ORIGINAL_PDB_DIR, f"{base_sample_name}.pdb")  # Original PDB containing chain B.
    packed_dir = os.path.join(sample_path, "packed")         # Side-chain-packed PDB directory.
    merge_dir = os.path.join(DEST_ROOT, sample_name, "merge_pdb")       # Output directory.
    os.makedirs(merge_dir, exist_ok=True)                                # Create the output directory.

    # Check that the original PDB exists.
    if not os.path.exists(original_pdb):
        msg = f"Sample {sample_name}: original PDB missing: {original_pdb}"
        print(f"[ERROR] {msg}")
        logging.error(msg)
        return

    # Find all side-chain-packed PDB files for this sample.
    packed_pdbs = glob.glob(os.path.join(packed_dir, f"*.pdb"))
    if not packed_pdbs:
        msg = f"Sample {sample_name}: no side-chain-packed PDB files found in {packed_dir}"
        print(f"[WARNING] {msg}")
        logging.warning(msg)
        return

    # Process each side-chain-packed PDB.
    for packed_pdb in packed_pdbs:
        packed_basename = os.path.basename(packed_pdb)
        output_pdb = os.path.join(merge_dir, packed_basename)  # Preserve the packed PDB basename in the output.

        # Skip files already processed.
        if os.path.exists(output_pdb):
            print(f"[SKIP] Aligned output already exists: {output_pdb}")
            continue

        try:
            # Clear the PyMOL session.
            cmd.reinitialize()

            # Load the original reference and the packed structure to align.
            cmd.load(original_pdb, "original")    # Original structure: retain other chains and use chain B as the reference.
            cmd.load(packed_pdb, "packed")        # Packed structure containing chain B to align.

            # Align packed chain B to original chain B.
            # Use chain B CA atoms for backbone alignment.
            alignment = cmd.align(
                "packed and chain B and name CA",  # Mobile chain: packed chain B.
                "original and chain B and name CA" # Reference chain: original chain B.
            )
            print(f"[ALIGN] {sample_name} - Alignment RMSD: {alignment[0]:.2f}")

            # Remove original chain B before replacement.
            cmd.remove("original and chain B")

            # Merge original non-B chains with packed chain B.
            cmd.create("merged", "original or (packed and chain B)")

            # Save the merged structure.
            cmd.save(output_pdb, "merged")

            # Log successful processing.
            msg = f"Sample {sample_name} processed successfully: {packed_basename}"
            print(f"[OK] {msg}")
            logging.info(msg)

        except Exception as e:
            # Log errors.
            msg = f"Sample {sample_name}: failed to process {packed_basename}: {str(e)}"
            print(f"[ERROR] {msg}")
            logging.error(msg)

# ============ 4. Main loop over samples ============
if __name__ == "__main__":
    # Find sample directories, e.g. cyclic_T100_018_0_g_B_chain_emref_1.
    sample_dirs = [
        d for d in glob.glob(os.path.join(PACKED_ROOT, "*"))
        if os.path.isdir(d) and "cyclic_" in os.path.basename(d)
    ]

    total = len(sample_dirs)
    print(f"[INFO] Found {total} sample directories; starting processing...")
    logging.info(f"Starting processing of {total} samples")

    for i, sample_dir in enumerate(sample_dirs, 1):
        sample_name = os.path.basename(sample_dir)
        print(f"\n===== Processing sample {i}/{total}: {sample_name} =====")
        process_sample(sample_name, sample_dir)

    # Quit PyMOL.
    cmd.quit()
    print("\n[INFO] Processing completed! Results saved in:", DEST_ROOT)
    logging.info("Processing completed")
