#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose:
Batch preparation of RFdiffusion-generated cyclic-peptide backbones:
1. Replace chain A in each structure with chain A from the reference PDB.
2. Align to the reference and save as *_g.pdb.
3. Default output directory: 02_rf_diffusion_gai/cyclic/.
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

# ============ 1. Path configuration ============
### Base path.
BASE_DIR = f"{WORKFLOW}"
### Reference protein structure.
reference_pdb = os.environ.get("LSP_REFERENCE_PDB") or os.environ.get("LSP_TARGET_PDB")
if not reference_pdb:
    raise SystemExit("Set LSP_TARGET_PDB or LSP_REFERENCE_PDB to the receptor chain A PDB")
REF_PDB = env_path("LSP_REFERENCE_PDB", reference_pdb)
### Input directory.
INPUT_DIR = env_path("LSP_INPUT_DIR", f"{WORKFLOW}/01_rf_diffusion/cyclic")
### Output directory.
OUTPUT_DIR = env_path("LSP_OUTPUT_DIR", f"{WORKFLOW}/02_rf_diffusion_gai/cyclic")
LOG_FILE = os.path.join(OUTPUT_DIR, "pymol_convert.log")

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============ 2. Initialize PyMOL ============
pymol.finish_launching(['pymol', '-cq'])  # Run headless (-c) and quiet (-q).

# ============ 3. Scan input files ============
pdb_files = sorted(glob.glob(f"{INPUT_DIR}/*.pdb"))
print(f"[INFO] Found {len(pdb_files)} input files")
if not pdb_files:
    raise SystemExit(f"No PDB inputs: {INPUT_DIR}")

# ============ 4. Main loop ============
with open(LOG_FILE, "a") as log:
    for pdb_path in pdb_files:
        pdb_name = os.path.basename(pdb_path)
        output_path = os.path.join(OUTPUT_DIR, pdb_name.replace(".pdb", "_g.pdb"))

        # Skip existing output files.
        if os.path.exists(output_path):
            print(f"[SKIP] Already exists: {output_path}")
            continue

        try:
            # Load reference and target structures.
            cmd.reinitialize()  # Clear the session.
            cmd.load(REF_PDB, "ref")
            cmd.load(pdb_path, "complex")

            # Align chain A.
            cmd.align("complex and chain A", "ref and chain A")

            # Remove the original chain A from the complex.
            cmd.remove("complex and chain A")

            # Copy reference chain A into the complex.
            cmd.create("complex", "complex or (ref and chain A)")

            # Save the structure.
            cmd.save(output_path, "complex")

            log.write(f"Success {pdb_name}\n")
            print(f"[OK] {pdb_name}")

        except Exception as e:
            log.write(f"Failed {pdb_name}: {e}\n")
            print(f"[ERROR] {pdb_name}: {e}")

# ============ 5. Cleanup ============
cmd.quit()
print("[INFO] Processing completed")
