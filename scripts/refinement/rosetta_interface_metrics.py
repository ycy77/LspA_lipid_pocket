
# Shared path settings; scripts can be invoked from any working directory.
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "lib"))
from lsp_paths import WORKFLOW, DATA, ANALYSIS, RF_ROOT, env_path

import os
import csv
import multiprocessing
from multiprocessing import Pool
import pyrosetta
from pyrosetta import *
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.core.kinematics import MoveMap

# Initialize each worker process.
def init_worker():
    # Initialize PyRosetta separately in every process.
    pyrosetta.init()

# Calculate metrics for one structure.
def process_pdb(args):
    backbone, best_model, pdb_path = args
    try:
        # Load the structure.
        pose = pyrosetta.Pose()
        pyrosetta.pose_from_file(pose, pdb_path)

        # Create the score function.
        scorefxn = ScoreFunctionFactory.create_score_function("ref2015")

        # Create a MoveMap permitting backbone, side-chain and jump degrees of freedom.
        movemap = MoveMap()
        movemap.set_bb(True)  # Allow backbone motion.
        movemap.set_chi(True)  # Allow side-chain motion.
        movemap.set_jump(True)  # Allow inter-chain rigid-body motion.

        # Minimize energy.
        min_mover = MinMover()
        min_mover.movemap(movemap)
        min_mover.score_function(scorefxn)
        min_mover.min_type("dfpmin_armijo_nonmonotone")
        min_mover.tolerance(0.1)
        min_mover.max_iter(1000)
        min_mover.apply(pose)

        # Create an InterfaceAnalyzerMover for interface metrics.
        interface_mover = InterfaceAnalyzerMover()
        interface_mover.set_scorefunction(scorefxn)
        interface_mover.set_interface("A_B")  # Analyze the interface between chains A and B.
        interface_mover.apply(pose)

        # Get dSASA (change in solvent-accessible surface area at the interface).
        dsasa = interface_mover.get_interface_delta_sasa()

        # Get ddG (interface energy change).
        ddg = interface_mover.get_interface_dG()

        # Get BUH (buried unsatisfied hydrogen-bond count).
        buh = interface_mover.get_interface_delta_hbond_unsat()

        # Legacy CSV name CYCttE: this is the total ref2015 score, not a cycle-specific energy.
        cyctte = scorefxn(pose)

        return (backbone, best_model, dsasa, ddg, buh, cyctte, True, "")
    except Exception as e:
        error_msg = str(e)
        return (backbone, best_model, "N/A", "N/A", "N/A", "N/A", False, error_msg)

# Main function.
def main():
    # Read the summary TSV file.
    tsv_file = env_path("LSP_SUMMARY_TSV", f"{DATA}/cyclic_youhua_1/final_summary.tsv")
    pdb_files = []

    with open(tsv_file, "r", encoding='utf-8', errors='ignore') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            backbone = row['Backbone']
            best_model = row['Best_Model']
            # Construct the full path.
            pdb_path = os.path.join(env_path("LSP_PDB_ROOT", f"{DATA}/cyclic_youhua_1"), backbone, os.environ.get("LSP_MODEL_SUBDIR", "minimized"), best_model)
            if os.path.exists(pdb_path):
                pdb_files.append((backbone, best_model, pdb_path))

    print(f"Found {len(pdb_files)} PDB files")

    # Define the output file.
    output_file = env_path("LSP_OUTPUT_CSV", f"{ANALYSIS}/pyrosetta_youhua_results_minimized2.csv")
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    if not pdb_files:
        raise SystemExit("No structures matched the summary; check LSP_PDB_ROOT and LSP_MODEL_SUBDIR")

    # Write the header.
    with open(output_file, "w") as f:
        f.write("Backbone,Best_Model,dSASA,ddG,BUH,CYCttE\n")

    # Prepare worker arguments.
    params = [(backbone, best_model, pdb_path) for backbone, best_model, pdb_path in pdb_files]

    # Process structures with multiple workers.
    num_processes = max(1, min(multiprocessing.cpu_count(), int(os.environ.get("LSP_WORKERS", "4"))))  # Limit workers to LSP_WORKERS (default 4) and available CPUs.
    print(f"Using {num_processes} worker processes")

    with Pool(processes=num_processes, initializer=init_worker) as pool:
        results = pool.imap(process_pdb, params, chunksize=1)

        # Process results.
        for i, result in enumerate(results):
            backbone, best_model, dsasa, ddg, buh, cyctte, success, error_msg = result
            progress = (i + 1) / len(pdb_files) * 100

            if success:
                print(f"[{progress:.1f}%] Success: {backbone} - dSASA={dsasa:.2f}, ddG={ddg:.2f}, BUH={buh}, CYCttE={cyctte:.2f}")
                # Write successful results.
                with open(output_file, "a") as f:
                    f.write(f"{backbone},{best_model},{dsasa},{ddg},{buh},{cyctte}\n")
            else:
                print(f"[{progress:.1f}%] Failed: {backbone} - {error_msg}")
                # Write failed results.
                with open(output_file, "a") as f:
                    f.write(f"{backbone},{best_model},N/A,N/A,N/A,N/A\n")

    print(f"\nAnalysis completed; results saved to {output_file}")
    return output_file

if __name__ == "__main__":
    main()
