#!/usr/bin/env bash
set -euo pipefail
LSP_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
source "$LSP_ROOT/config/paths.sh"
export LC_ALL=C.UTF-8  # Use UTF-8 terminal output.

# --------------------------
# Batch cyclic-peptide docking configuration.
# --------------------------
CONDA_ENV="${LSP_HADDOCK_ENV:-}"                             # Conda environment containing haddock3.
TARGET_PDB="${LSP_TARGET_PDB:?Set LSP_TARGET_PDB to the receptor chain A PDB}"  # Receptor chain A PDB.
B_CHAINS_FOLDER="${LSP_INPUT_DIR:-${LSP_WORKFLOW_DIR}/05_peptide_chains/cyclic}"            # Directory containing cyclic-peptide chain B PDB files.
HOTSPOT_LIST=(129 130 132 136 137 138 140)  # Hotspot residues on chain A.
OUTPUT_ROOT="${LSP_OUTPUT_DIR:-${LSP_DATA_DIR}/xiao_haddock3_cyclic_batch}"       # Root directory for batch outputs.
CPU_PER_TASK="${LSP_CPUS_PER_TASK:-4}"                              # Cores per task (at least 4 recommended for cyclic peptides).
MAX_CONCURRENT="${LSP_MAX_CONCURRENT:-1}"                             # Concurrent tasks; total cores = MAX_CONCURRENT * CPU_PER_TASK.
[[ "$CPU_PER_TASK" =~ ^[1-9][0-9]*$ && "$MAX_CONCURRENT" =~ ^[1-9][0-9]*$ ]] || { echo 'CPU/concurrency must be positive integers' >&2; exit 2; }
SLEEP_WHEN_FULL=2                          # Seconds to wait when all task slots are occupied.
# Cyclic-peptide parameters for preserving ring connectivity.
CLUST_CUTOFF=0.8                            # Clustering cutoff for cyclic-peptide structural similarity.
TOPOAA_CYCLIC_PEPT="true"                   # Enable cyclic-peptide recognition.
TOPOAA_CYCLIC_DIST="1.8"                    # Peptide-bond detection cutoff (1.8 angstrom).
TOPOAA_DISULPHIDE_DIST="3.0"                # Disulfide detection cutoff; use the British spelling required by HADDOCK.
FLEXREF_W_VDW="0.6"                         # Reduced van der Waals weight for cyclic-peptide refinement.
# --------------------------

# --------------------------
# Initialize task counters and the global score file.
# --------------------------
success_count=0
fail_count=0
total_task=0
GLOBAL_SCORE_FILE="${OUTPUT_ROOT}/global_score_summary.tsv"  # Global score summary.


# --------------------------
# Check the environment.
# --------------------------
# Use the active environment, or an explicitly requested Conda environment.
if [[ -n "$CONDA_ENV" ]]; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV"
fi
if [[ "${LSP_PREPARE_ONLY:-0}" != 1 ]]; then
    command -v haddock3 >/dev/null || { echo "Activate a HADDOCK3 environment first" >&2; exit 1; }
fi
# --------------------------
# Validate batch input files and directories.
# --------------------------
# 1. Check the receptor.
[ -f "${TARGET_PDB}" ] || {
    echo "ERROR Receptor not found: ${TARGET_PDB}"
    exit 1
}

# 2. Check the chain B directory and file count.
[ -d "${B_CHAINS_FOLDER}" ] || {
    echo "ERROR Chain B directory not found: ${B_CHAINS_FOLDER}"
    exit 1
}
shopt -s nullglob
pdb_list=("${B_CHAINS_FOLDER}"/*.pdb)
total_task=${#pdb_list[@]}
shopt -u nullglob
[ "${total_task}" -eq 0 ] && {
    echo "ERROR No .pdb files in ${B_CHAINS_FOLDER}!"
    exit 1
}

# 3. Check resources to avoid oversubscription.
if command -v nproc >/dev/null 2>&1; then
    avail_cpu="$(nproc)"
    need_cpu=$((MAX_CONCURRENT * CPU_PER_TASK))
    if [ "${need_cpu}" -gt "${avail_cpu}" ]; then
        echo "WARNING: ${need_cpu} cores requested, but only ${avail_cpu} available; performance may suffer!"
        exit 1
    fi
fi

# 4. Create the output directory and global score header.
mkdir -p "${OUTPUT_ROOT}"
real_output_root="$(cd -- "$OUTPUT_ROOT" && pwd)"
[[ -f "$GLOBAL_SCORE_FILE" ]] || printf "Model\tHADDOCK_Score\tStatus\n" > "$GLOBAL_SCORE_FILE"  # Header.

# Initialize logging.
echo "========================================"
echo " Starting batch cyclic-peptide docking with ring-preserving restraints"
echo "========================================"
echo " Tasks: ${total_task} peptides | Concurrent: ${MAX_CONCURRENT} | Cores per task: ${CPU_PER_TASK}"
echo " Output root: ${real_output_root}"
echo " Cyclic-peptide settings: covalent ring topology, fixed terminal residues, vdW weight 0.6"
echo "========================================"

# --------------------------
# Dock one cyclic peptide and organize its results.
# --------------------------
run_cyclic_batch() {
    local b_pdb_path="$1"
    local b_basename="$(basename -- "$b_pdb_path" .pdb)"  # Ligand basename, e.g. cyclic_T100_002.
    local task_dir="${real_output_root}/${b_basename}"   # Directory for this task.
    local log_path="${task_dir}/haddock3.log"
    local status_file="${task_dir}/.status"
    local task_tag="[$b_basename]"
    local b_chain_max_residue  # Maximum cyclic-peptide residue number.
    local emref1_pdb=""        # Path to emref_1.pdb.
    local emref1_score="NA"    # Score of emref_1.

    # --------------------------
    # 1. Validate the cyclic-peptide PDB.
    # --------------------------
    # Check chain B ATOM records for empty or invalid input.
    local b_atom_count=$(awk '$1 == "ATOM" && $5 == "B" {print $0}' "${b_pdb_path}" | wc -l)
    if [ "${b_atom_count}" -lt 10 ]; then
        echo "${task_tag} ERROR Too few ATOM records (${b_atom_count}); skipping!"
        echo -e "${b_basename}\t${emref1_score}\tPDB_ERROR" >> "${GLOBAL_SCORE_FILE}"
        fail_count=$((fail_count + 1))
        return 0
    fi

    # Count peptide residues, excluding insertion codes.
    echo "${task_tag} Counting cyclic-peptide residues..."
    b_chain_max_residue=$(awk '
        $1 == "ATOM" && $5 == "B" {
            clean_resid = $6; gsub(/[^0-9]/, "", clean_resid);
            if (clean_resid ~ /^[0-9]+$/ && clean_resid > 0) print clean_resid
        }
    ' "${b_pdb_path}" | sort -n | tail -n 1)
    # Require at least three residues in the cyclic peptide.
    if [ -z "${b_chain_max_residue}" ] || [ "${b_chain_max_residue}" -lt 3 ]; then
        echo "${task_tag} ERROR Invalid residue count (<3); skipping!"
        echo -e "${b_basename}\t${emref1_score}\tRESID_ERROR" >> "${GLOBAL_SCORE_FILE}"
        fail_count=$((fail_count + 1))
        return 0
    fi
    echo "${task_tag} OK Cyclic-peptide residues: 1 to ${b_chain_max_residue} (${b_chain_max_residue} total)"

    # --------------------------
    # 2. Prepare the task directory.
    # --------------------------
    if [ -d "${task_dir}" ]; then
        echo "${task_tag} SKIP existing directory (use a new output directory to rerun)."
        return 0
    fi
    mkdir -p "${task_dir}/data" || {
        echo "${task_tag} ERROR Cannot create task directory"
        echo -e "${b_basename}\t${emref1_score}\tDIR_ERROR" >> "${GLOBAL_SCORE_FILE}"
        fail_count=$((fail_count + 1))
        return 0
    }

    # --------------------------
    # 3. Copy required inputs while preserving originals.
    # --------------------------
    cp -- "${b_pdb_path}" "${task_dir}/ligand_B.pdb" || {
        echo "${task_tag} ERROR Failed to copy cyclic peptide"
        echo -e "${b_basename}\t${emref1_score}\tCOPY_B_ERROR" >> "${GLOBAL_SCORE_FILE}"
        fail_count=$((fail_count + 1))
        return 0
    }
    cp -- "${TARGET_PDB}" "${task_dir}/receptor_A.pdb" || {
        echo "${task_tag} ERROR Failed to copy receptor"
        echo -e "${b_basename}\t${emref1_score}\tCOPY_TARGET_ERROR" >> "${GLOBAL_SCORE_FILE}"
        fail_count=$((fail_count + 1))
        return 0
    }

    # --------------------------
    # 4. Generate ring-closure and hotspot restraints.
    # --------------------------
    local ambig_tbl="${task_dir}/ambig.tbl"
    cat > "${ambig_tbl}" <<EOF
# Cyclic-peptide docking restraints for ring closure and binding hotspots.
# 1. Ring closure: first-residue N to last-residue C (1.36 +/- 0.05 angstrom).
assign (resid 1 and segid B and name N) (resid ${b_chain_max_residue} and segid B and name C) 1.36 0.05 0.05
# 2. Hotspots: chain A hotspots to chain B heavy atoms (2.0 +/- 1.0 angstrom).
EOF
    # Add hotspot restraints.
    for resid in "${HOTSPOT_LIST[@]}"; do
        echo "assign (resid ${resid} and segid A and not hydrogen) (segid B and not hydrogen) 2.0 1.0 3.0" >> "${ambig_tbl}"
    done
    echo "${task_tag} OK Restraints written: ${ambig_tbl} (1 ring closure + ${#HOTSPOT_LIST[@]} hotspots)"

    # --------------------------
    # 5. Generate the cyclic-peptide config.toml.
    # --------------------------
    local flex_start=2
    local flex_end=$((b_chain_max_residue - 1))  # Exclude terminal residues from the flexible region.
    cat > "${task_dir}/config.toml" <<EOF
# HADDOCK3 cyclic-peptide docking with ring-preserving restraints.
run_dir = "./"
ncores = ${CPU_PER_TASK}
molecules = ["receptor_A.pdb", "ligand_B.pdb"]

# Stage 1: topoaa (generate covalent cyclic topology).
[topoaa]
autohis = false
cyclicpept = ${TOPOAA_CYCLIC_PEPT}          # Global cyclic-peptide switch.
cyclicpept_dist = ${TOPOAA_CYCLIC_DIST}      # Peptide-bond detection cutoff.
disulphide_dist = ${TOPOAA_DISULPHIDE_DIST}  # Disulfide cutoff (British spelling).
# Terminal settings for noncyclic chain A.
mol1.charged_nter = false
mol1.charged_cter = false
mol1.5_phosphate = false
# Cyclic-peptide settings for chain B.
mol2.charged_nter = false                    # No free N terminus.
mol2.charged_cter = false                    # No free C terminus.
mol2.5_phosphate = false
mol2.cyclicpept = true                       # Explicitly mark chain B as cyclic.

# Stage 2: rigidbody (sample binding orientations).
[rigidbody]
mode = "local"
sampling = 300                             # Sample diverse conformations.
ambig_fname = "ambig.tbl"
iniseed = 12345                              # Fixed seed for reproducibility.
w_vdw = 0.8                                  # Reduce the van der Waals weight.
w_elec = 0.2
w_air = 2.0                                  # Increase the hotspot restraint weight.

# Stage 3: clustfcc (cluster rigid-body poses).
[clustfcc]
clust_cutoff = ${CLUST_CUTOFF}
min_population = 5                           # Filter isolated conformations.
contact_distance_cutoff = 8.0
plot_matrix = false

# Stage 4: seletopclusts (select poses after rigid-body docking).
[seletopclusts]
top_clusters = 10
top_models = 1
sortby = "score"                             # Select by score.

# Stage 5: flexref (flexible refinement with ring restraints).
[flexref]
mode = "local"
max_nmodels = 10
ambig_fname = "ambig.tbl"
sampling_factor = 1
tolerance = 1
w_vdw = ${FLEXREF_W_VDW}                     # Use a van der Waals weight of 0.6.
w_elec = 0.2
w_air = 2.5                                  # Increase restraint weight to preserve ring closure.
fle_sta_1 = ${flex_start}                    # Flexible region starts at residue 2.
fle_end_1 = ${flex_end}                      # Flexible region ends before the last residue.
fle_seg_1 = "B"                              # Only chain B is flexible.

# Stage 6: clustfcc (cluster flexibly refined poses).
[clustfcc]
clust_cutoff = ${CLUST_CUTOFF}
min_population = 2
contact_distance_cutoff = 8.0
plot_matrix = false

# Stage 7: seletopclusts (select refined poses).
[seletopclusts]
top_clusters = 5
top_models = 1
sortby = "score"

# Stage 8: emref (energy minimization).
[emref]
max_nmodels = 5
ambig_fname = "ambig.tbl"
mode = "local"
tolerance = 2
fle_sta_1 = ${flex_start}                    # Use the same flexible region as flexref.
fle_end_1 = ${flex_end}
fle_seg_1 = "B"
w_vdw = ${FLEXREF_W_VDW}                     # Keep the van der Waals weight consistent.
w_elec = 0.2
w_air = 2.5

# Stage 9: capri (compute metrics and rank conformations).
[caprieval]
ncores = ${CPU_PER_TASK}          # Match the core count to the per-task allocation.
max_cpus = ${CPU_PER_TASK}        # Maximum CPU count equals ncores.
sortby = "score"
sort_ascending = true             # Sort ascending: lower scores rank first.
clt_threshold = 1           # Default clustering threshold.
EOF

    if [[ "${LSP_PREPARE_ONLY:-0}" == 1 ]]; then
        printf 'PREPARED\n' > "$status_file"
        echo "${task_tag} Prepared configuration; no docking launched."
        return 0
    fi

    # --------------------------
    # 6. Run the docking task.
    # --------------------------
    echo "${task_tag} Starting docking (log: ${log_path})"
    (
        cd "${task_dir}" || exit 1
        export OMP_NUM_THREADS="${CPU_PER_TASK}"  # Limit OpenMP threads to the per-task core allocation.
        if haddock3 config.toml --restart 0 > "${log_path}" 2>&1; then
            printf "SUCCESS\n0\n" > "${status_file}"
        else
            printf "FAIL\nhaddock_failed\n" > "${status_file}"
        fi
    )

    # --------------------------
    # 7. Postprocess: decompress emref_1.pdb and extract its score.
    # --------------------------
    if [ -f "${status_file}" ] && [ "$(head -n1 "${status_file}")" = "SUCCESS" ]; then
        # 7.1. Decompress emref_1.pdb.gz into the shared output directory.
        local emref_gz="${task_dir}/7_emref/emref_1.pdb.gz"
        # Shared output directory.
        BEST_PDB_DIR="${OUTPUT_ROOT}/best_pdb_all"
        mkdir -p "${BEST_PDB_DIR}"

        if [ -f "${emref_gz}" ]; then
            # Decompress and rename as ligand_basename_emref_1.pdb.
            gunzip -c "${emref_gz}" > "${BEST_PDB_DIR}/${b_basename}_emref_1.pdb"
            emref1_pdb="${BEST_PDB_DIR}/${b_basename}_emref_1.pdb"
            echo "${task_tag} OK Best PDB extracted to: ${emref1_pdb}"
        else
            echo "${task_tag} WARNING emref_1.pdb.gz not found!"
        fi

        # 7.2. Extract the emref_1 score from capri_ss.tsv.
        local capri_tsv="${task_dir}/8_caprieval/capri_ss.tsv"
        if [ -f "${capri_tsv}" ]; then
            # Match the emref_1.pdb row with awk and extract column 4 (score).
            emref1_score=$(awk -F'\t' '$1 ~ /emref_1.pdb/ {print $4; exit}' "${capri_tsv}")
            # Validate that the score is numeric.
            if [[ "${emref1_score}" =~ ^-?[0-9]+\.[0-9]+$ ]]; then
                echo "${task_tag} OK emref_1 score：${emref1_score}"
            else
                emref1_score="NA"
                echo "${task_tag} WARNING Score extraction failed (nonnumeric value)"
            fi
        else
            echo "${task_tag} WARNING capri_ss.tsv not found!"
        fi

        # 7.3. Mark success and append the global score.
        success_count=$((success_count + 1))
        echo -e "${b_basename}\t${emref1_score}\tSUCCESS" >> "${GLOBAL_SCORE_FILE}"
        echo "${task_tag} Task completed successfully!"

    else
        # Record failure status.
        fail_count=$((fail_count + 1))
        echo -e "${b_basename}\t${emref1_score}\tFAIL" >> "${GLOBAL_SCORE_FILE}"
        echo "${task_tag} ERROR Task failed (see log: ${log_path})"
    fi
}

# --------------------------
# Main loop: batch scheduling and concurrency control.
# --------------------------
echo " Submitting ${total_task} cyclic-peptide tasks..."
pids=()
for bpath in "${pdb_list[@]}"; do
    # Wait while active tasks reach MAX_CONCURRENT.
    while true; do
        current_jobs=$(jobs -rp | wc -l)
        if [ "${current_jobs}" -lt "${MAX_CONCURRENT}" ]; then
            break
        fi
        echo " Active tasks: ${current_jobs}/${MAX_CONCURRENT}; waiting ${SLEEP_WHEN_FULL} seconds..."
        sleep "${SLEEP_WHEN_FULL}"
    done

    # Launch the task in the background.
    run_cyclic_batch "${bpath}" &
    pids+=("$!")
    :
done

# --------------------------
# Wait for all tasks and report final counts.
# --------------------------
echo "-----------------------------------------"
echo " All tasks submitted; waiting for background tasks to finish..."
child_failed=0
for pid in "${pids[@]}"; do
    wait "$pid" || child_failed=1
done

# Background functions cannot update parent counters. Count persisted statuses.
success_count=$(find "$real_output_root" -mindepth 2 -maxdepth 2 -name .status -exec head -n1 {} \; | awk '$0=="SUCCESS" {n++} END {print n+0}')
fail_count=$(awk -F '\t' 'NR>1 && $3!="SUCCESS" {n++} END {print n+0}' "$GLOBAL_SCORE_FILE")
echo "Completed: success=$success_count failed=$fail_count; summary=$GLOBAL_SCORE_FILE"
[[ "$fail_count" -eq 0 && "$child_failed" -eq 0 ]]
