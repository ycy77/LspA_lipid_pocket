# LspA peptide inhibitor design

LspA cyclic peptide inhibitor design workflow, including RF backbone generation, LigandMPNN sequence design, HADDOCK3 docking, PRODIGY scoring, structural refinement, and analysis plots.

## Design parameters

| Stage | Cyclic peptide length | CA / frame noise | Diffusion steps T | Backbones |
|---|---|---|---|---|
| Initial design | 8–15 aa | 0.1 / 0.1 | 50–200 in steps of 10 | 100 per T; 1,600 total |
| Second-round design | 8–12 aa | 0.2 / 0.2 | 20, 30, 40, 50 | 100 per T; 400 total |

LigandMPNN generates 10 sequences per backbone. Presets are defined in `config/design_presets.json` and can be overridden on the command line.

## Directory structure

```text
LspA_lipid_pocket/
├── README.md
├── run.sh                       # Common entry point
├── requirements-analysis.txt    # Analysis and plotting dependencies
├── config/                      # Parameters and path configuration
├── lib/                         # Shared path utilities
└── scripts/
    ├── design/                  # RF backbone generation
    ├── preparation/             # Receptor replacement, chain B extraction and merging
    ├── sequence/                # LigandMPNN and sequence selection
    ├── docking/                 # HADDOCK3 docking
    ├── scoring/                 # PRODIGY scoring
    ├── refinement/              # OpenMM minimization and Rosetta interface metrics
    └── analysis/                # Composition, logos, RMSD, pockets and score distributions
```

## Environment and inputs

Use Python >=3.9 and Bash. Activate the appropriate tool environment before each stage:

| Stage | Dependencies |
|---|---|
| Backbone generation | RFdiffusion/RFpeptide with `inference.cyclic` support, model weights, target adjacency and a scaffold library |
| Sequence design | LigandMPNN with membrane-model and side-chain-model weights |
| Structure preparation | PyMOL; standalone chain extraction requires only the Python standard library |
| Docking and scoring | HADDOCK3/CNS and command-line PRODIGY |
| Structural refinement | OpenMM/PDBFixer and PyRosetta |
| Analysis and plotting | `pip install -r requirements-analysis.txt` |

Store structures, model weights and results outside the code directory and specify their paths in the configuration.

```bash
cd LspA_lipid_pocket
# Create a configuration on a new machine; edit local.env directly if it already exists.
cp config/example.env config/local.env
```

Edit `config/local.env` and set at least:

```bash
export RFDIFFUSION_DIR=/path/to/RFdiffusion
export LIGANDMPNN_DIR=/path/to/LigandMPNN
export LSP_TARGET_PDB=/path/to/LspA_chain_A.pdb
export LSP_ADJ_FILE=/path/to/LspA_adj.pt
export LSP_SCAFFOLD_DIR=/path/to/ppi_scaffolds
```

`LSP_TARGET_PDB` should contain receptor chain A with its original residue numbering. Optionally set `LSP_REFERENCE_PDB` as the reference for receptor replacement. For the second round, use `--target` to specify a template complex containing chains A and B. Set paths in `local.env` for your installation.

Default outputs are `results/workflow`, `results/docking` and `results/analysis`; override them with `LSP_WORKFLOW_DIR`, `LSP_DATA_DIR` and `LSP_ANALYSIS_DIR`. Run the examples below from the repository root using the default output paths. `local.env` and generated results are excluded by `.gitignore`.

## Initial design

```bash
# Show help and preview a generation command.
bash run.sh --help
bash run.sh design/generate_backbones.py --t-values 50 --count-per-t 1 --dry-run

# Generate 1,600 backbones of 8–15 aa with noise 0.1.
bash run.sh design/generate_backbones.py --mode cyclic

# Replace receptor chain A.
bash run.sh preparation/replace_receptor.py

# Design 10 sequences per backbone and pack side chains.
bash run.sh sequence/design_sequences.py \
  --input-dir results/workflow/02_rf_diffusion_gai/cyclic \
  --output-dir results/workflow/03_MPNN/cyclic

# Select the sequence with the highest overall_confidence for each backbone.
bash run.sh sequence/select_best_packed.py \
  --input-dir results/workflow/03_MPNN/cyclic \
  --output-dir results/workflow/04_select_MPNN/cyclic

# Extract chain B and renumber residues consecutively.
bash run.sh preparation/extract_chain.py \
  --input results/workflow/04_select_MPNN/cyclic --pattern '*/packed/*.pdb' \
  --output-dir results/workflow/05_peptide_chains/cyclic --chain B --renumber

# Run docking.
bash run.sh docking/run_haddock3_batch.sh

# Score docked structures.
bash run.sh scoring/prodigy_batch.py \
  --input-dir results/docking/xiao_haddock3_cyclic_batch/best_pdb_all \
  --summary results/docking/xiao_haddock3_cyclic_batch/global_score_summary.tsv \
  --output-dir results/docking/xiao_haddock3_cyclic_batch
```

## Second-round design

Use a selected first-round complex as the template to generate 8–12 aa cyclic peptides with noise 0.2.

```bash
bash run.sh design/generate_backbones.py --mode refine \
  --target /path/to/selected_complex.pdb

LSP_INPUT_DIR="$PWD/results/workflow/01_rf_diffusion/cyclic_refine" \
LSP_OUTPUT_DIR="$PWD/results/workflow/02_rf_diffusion_gai/cyclic_refine" \
  bash run.sh preparation/replace_receptor.py

bash run.sh sequence/design_sequences.py \
  --input-dir results/workflow/02_rf_diffusion_gai/cyclic_refine \
  --output-dir results/workflow/03_MPNN/cyclic_refine

bash run.sh sequence/select_best_packed.py \
  --input-dir results/workflow/03_MPNN/cyclic_refine \
  --output-dir results/workflow/04_select_MPNN/cyclic_refine

bash run.sh preparation/extract_chain.py \
  --input results/workflow/04_select_MPNN/cyclic_refine --pattern '*/packed/*.pdb' \
  --output-dir results/workflow/05_peptide_chains/cyclic_refine --chain B --renumber

LSP_INPUT_DIR="$PWD/results/workflow/05_peptide_chains/cyclic_refine" \
LSP_OUTPUT_DIR="$PWD/results/docking/cyclic_refine" \
  bash run.sh docking/run_haddock3_batch.sh

bash run.sh scoring/prodigy_batch.py \
  --input-dir results/docking/cyclic_refine/best_pdb_all \
  --summary results/docking/cyclic_refine/global_score_summary.tsv \
  --output-dir results/docking/cyclic_refine
```

## Optional refinement and grouped scoring

Merge all second-round packed chain B structures into the reference complex, minimize, and select by backbone:

```bash
LSP_OUTPUT_DIR="$PWD/results/workflow/merged/cyclic_refine" \
  bash run.sh preparation/merge_packed_peptides.py

bash run.sh refinement/minimize_cyclic_openmm.py \
  --input-dir results/workflow/merged/cyclic_refine \
  --output-dir results/docking/refined_openmm

bash run.sh scoring/prodigy_batch.py \
  --input-dir results/docking/refined_openmm \
  --pattern '*/merge_pdb/minimized/*.pdb' --grouped \
  --output-dir results/docking/refined_scores

LSP_SUMMARY_TSV="$PWD/results/docking/refined_scores/final_summary.tsv" \
LSP_PDB_ROOT="$PWD/results/docking/refined_openmm" \
LSP_MODEL_SUBDIR="merge_pdb/minimized" \
  bash run.sh refinement/rosetta_interface_metrics.py

bash run.sh analysis/plot_rosetta_distributions.py
```

OpenMM preserves relative input paths and appends `minimized/`; place the output directory outside the input directory. Grouped scoring produces `final_all.tsv`, `final_summary.tsv` and `best_structures/`.

## Analysis scripts

All analysis outputs default to `results/analysis`.

| Script | Main input variables |
|---|---|
| `analysis/plot_peptide_composition.py` | `LSP_INPUT_DIR`: PDB directory; use a new `LSP_ANALYSIS_DIR` when changing inputs |
| `analysis/plot_sequence_logos.py` | `LSP_INPUT_DIR`: PDB directory |
| `analysis/calculate_peptide_rmsd.py` | `LSP_INITIAL_DIR`: sample/packed/unique PDB; `LSP_DOCKED_DIR`: directory of `*_best.pdb` files matched by original sample names |
| `analysis/analyze_binding_pockets.py` | `LSP_INPUT_DIR`: docked complexes; `LSP_SCORES_TSV`: score table |
| `analysis/compare_prodigy_by_length.py` | `LSP_INITIAL_TSV`, `LSP_OPTIMIZED_TSV`: paired scores for the same backbones before and after optimization |
| `refinement/rosetta_interface_metrics.py` | `LSP_SUMMARY_TSV`, `LSP_PDB_ROOT`, `LSP_MODEL_SUBDIR`; `LSP_WORKERS` controls worker count |
| `analysis/plot_rosetta_distributions.py` | `LSP_INPUT_CSV`: CSV containing dSASA/ddG/BUH/CYCttE columns |

```bash
LSP_INPUT_DIR=/path/to/complexes \
LSP_ANALYSIS_DIR="$PWD/results/composition" \
  bash run.sh analysis/plot_peptide_composition.py
```

## Execution controls

- Override RF parameters with `--t-values`, `--count-per-t`, `--noise` and `--contigs`; `--dry-run` previews commands. Use `--mode linear` for linear controls, which do not enter the cyclic-peptide docking workflow.
- HADDOCK defaults to 4 cores per task and 1 concurrent task. Adjust with `LSP_CPUS_PER_TASK` and `LSP_MAX_CONCURRENT`.
- `LSP_PREPARE_ONLY=1` generates docking configurations only. Use a separate `LSP_OUTPUT_DIR` and a fresh directory for production runs.
- RF skips existing nonempty `_0.pdb` files; MPNN and HADDOCK skip existing task directories. Choose a new output directory to rerun a calculation.
- Set `LSP_PYTHON` to select the entry-point interpreter; activate the appropriate environment for each stage.
