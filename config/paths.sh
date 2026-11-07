#!/usr/bin/env bash
# Sourced by run.sh or a shell entry point. Do not store credentials here.
LSP_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
if [[ -f "${LSP_CONFIG:-$LSP_ROOT/config/local.env}" ]]; then
    set -a
    source "${LSP_CONFIG:-$LSP_ROOT/config/local.env}"
    set +a
fi
export LSP_ROOT
export LSP_WORKFLOW_DIR="${LSP_WORKFLOW_DIR:-$LSP_ROOT/results/workflow}"
export LSP_DATA_DIR="${LSP_DATA_DIR:-$LSP_ROOT/results/docking}"
export LSP_ANALYSIS_DIR="${LSP_ANALYSIS_DIR:-$LSP_ROOT/results/analysis}"
export RFDIFFUSION_DIR="${RFDIFFUSION_DIR:-$LSP_ROOT/external/RFdiffusion}"
export LIGANDMPNN_DIR="${LIGANDMPNN_DIR:-$RFDIFFUSION_DIR/LigandMPNN}"
export MPLBACKEND=Agg
