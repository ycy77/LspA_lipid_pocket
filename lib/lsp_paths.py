"""Portable defaults. Local absolute paths belong in ignored config/local.env."""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
def env_path(name, default):
    return str(Path(os.environ.get(name) or str(default)).expanduser().resolve())

WORKFLOW = env_path("LSP_WORKFLOW_DIR", ROOT / "results/workflow")
DATA = env_path("LSP_DATA_DIR", ROOT / "results/docking")
ANALYSIS = env_path("LSP_ANALYSIS_DIR", ROOT / "results/analysis")
RF_ROOT = env_path("RFDIFFUSION_DIR", ROOT / "external/RFdiffusion")
