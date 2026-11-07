#!/usr/bin/env python3
"""RF backbone generation for cyclic design, local refinement and linear controls."""
import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode', choices=['cyclic', 'refine', 'linear'], default='cyclic')
    p.add_argument('--presets', type=Path, default=ROOT/'config/design_presets.json')
    p.add_argument('--target', type=Path)
    p.add_argument('--output-dir', type=Path)
    p.add_argument('--rf-dir', type=Path, default=Path(os.getenv('RFDIFFUSION_DIR', ROOT/'external/RFdiffusion')))
    p.add_argument('--checkpoint', type=Path)
    p.add_argument('--adj-file', type=Path, default=os.getenv('LSP_ADJ_FILE'))
    p.add_argument('--scaffold-dir', type=Path, default=os.getenv('LSP_SCAFFOLD_DIR'))
    p.add_argument('--t-values', type=int, nargs='+')
    p.add_argument('--count-per-t', type=int)
    p.add_argument('--noise', type=float, help='Override both CA and frame noise scales')
    p.add_argument('--contigs')
    p.add_argument('--hotspots', default='A125,A129,A130,A132,A136,A137,A140,A141,A144')
    p.add_argument('--dry-run', action='store_true', help='Print commands; do not launch RF')
    a = p.parse_args()
    cfg = json.loads(a.presets.read_text())[a.mode]
    if a.mode == 'refine':
        target = a.target or (Path(os.environ['LSP_REFINE_TEMPLATE']) if os.getenv('LSP_REFINE_TEMPLATE') else None)
        if target is None:
            p.error('refine requires --target or LSP_REFINE_TEMPLATE (selected A+B complex)')
    else:
        target = a.target or (Path(os.environ['LSP_TARGET_PDB']) if os.getenv('LSP_TARGET_PDB') else None)
        if target is None:
            p.error('provide --target or LSP_TARGET_PDB (LspA chain A PDB)')
    if not a.adj_file or not a.scaffold_dir:
        p.error('provide --adj-file/--scaffold-dir or LSP_ADJ_FILE/LSP_SCAFFOLD_DIR')
    rf = a.rf_dir.resolve()
    checkpoint = (a.checkpoint or rf/'models/Complex_base_ckpt.pt').resolve()
    out = (a.output_dir or Path(os.getenv('LSP_WORKFLOW_DIR', ROOT/'results/workflow'))/'01_rf_diffusion'/cfg['prefix']).resolve()
    count = cfg['count_per_t'] if a.count_per_t is None else a.count_per_t
    ts = cfg['t_values'] if a.t_values is None else a.t_values
    if count < 1 or any(t < 1 for t in ts):
        p.error('counts and T values must be positive')
    if not a.dry_run:
        for path in [rf/'scripts/run_inference.py', target, checkpoint, a.adj_file, a.scaffold_dir]:
            if not path.exists(): p.error(f'Missing input: {path}')
        out.mkdir(parents=True, exist_ok=True)
    failures = 0
    for t in ts:
        for i in range(1, count+1):
            prefix = out/f'{cfg["prefix"]}_T{t}_{i:03d}'
            # RF appends the design index (_0) to a supplied output prefix.
            done = Path(str(prefix)+'_0.pdb')
            if done.is_file() and done.stat().st_size:
                print(f'SKIP {done}')
                continue
            noise = cfg['noise'] if a.noise is None else a.noise
            cmd = [sys.executable, str(rf/'scripts/run_inference.py'),
                f'inference.input_pdb={target.resolve()}',
                f'contigmap.contigs=[{a.contigs or cfg["contigs"]}]',
                f'inference.output_prefix={prefix}', f'inference.ckpt_override_path={checkpoint}',
                f'diffuser.T={t}', 'inference.num_designs=1', 'inference.design_startnum=0',
                'inference.deterministic=False', f'ppi.hotspot_res=[{a.hotspots}]',
                'potentials.guiding_potentials=["type:monomer_contacts,weight:0.05", "type:interface_ncontacts,weight:0.3"]',
                f'scaffoldguided.target_adj={a.adj_file.resolve()}',
                f'scaffoldguided.scaffold_dir={a.scaffold_dir.resolve()}',
                f'denoiser.noise_scale_ca={noise}', f'denoiser.noise_scale_frame={noise}']
            if cfg['cyclic']: cmd += ['inference.cyclic=True']
            if 'length' in cfg: cmd += [f'contigmap.length={cfg["length"]}']
            print(shlex.join(cmd), flush=True)
            if not a.dry_run:
                with Path(str(prefix)+'.log').open('w') as log:
                    result = subprocess.run(cmd, cwd=rf, stdout=log, stderr=subprocess.STDOUT)
                failures += int(result.returncode != 0 or not done.exists())
    return 1 if failures else 0

if __name__ == '__main__':
    raise SystemExit(main())
