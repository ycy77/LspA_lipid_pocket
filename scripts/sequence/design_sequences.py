#!/usr/bin/env python3
"""LigandMPNN membrane-aware chain-B sequence design with side-chain packing."""
import argparse
import os
from pathlib import Path
import shlex
import subprocess
import sys

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-dir', type=Path, required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--mpnn-dir', type=Path, default=os.getenv('LIGANDMPNN_DIR'))
    p.add_argument('--sequences', type=int, default=10)
    p.add_argument('--chains', default='B')
    p.add_argument('--dry-run', action='store_true')
    a=p.parse_args()
    if a.mpnn_dir is None: p.error('provide --mpnn-dir or LIGANDMPNN_DIR')
    if a.sequences < 1: p.error('--sequences must be positive')
    mpnn=a.mpnn_dir.resolve()
    model=mpnn/'model_params/global_label_membrane_mpnn_v_48_020.pt'
    sc=mpnn/'model_params/ligandmpnn_sc_v_32_002_16.pt'
    pdbs=sorted(a.input_dir.glob('*.pdb'))
    if not pdbs: p.error(f'No PDB files in {a.input_dir}')
    if not a.dry_run:
        for path in (mpnn/'run.py',model,sc):
            if not path.is_file(): p.error(f'Missing tool/model: {path}')
    failures=0
    for pdb in pdbs:
        out=a.output_dir.resolve()/pdb.stem
        if out.exists() and any(out.iterdir()):
            print(f'SKIP existing sample: {out} (use a new output directory to rerun)')
            continue
        cmd=[sys.executable,str(mpnn/'run.py'),'--model_type','global_label_membrane_mpnn',
            '--checkpoint_global_label_membrane_mpnn',str(model),'--checkpoint_path_sc',str(sc),
            '--pdb_path',str(pdb.resolve()),'--out_folder',str(out),'--chains_to_design',a.chains,
            '--global_transmembrane_label','1','--batch_size',str(a.sequences),'--number_of_batches','1',
            '--ligand_mpnn_use_side_chain_context','1','--pack_side_chains','1',
            '--number_of_packs_per_design','1','--pack_with_ligand_context','1','--verbose','1']
        print(shlex.join(cmd),flush=True)
        if not a.dry_run:
            out.mkdir(parents=True,exist_ok=True)
            with (out/'run.log').open('w') as log:
                result=subprocess.run(cmd,cwd=mpnn,stdout=log,stderr=subprocess.STDOUT)
            failures+=int(result.returncode!=0)
    return 1 if failures else 0

if __name__=='__main__': raise SystemExit(main())
