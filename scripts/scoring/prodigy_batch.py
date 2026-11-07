#!/usr/bin/env python3
"""Score PDB complexes with PRODIGY and optionally select one per backbone."""
import argparse
import csv
import math
from pathlib import Path
import re
import shutil
import subprocess

def parse_scores(text):
    values=[]
    for label in ['Predicted binding affinity','Predicted dissociation constant']:
        line=next((line for line in text.splitlines() if label in line),None)
        if line is None: raise ValueError(f'Missing {label}')
        token=line.rsplit(':',1)[-1].strip().split()[0]
        value=float(token)
        if not math.isfinite(value): raise ValueError(f'Non-finite {label}')
        values.append(value)
    return values

def write_tsv(path, header, rows):
    with path.open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(header);w.writerows(rows)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-dir',required=True,type=Path)
    p.add_argument('--output-dir',required=True,type=Path)
    p.add_argument('--pattern',default='*.pdb',help='Use */minimized/*.pdb or */packed/*.pdb for grouped selection')
    p.add_argument('--grouped',action='store_true',help='Group by first directory below input-dir')
    p.add_argument('--summary',type=Path,help='HADDOCK global_score_summary.tsv to retain scores/statuses')
    p.add_argument('--prodigy',default='prodigy')
    p.add_argument('--chains',nargs=2,default=['A','B'])
    a=p.parse_args()
    pdbs=sorted(a.input_dir.glob(a.pattern))
    if not pdbs: p.error(f'No PDBs match {a.pattern}')
    a.output_dir.mkdir(parents=True,exist_ok=True)
    logs=a.output_dir/'logs'; logs.mkdir(exist_ok=True)
    scores={}
    if a.summary:
        # Accept original Chinese header and the normalized English header.
        with a.summary.open() as f:
            for row in csv.reader(f,delimiter='\t'):
                if len(row)>=3 and row[2]=='SUCCESS': scores[row[0]]=row[1]
    rows=[]; failures=[]; best={}
    for index,pdb in enumerate(pdbs):
        relative=pdb.relative_to(a.input_dir)
        if a.grouped and len(relative.parts)<2: p.error('--grouped requires one subdirectory per backbone')
        group=relative.parts[0] if a.grouped else pdb.stem.removesuffix('_emref_1')
        if a.summary and group not in scores: continue
        result=subprocess.run([a.prodigy,str(pdb.resolve()),'--selection',*a.chains],text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        (logs/f'{index:06d}_{pdb.stem}.log').write_text(result.stdout)
        try:
            if result.returncode: raise ValueError(f'exit={result.returncode}')
            dg,kd=parse_scores(result.stdout)
        except (ValueError,IndexError) as exc:
            failures.append([str(relative),str(exc)]);continue
        if a.grouped:
            rows.append([group,pdb.name,dg,kd])
            if group not in best or dg<best[group][1]: best[group]=(pdb,dg,kd)
        else:
            rows.append([group,scores.get(group,''),'SUCCESS',dg,kd])
    if a.grouped:
        write_tsv(a.output_dir/'final_all.tsv',['Backbone','Model','ΔG(kcal/mol)','Kd'],rows)
        best_dir=a.output_dir/'best_structures';best_dir.mkdir(exist_ok=True)
        summaries=[]
        for group,(pdb,dg,kd) in sorted(best.items()):
            shutil.copy2(pdb,best_dir/f'{group}_best.pdb')
            # Counts are from unique ATOM residues in the peptide chain.
            residues={line[22:27] for line in pdb.read_text().splitlines() if line.startswith('ATOM  ') and line[21:22]==a.chains[1]}
            summaries.append([group,pdb.name,dg,kd,len(residues)])
        write_tsv(a.output_dir/'final_summary.tsv',['Backbone','Best_Model','ΔG(kcal/mol)','Kd','B_chain_AA_count'],summaries)
        write_tsv(a.output_dir/'final_summary_with_AA.tsv',['Backbone','Best_Model','ΔG(kcal/mol)','Kd','B_chain_AA_count'],summaries)
    else:
        write_tsv(a.output_dir/'global_score_with_prodigy.tsv',['Model','HADDOCK_Score','Status','ΔG(kcal/mol)','Kd(M)'],rows)
    write_tsv(a.output_dir/'failures.tsv',['PDB','Error'],failures)
    print(f'Scored {len(rows)} structures; failed={len(failures)}')
    return 1 if failures or not rows else 0

if __name__=='__main__': raise SystemExit(main())
