#!/usr/bin/env python3
"""Extract one PDB chain; optional renumbering prepares HADDOCK chain B as 1..N."""
import argparse
from pathlib import Path

def extract(source, chain, renumber):
    result=[]; residues={}
    for line in source.read_text().splitlines():
        if line[:6] not in ('ATOM  ','HETATM') or len(line)<27 or line[21]!=chain: continue
        if renumber:
            key=line[22:27]
            residues.setdefault(key,len(residues)+1)
            line=line[:22]+f'{residues[key]:4d} '+line[27:]
        result.append(line)
    if not result: raise ValueError(f'No chain {chain} in {source}')
    return '\n'.join(result+['TER','END'])+'\n'

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--chain',default='B')
    p.add_argument('--renumber',action='store_true')
    p.add_argument('--pattern',default='*.pdb',help='Glob relative to a directory; e.g. */packed/*.pdb')
    a=p.parse_args()
    if len(a.chain)!=1: p.error('--chain must be one character')
    files=[a.input] if a.input.is_file() else sorted(a.input.glob(a.pattern))
    if not files: p.error('No input PDBs')
    names=[f.stem for f in files]
    if len(names)!=len(set(names)): p.error('Duplicate basenames would overwrite files; narrow --pattern')
    a.output_dir.mkdir(parents=True,exist_ok=True)
    for src in files:
        dest=a.output_dir/f'{src.stem}_{a.chain}_chain.pdb'
        dest.write_text(extract(src,a.chain,a.renumber))
    print(f'Extracted {len(files)} chains to {a.output_dir}')

if __name__=='__main__': main()
