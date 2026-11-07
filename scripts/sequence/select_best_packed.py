#!/usr/bin/env python3
"""Select the highest LigandMPNN overall_confidence for each backbone."""
import argparse
import csv
from pathlib import Path
import re
import shutil

def fasta_records(path):
    header=None
    seq=[]
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            if header is not None: yield header,''.join(seq)
            header=line[1:]; seq=[]
        elif line.strip(): seq.append(line.strip())
    if header is not None: yield header,''.join(seq)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-dir',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--receptor-fasta',type=Path)
    a=p.parse_args()
    fastas=sorted(a.input_dir.glob('*/seqs/*.fa'))
    if not fastas: p.error('No */seqs/*.fa files found')
    a.output_dir.mkdir(parents=True,exist_ok=True)
    rows=[]; failures=[]
    for fasta in fastas:
        sample=fasta.parent.parent
        candidates=[]
        for header,seq in fasta_records(fasta):
            ident=re.search(r'(?:^|[,\s])id=(\d+)(?:,|\s|$)',header)
            conf=re.search(r'overall_confidence=([\d.eE+-]+)',header)
            if ident and conf: candidates.append((float(conf[1]),int(ident[1]),header,seq))
        if not candidates:
            failures.append(str(fasta)); continue
        score,ident,header,seq=max(candidates,key=lambda v:(v[0],v[1]))
        packed=sorted((sample/'packed').glob(f'{sample.name}_packed_{ident}_*.pdb'))
        if not packed:
            failures.append(f'No packed PDB for {sample.name} id={ident}'); continue
        dest=a.output_dir/sample.name
        (dest/'seqs').mkdir(parents=True,exist_ok=True)
        (dest/'packed').mkdir(exist_ok=True)
        receptor=a.receptor_fasta.read_text().rstrip()+'\n' if a.receptor_fasta else ''
        (dest/'seqs'/f'{sample.name}_best.fa').write_text(receptor+f'>{header}\n{seq}\n')
        for pdb in packed: shutil.copy2(pdb,dest/'packed'/pdb.name)
        rows.append([sample.name,ident,score,len(packed)])
    with (a.output_dir/'selection.tsv').open('w') as f:
        writer=csv.writer(f,delimiter='\t'); writer.writerow(['Backbone','id','overall_confidence','packed_count']); writer.writerows(rows)
    for error in failures: print('ERROR:',error)
    print(f'Selected {len(rows)} backbones; failed={len(failures)}')
    return bool(failures)

if __name__=='__main__': raise SystemExit(main())
