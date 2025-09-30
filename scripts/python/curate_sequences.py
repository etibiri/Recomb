#!/usr/bin/env python3
"""Filter sequences by length, N content, and deduplicate nearly-identical records.
Outputs a cleaned FASTA preserving record IDs.
"""
import sys, hashlib
from Bio import SeqIO
import click

@click.command()
@click.option('--in', 'inf', required=True, type=click.Path(exists=True))
@click.option('--out', 'outf', required=True)
@click.option('--min-len', type=int, required=True)
@click.option('--max-N-frac', type=float, default=0.02)
def main(inf, outf, min_len, max_n_frac):
    seen = set()
    out = []
    for rec in SeqIO.parse(inf, 'fasta'):
        seq = str(rec.seq).upper()
        if len(seq) < min_len:
            continue
        if (seq.count('N')/len(seq)) > max_n_frac:
            continue
        # hash to deduplicate
        h = hashlib.md5(seq.encode()).hexdigest()
        if h in seen:
            continue
        seen.add(h)
        out.append(rec)
    SeqIO.write(out, outf, 'fasta')
    print(f"[INFO] Wrote {len(out)} records → {outf}")

if __name__ == '__main__':
    main()

