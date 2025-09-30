#!/usr/bin/env python3
"""
Mask alignment columns with high gap fraction and trim terminal gap-heavy regions.
- Reads FASTA alignment
- Drops columns with gap fraction > threshold
- Writes a new FASTA alignment built from fresh SeqRecords (no in-place mutation)
"""
import sys
import click
import numpy as np
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

@click.command()
@click.option('--in', 'inf', required=True, help="Input alignment (FASTA)")
@click.option('--out', 'outf', required=True, help="Output masked alignment (FASTA)")
@click.option('--mask-gap-frac', type=float, default=0.5, show_default=True)
def main(inf, outf, mask_gap_frac):
    try:
        aln = AlignIO.read(inf, 'fasta')  # MultipleSeqAlignment
    except Exception as e:
        print(f"[FATAL] Cannot read alignment '{inf}': {e}", file=sys.stderr)
        sys.exit(1)

    if len(aln) == 0:
        print("[ERROR] Empty alignment.", file=sys.stderr)
        sys.exit(2)

    # Convert to array
    arr = np.array([list(str(rec.seq)) for rec in aln], dtype='<U1')
    L = arr.shape[1]
    gaps = (arr == '-') | (arr == 'N')  # optionally also mask all-N columns if desired
    gap_frac = gaps.mean(axis=0)

    keep = gap_frac <= mask_gap_frac
    kept_cols = int(keep.sum())

    if kept_cols == 0:
        # Write an empty alignment with same IDs (or fall back to original)
        # Safer: write the original but warn
        AlignIO.write(aln, outf, 'fasta')
        print(f"[WARN] All columns masked (threshold={mask_gap_frac}). Wrote original alignment unchanged → {outf}")
        sys.exit(0)

    arr2 = arr[:, keep]
    # Build new SeqRecords cleanly
    new_records = []
    for i, rec in enumerate(aln):
        new_seq = ''.join(arr2[i, :].tolist())
        new_records.append(SeqRecord(Seq(new_seq), id=rec.id, description=rec.description))

    new_aln = MultipleSeqAlignment(new_records)
    AlignIO.write(new_aln, outf, 'fasta')
    kept_pct = 100.0 * kept_cols / L
    print(f"[INFO] Kept {kept_cols}/{L} columns ({kept_pct:.1f}%) → {outf}")

if __name__ == '__main__':
    main()
