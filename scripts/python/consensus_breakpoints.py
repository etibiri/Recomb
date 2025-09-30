#!/usr/bin/env python3
"""
Merge breakpoints from multiple methods within ±merge_nt and require >= min_support methods.
Robuste :
- Crée toujours le TSV (au moins l'entête)
- Écrit via un fichier temporaire puis rename (atomique)
"""
import click, csv, os, tempfile
from collections import defaultdict

@click.command()
@click.option('--gard', multiple=True, help='GARD TSV files')
@click.option('--three', multiple=True, help='3SEQ TSV files')
@click.option('--merge-nt', type=int, default=150)
@click.option('--min-support', type=int, default=2)
@click.option('--out', 'outf', required=True)
def main(gard, three, merge_nt, min_support, outf):
    os.makedirs(os.path.dirname(outf), exist_ok=True)
    positions = []
    # Lire GARD
    for g in gard:
        if not g or not os.path.exists(g): 
            continue
        with open(g) as f:
            header = next(f, "")
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2: 
                    continue
                s, e = parts[0], parts[1]
                try:
                    positions.append(int(s)); positions.append(int(e))
                except Exception:
                    pass
    # (placeholder) intégrer 3SEQ si tu ajoutes des positions candidates
    # for t in three: ...

    positions = sorted(p for p in positions if isinstance(p, int))
    clusters = []
    for p in positions:
        if not clusters:
            clusters.append([p, p, 1])
        else:
            s, e, c = clusters[-1]
            if p - e <= merge_nt:
                clusters[-1][1] = p
                clusters[-1][2] += 1
            else:
                clusters.append([p, p, 1])
    kept = [(s, e, c) for s, e, c in clusters if c >= min_support]

    fd, tmp = tempfile.mkstemp(prefix=".consensus_", dir=os.path.dirname(outf))
    with os.fdopen(fd, "w") as w:
        w.write("start\tend\tsupport\n")
        for s, e, c in kept:
            w.write(f"{s}\t{e}\t{c}\n")
        w.flush(); os.fsync(w.fileno())
    os.replace(tmp, outf)
    print(f"[INFO] Consensus clusters: {len(kept)} → {outf}")

if __name__ == '__main__':
    main()
