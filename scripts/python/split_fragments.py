#!/usr/bin/env python3
import click, os, sys
from pathlib import Path

def read_fasta(p):
    hdr=None; buf=[]; seqs=[]
    with open(p) as fh:
        for ln in fh:
            ln=ln.rstrip()
            if not ln: continue
            if ln.startswith(">"):
                if hdr is not None: seqs.append((hdr,"".join(buf)))
                hdr=ln[1:].split()[0]; buf=[]
            else: buf.append(ln)
    if hdr is not None: seqs.append((hdr,"".join(buf)))
    return seqs

def parse_breaks(tsv_path):
    # Accepte soit bloc #DETAIL (col 'position'), soit liste "pos1;pos2;..." en #SUMMARY
    bps=set()
    with open(tsv_path) as fh:
        for ln in fh:
            if ln.startswith("#DETAIL"):
                parts=ln.strip().split("\t")
                if len(parts)>=4:
                    try: bps.add(int(parts[3]))
                    except: pass
            elif ln.startswith("#SUMMARY") and "breaks" in ln:
                parts=ln.strip().split("\t")
                # header: #SUMMARY alignment gard_json n_breakpoints breaks best_aicc
                if len(parts)>=5 and parts[4]:
                    for t in parts[4].split(";"):
                        t=t.strip()
                        if t:
                            try: bps.add(int(t))
                            except: pass
    return sorted(x for x in bps if x>0)

@click.command()
@click.option("--aln", "aln_path", required=True, help="MSA nucléotidique masqué (FASTA)")
@click.option("--breaks", "breaks_tsv", required=True, help="TSV des cassures (gard_to_breaks / consensus)")
@click.option("--outdir", required=True, help="Répertoire racine des fragments (ex: results/fragments)")
@click.option("--label", default="dnaA", show_default=True, help="Sous-dossier (dnaA/dnaB)")
def main(aln_path, breaks_tsv, outdir, label):
    seqs = read_fasta(aln_path)
    if not seqs:
        print("[ERROR] alignment vide", file=sys.stderr); sys.exit(2)
    L = len(seqs[0][1])
    if any(len(s)!=L for _,s in seqs):
        print("[ERROR] alignment non colinéaire", file=sys.stderr); sys.exit(2)

    bps = parse_breaks(breaks_tsv)
    # bornes internes 1..L-1
    bps = [bp for bp in bps if 1 <= bp < L]

    # Définir segments [start..end] inclusifs
    segments=[]
    prev=1
    for bp in bps + [L]:
        segments.append((prev,bp))
        prev=bp+1

    outroot = Path(outdir) / label
    outroot.mkdir(parents=True, exist_ok=True)

    for i,(a,b) in enumerate(segments, start=1):
        frag = outroot / f"fragment_{i:03d}.fasta"
        with open(frag,"w") as o:
            for h,seq in seqs:
                sub = seq[a-1:b]
                o.write(f">{h}|{a}-{b}\n{sub}\n")

    print(f"[OK] {len(segments)} fragments écrits → {outroot}")

if __name__ == "__main__":
    main()
