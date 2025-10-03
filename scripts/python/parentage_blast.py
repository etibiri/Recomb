#!/usr/bin/env python3
import os, glob, subprocess, tempfile, shutil, random, argparse, sys

def read_fasta(path):
    seqs=[]; hdr=None; buf=[]
    with open(path) as fh:
        for line in fh:
            line=line.rstrip()
            if not line: 
                continue
            if line.startswith(">"):
                if hdr is not None:
                    seqs.append((hdr,"".join(buf)))
                hdr=line[1:].split()[0]; buf=[]
            else:
                buf.append(line)
        if hdr is not None:
            seqs.append((hdr,"".join(buf)))
    return seqs

def read_breaks(consensus_path):
    bps=set()
    if not os.path.isfile(consensus_path):
        return []
    with open(consensus_path) as fh:
        _ = fh.readline()  # header
        for line in fh:
            for tok in line.replace(",","\t").split():
                try:
                    v=int(tok)
                    if v>0: bps.add(v)
                except:
                    pass
    return sorted(bps)

def slice_columns(msa, start, end):
    s = start-1; e = end
    return [(sid, seq[s:e]) for sid,seq in msa]

def consensus_of_slice(sliced):
    if not sliced:
        return ""
    L=len(sliced[0][1]); out=[]
    for i in range(L):
        freq={}
        for sid,s in sliced:
            c=s[i]
            freq[c]=freq.get(c,0)+1
        out.append(sorted(freq.items(), key=lambda kv:(-kv[1],kv[0]))[0][0])
    return "".join(out)

def build_blast_db_if_any(refs_dir, tmpdir):
    pats = ("*.fa","*.fasta","*.fna")
    fa_list=[]
    for p in pats:
        fa_list.extend(sorted(glob.glob(os.path.join(refs_dir,p))))
    if not fa_list:
        return None
    cat=os.path.join(tmpdir,"refs.cat.fasta")
    with open(cat,"w") as out:
        for fp in fa_list:
            with open(fp) as fh:
                out.write(fh.read())
    db=os.path.join(tmpdir,"refdb")
    subprocess.run(["makeblastdb","-in",cat,"-dbtype","nucl","-out",db],
                   check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return db

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--msa", required=True)
    ap.add_argument("--consensus", required=True)
    ap.add_argument("--refs", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--threads", type=int, default=1)
    args = ap.parse_args()

    random.seed(args.seed)

    msa=read_fasta(args.msa)
    if not msa:
        print("Masked alignment is empty.", file=sys.stderr)
        sys.exit(2)
    L=len(msa[0][1])

    bps=read_breaks(args.consensus)
    cut_points=[bp for bp in bps if 1 <= bp < L]

    segments=[]
    prev=1
    for bp in (cut_points + [L]):
        segments.append((prev,bp))
        prev=bp+1

    tmpdir=tempfile.mkdtemp(prefix="blast_parentage_")
    db=None
    try:
        have_refs=os.path.isdir(args.refs) and any(
            glob.glob(os.path.join(args.refs, p)) for p in ("*.fa","*.fasta","*.fna")
        )
        if have_refs:
            try:
                db=build_blast_db_if_any(args.refs,tmpdir)
            except Exception as e:
                print("[WARN] build DB:", e, file=sys.stderr)

        rows=[]
        for idx,(a,b) in enumerate(segments, start=1):
            sliced=slice_columns(msa,a,b)
            cons=consensus_of_slice(sliced)
            qfa=os.path.join(tmpdir, "seg_%03d.fa" % idx)
            with open(qfa,"w") as o:
                o.write(">seg_%03d_%d_%d\n" % (idx,a,b))
                o.write(cons+"\n")

            major="NA"; minor="NA"; score="NA"
            if db and cons:
                cmd=["blastn","-query",qfa,"-db",db,
                     "-outfmt","6 qseqid sseqid pident length evalue bitscore qcovs",
                     "-max_target_seqs","5",
                     "-num_threads", str(args.threads)]
                try:
                    cp=subprocess.run(cmd,check=True,text=True,capture_output=True)
                    lines=[ln for ln in cp.stdout.splitlines() if ln.strip()]
                    if lines:
                        best=lines[0].split("\t")
                        if len(best)>=7:
                            major=best[1]; score=best[6]
                        if len(lines)>1:
                            second=lines[1].split("\t")
                            if len(second)>=2:
                                minor=second[1]
                except Exception as e:
                    print("[WARN] BLAST seg %d: %s" % (idx, str(e)), file=sys.stderr)

            rows.append(("%d-%d" % (a,b), major, minor, score))

        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        with open(args.out,"w") as out:
            out.write("breakpoint_interval\tmajor_parent\tminor_parent\tscore\n")
            for r in rows:
                out.write("\t".join(r)+"\n")
    finally:
        try:
            shutil.rmtree(tmpdir, ignore_errors=True)
        except:
            pass

if __name__ == "__main__":
    main()
