################################################################################
# Recombination Analysis — CMD DNA-A
# Author: Ezechiel B. TIBIRI
# Reproducible, offline-ready Snakemake workflow.
# Uses existing Conda envs (envs/) and provided scripts (scripts/).
################################################################################

import os

configfile: "config.yaml"

# --- Shortcuts from config -----------------------------------------------------
DATA_FASTA   = config["inputs"]["dnaA_fasta"]
METADATA_TSV = config["inputs"]["metadata"]
REFS_DIR     = config["inputs"]["refs_dir"]

SEEDS        = config["seeds"]
MAFFT_OPTS   = config["align"]
IQTREE_OPTS  = config["tools"]["iqtree_opts"]
HYPHY_OPTS   = config["tools"]["hyphy_opts"]
MASKING      = config["masking"]

RES_DEFAULT  = config["resources"]["default"]
RES_HEAVY    = config["resources"]["heavy"]

RESULTS   = "results"
ALIGN_DIR = f"{RESULTS}/alignments"
RECOMB_DIR= f"{RESULTS}/recombination"
FRAG_DIR  = f"{RESULTS}/fragments"
TREE_DIR  = f"{RESULTS}/phylogeny/iqtree"
NET_DIR   = f"{RESULTS}/network"
PARENT_DIR= f"{RESULTS}/parentage"
SUM_DIR   = f"{RESULTS}/summary"
LOGDIR    = f"{RESULTS}/logs"
BMDIR     = f"{RESULTS}/benchmarks"
#HYPHY_ANALYSES = "external/hyphy-analyses"
#GARD_BF_PATH = f"{HYPHY_ANALYSES}/GARD/GARD.bf" 



# ensure base directories exist even on dry-run
for d in [RESULTS, ALIGN_DIR, RECOMB_DIR, FRAG_DIR, TREE_DIR, NET_DIR, PARENT_DIR, SUM_DIR, LOGDIR, BMDIR]:
    os.makedirs(d, exist_ok=True)

# --- Helpers ------------------------------------------------------------------
def log_of(rule_name):   return f"{LOGDIR}/{rule_name}.log"
def bench_of(rule_name): return f"{BMDIR}/{rule_name}.tsv"

# --- Rule: all ----------------------------------------------------------------
rule all:
    message: "Aggregating final outputs for CMD DNA-A recombination analysis."
    input:
        # Curation
        "results/curation/dnaA.cleaned.fasta",
        # Alignments
        f"{ALIGN_DIR}/dnaA.aln.fasta",
        f"{ALIGN_DIR}/dnaA.masked.fasta",
        # Recombination (GARD → JSON → TSV) + consensus
        f"{RECOMB_DIR}/dnaA.gard.json",
        f"{RECOMB_DIR}/dnaA.gard_breakpoints.tsv",
        f"{RECOMB_DIR}/events_all.tsv",
        f"{RECOMB_DIR}/consensus_breakpoints.tsv",
        # Parentage
        f"{PARENT_DIR}/parent_calls.tsv",
        # Phylogeny
        f"{TREE_DIR}/dnaA.treefile",
        f"{TREE_DIR}/dnaA.iqtree",
        f"{TREE_DIR}/dnaA.log",
        # Network
        f"{NET_DIR}/neighborNet.pdf",
        f"{NET_DIR}/neighborNet.png",
        f"{NET_DIR}/taxa_colors.csv",
        # Summary
        f"{SUM_DIR}/recombination_summary.tsv"
    threads: 1
    resources:
        mem_mb=RES_DEFAULT["mem_mb"],
        runtime=RES_DEFAULT["runtime_min"]
    log: log_of("all")
    benchmark: bench_of("all")
    conda: "envs/base-python.yaml"
    shell:
        r"""
        echo "[OK] Final targets exist." > {log}
        """

# --- Sequence curation & normalization ----------------------------------------
rule curate_sequences:
    input:
        fasta = DATA_FASTA
    output:
        cleaned = "results/curation/dnaA.cleaned.fasta"
    params:
        seed   = SEEDS["cleaning"],
        length = SEEDS["min_len"],
        nqc    = SEEDS["max_N_frac"]
    threads: 2
    resources:
        mem_mb=RES_DEFAULT["mem_mb"],
        runtime=RES_DEFAULT["runtime_min"]
    log: log_of("curate_sequences")
    benchmark: bench_of("curate_sequences")
    conda: "envs/base-python.yaml"
    shell:
        r"""
        mkdir -p results/curation
        python scripts/python/curate_sequences.py \
          --in {input.fasta} \
          --out {output.cleaned} \
          --min-len {params.length} \
          --max-N-frac {params.nqc} \
          > {log} 2>&1
        """

# --- MAFFT alignment -----------------------------------------------------------
rule mafft_align:
    input:  "results/curation/dnaA.cleaned.fasta"
    output: f"{ALIGN_DIR}/dnaA.aln.fasta"
    params: opts = MAFFT_OPTS["mafft_opts"]
    threads: 8
    resources:
        mem_mb=RES_HEAVY["mem_mb"],
        runtime=RES_HEAVY["runtime_min"]
    log: log_of("mafft_align")
    benchmark: bench_of("mafft_align")
    conda: "envs/mafft.yaml"
    shell:
        r"""
        mkdir -p {ALIGN_DIR}
        mafft --thread {threads} {params.opts} {input} > {output} 2> {log}
        """

# --- Mask alignment ------------------------------------------------------------
rule mask_alignment:
    input:  aln = f"{ALIGN_DIR}/dnaA.aln.fasta"
    output: masked = f"{ALIGN_DIR}/dnaA.masked.fasta"
    params: gap = MASKING["gap_threshold"]
    threads: 2
    resources:
        mem_mb=RES_DEFAULT["mem_mb"],
        runtime=RES_DEFAULT["runtime_min"]
    log: log_of("mask_alignment")
    benchmark: bench_of("mask_alignment")
    conda: "envs/base-python.yaml"
    shell:
        r"""
        python scripts/python/mask_alignment.py \
          --in {input.aln} \
          --out {output.masked} \
          --mask-gap-frac {params.gap} \
          > {log} 2>&1
        """

# --- HyPhy GARD (produit un JSON) ---------------------------------------------
# Déclare ce chemin en haut de ton Snakefile avec tes autres constantes :
# en haut du Snakefile, avec tes autres constantes
HYPHY_ANALYSES = "external/hyphy-analyses"

rule hyphy_gard:
    input:
        aln = f"{ALIGN_DIR}/dnaA.masked.fasta",
        analyses = HYPHY_ANALYSES
    output:
        json = f"{RECOMB_DIR}/dnaA.gard.json"
    threads: 8
    log: log_of("hyphy_gard")
    benchmark: bench_of("hyphy_gard")
    conda: "envs/hyphy.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {RECOMB_DIR}
        TMPDIR=$(mktemp -d)
        OUTJSON="$TMPDIR/GARD.json"

        # 1) pointer HyPhy vers les analyses clonées
        export HYPHY_ANALYSES="{HYPHY_ANALYSES}"

        # 2) choisir le binaire (dans l'env conda)
        if command -v hyphy >/dev/null 2>&1; then
          H=hyphy
        elif command -v hyphy-avx >/dev/null 2>&1; then
          H=hyphy-avx
        else
          echo "[ERROR] HyPhy introuvable dans l'env conda" > {log}
          rm -rf "$TMPDIR"; exit 127
        fi

        echo "[which] $($H --version | head -n1)" > {log}
        echo "[env] HYPHY_ANALYSES=$HYPHY_ANALYSES" >> {log}

        # 3) GARD : ICI on passe un FICHIER en --output, pas un dossier
        if ! "$H" gard \
              --alignment {input.aln} \
              --type nucleotide \
              --output "$OUTJSON" \
              --threads {threads} >> {log} 2>&1; then
          echo "[ERROR] 'hyphy gard' a échoué" >> {log}
          rm -rf "$TMPDIR"; exit 1
        fi

        # 4) récupérer le JSON
        if [ -f "$OUTJSON" ]; then
          mv "$OUTJSON" {output.json}
        else
          cand=$(find "$TMPDIR" -maxdepth 2 -type f \( -name '*GARD*.json' -o -name 'GARD.json' \) | head -n1 || true)
          if [ -n "$cand" ]; then
            mv "$cand" {output.json}
          else
            echo "[ERROR] GARD.json introuvable dans $TMPDIR" >> {log}
            rm -rf "$TMPDIR"; exit 1
          fi
        fi

        rm -rf "$TMPDIR"
        """



# --- HyPhy sélection sur fragments (SLAC/FEL/FUBAR/MEME) ----------------------

rule make_fragments:
    input:
        msa = f"{ALIGN_DIR}/dnaA.masked.fasta",
        # tu peux choisir la source: consensus_breakpoints.tsv (idéal)…
        breaks = f"{RECOMB_DIR}/consensus_breakpoints.tsv"
        # … ou, si tu veux aller plus vite, {RECOMB_DIR}/dnaA.gard_breakpoints.tsv
    output:
        touch(f"{FRAG_DIR}/.done")
    threads: 1
    conda: "envs/base-python.yaml"
    log: log_of("make_fragments")
    benchmark: bench_of("make_fragments")
    shell:
        r"""
        mkdir -p {FRAG_DIR}
        python scripts/python/split_fragments.py \
          --aln {input.msa} \
          --breaks {input.breaks} \
          --outdir {FRAG_DIR} \
          --label dnaA \
          > {log} 2>&1
        touch {output}
        """


rule hyphy_batch_runner:
    input:
        fragments = FRAG_DIR,                  # <-- au lieu de l'aln
        treefile  = f"{TREE_DIR}/dnaA.treefile"  # optionnel mais malin
    output:
        touch("results/hyphy/.done")
    threads: 8
    conda: "envs/hyphy.yaml"
    log: log_of("hyphy_batch_runner")
    benchmark: bench_of("hyphy_batch_runner")
    shell:
        r"""
        python scripts/python/hyphy_batch_runner.py \
          --fragments {input.fragments} \
          --outdir results/hyphy \
          --tree-dir {TREE_DIR} \
          --build-tree \
          --cpu {threads} \
          --genetic-code 1 \
          --iqtree-model GTR+G \
          > {log} 2>&1
        touch {output}
        """



# --- Parse GARD JSON → TSV breakpoints ----------------------------------------
rule gard_to_breaks:
    input:
        json = f"{RECOMB_DIR}/dnaA.gard.json"
    output:
        tsv  = f"{RECOMB_DIR}/dnaA.gard_breakpoints.tsv"
    threads: 1
    resources:
        mem_mb=RES_DEFAULT["mem_mb"],
        runtime=RES_DEFAULT["runtime_min"]
    log: log_of("gard_to_breaks")
    benchmark: bench_of("gard_to_breaks")
    conda: "envs/base-python.yaml"
    shell:
        r"""
        mkdir -p {RECOMB_DIR}
        python scripts/python/gard_to_breaks.py \
          --in {input.json} \
          --out {output.tsv} \
          > {log} 2>&1
        """

# --- Combine events (scaffold) -------------------------------------------------
rule combine_events:
    input:  gard = f"{RECOMB_DIR}/dnaA.gard_breakpoints.tsv"
    output: combined = f"{RECOMB_DIR}/events_all.tsv"
    threads: 1
    resources:
        mem_mb=RES_DEFAULT["mem_mb"],
        runtime=RES_DEFAULT["runtime_min"]
    log: log_of("combine_events")
    benchmark: bench_of("combine_events")
    conda: "envs/base-python.yaml"
    shell:
        r"""
        # Placeholder: replace with combine_rdp_events.py when other detectors are added.
        cp {input.gard} {output.combined}
        echo "Combined events currently equal to GARD-only breaks." >> {log}
        """

# --- Consensus breakpoints -----------------------------------------------------
rule consensus_breakpoints:
    input:  combined = f"{RECOMB_DIR}/events_all.tsv"
    output: consensus = f"{RECOMB_DIR}/consensus_breakpoints.tsv"
    params: seed = SEEDS["consensus"]
    threads: 1
    resources:
        mem_mb=RES_DEFAULT["mem_mb"],
        runtime=RES_DEFAULT["runtime_min"]
    log: log_of("consensus_breakpoints")
    benchmark: bench_of("consensus_breakpoints")
    conda: "envs/base-python.yaml"
    shell:
        r"""
        python scripts/python/consensus_breakpoints.py \
          --in {input.combined} \
          --out {output.consensus} \
          --seed {params.seed} \
          > {log} 2>&1
        """

# --- Parentage assignment (BLAST on non-recombinant segments) -----------------
rule parentage_blast:
    input:
        msa       = f"{ALIGN_DIR}/dnaA.masked.fasta",
        consensus = f"{RECOMB_DIR}/consensus_breakpoints.tsv"
    output:
        calls = f"{PARENT_DIR}/parent_calls.tsv"
    params:
        seed     = SEEDS["blast"],
        refs_dir = REFS_DIR
    threads: 8
    resources:
        mem_mb=RES_HEAVY["mem_mb"],
        runtime=RES_HEAVY["runtime_min"]
    log: log_of("parentage_blast")
    benchmark: bench_of("parentage_blast")
    conda: "envs/blast.yaml"
    shell:
        r"""
        mkdir -p {PARENT_DIR}
        python - << 'PY' > {log} 2>&1
import os, glob, subprocess, tempfile, shutil, random, pandas as pd
random.seed({params.seed})

msa_path = "{input.msa}"
cons_path = "{input.consensus}"
refs_dir = "{params.refs_dir}"
out_tsv  = "{output.calls}"

def read_fasta(path):
    seqs=[]; hdr=None; buf=[]
    with open(path) as fh:
        for line in fh:
            line=line.rstrip()
            if not line: continue
            if line.startswith(">"):
                if hdr is not None: seqs.append((hdr,"".join(buf)))
                hdr=line[1:].split()[0]; buf=[]
            else:
                buf.append(line)
        if hdr is not None: seqs.append((hdr,"".join(buf)))
    return seqs

def read_breaks(consensus_path):
    # expect a TSV; try to pull numeric breakpoints from any column
    bps=set()
    with open(consensus_path) as fh:
        header = fh.readline()
        for line in fh:
            toks=line.replace(",","\t").split()
            for t in toks:
                try:
                    v=int(t)
                    if v>0: bps.add(v)
                except: pass
    bps=sorted(bps)
    return bps

def slice_columns(msa, start, end):
    s=start-1; e=end
    return [(sid, seq[s:e]) for sid,seq in msa]

def consensus_of_slice(sliced):
    if not sliced: return ""
    L=len(sliced[0][1]); out=[]
    for i in range(L):
        freq={}
        for sid,s in sliced:
            c=s[i]
            freq[c]=freq.get(c,0)+1
        out.append(sorted(freq.items(), key=lambda kv:(-kv[1],kv[0]))[0][0])
    return "".join(out)

def build_blast_db_if_any(refs_dir, tmpdir):
    fa_list=sorted(glob.glob(os.path.join(refs_dir,"*.fa"))+
                   glob.glob(os.path.join(refs_dir,"*.fasta"))+
                   glob.glob(os.path.join(refs_dir,"*.fna")))
    if not fa_list: return None
    cat=os.path.join(tmpdir,"refs.cat.fasta")
    with open(cat,"w") as out:
        for fp in fa_list:
            with open(fp) as fh: out.write(fh.read())
    db=os.path.join(tmpdir,"refdb")
    subprocess.run(["makeblastdb","-in",cat,"-dbtype","nucl","-out",db],
                   check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return db

msa=read_fasta(msa_path)
if not msa: raise SystemExit("Masked alignment is empty.")
L=len(msa[0][1])

bps=read_breaks(cons_path)
cut_points=[bp for bp in bps if 1 <= bp < L]

segments=[]
prev=1
for bp in cut_points + [L]:
    segments.append((prev,bp))
    prev=bp+1

tmpdir=tempfile.mkdtemp(prefix="blast_parentage_")
db=None
try:
    have_refs=os.path.isdir(refs_dir)
    if have_refs:
        try: db=build_blast_db_if_any(refs_dir,tmpdir)
        except Exception as e: print("[WARN] build DB:", e)

    rows=[]
    for idx,(a,b) in enumerate(segments, start=1):
        sliced=slice_columns(msa,a,b)
        cons=consensus_of_slice(sliced)
        qfa=os.path.join(tmpdir,f"seg_{idx:03d}.fa")
        with open(qfa,"w") as o:
            o.write(f">seg_{idx:03d}_{a}_{b}\n{cons}\n")

        major="NA"; minor="NA"; score="NA"
        if db and cons:
            cmd=["blastn","-query",qfa,"-db",db,
                 "-outfmt","6 qseqid sseqid pident length evalue bitscore qcovs",
                 "-max_target_seqs","5"]
            try:
                cp=subprocess.run(cmd,check=True,text=True,capture_output=True)
                lines=[ln for ln in cp.stdout.splitlines() if ln.strip()]
                if lines:
                    best=lines[0].split("\t")
                    if len(best)>=7:
                        major=best[1]; score=best[6]
                    if len(lines)>1:
                        second=lines[1].split("\t")
                        if len(second)>=2: minor=second[1]
            except Exception as e:
                print(f"[WARN] BLAST seg {idx}: {e}")

        rows.append((f"{a}-{b}", major, minor, score))

    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    with open(out_tsv,"w") as out:
        out.write("breakpoint_interval\tmajor_parent\tminor_parent\tscore\n")
        for r in rows: out.write("\t".join(r)+"\n")
finally:
    try: shutil.rmtree(tmpdir, ignore_errors=True)
    except: pass
PY
        """

# --- NeighborNet network (R) --------------------------------------------------
rule neighbornet:
    input:
        msa  = f"{ALIGN_DIR}/dnaA.masked.fasta",
        meta = METADATA_TSV
    output:
        pdf    = f"{NET_DIR}/neighborNet.pdf",
        png    = f"{NET_DIR}/neighborNet.png",
        colors = f"{NET_DIR}/taxa_colors.csv"
    params:
        max_tips = config["network"]["max_tips"],
        per_cap  = config["network"]["per_species_cap"],
        seed     = SEEDS["network"],
        pdf_w    = config["network"]["pdf_w"],
        pdf_h    = config["network"]["pdf_h"]
    threads: 2
    resources:
        mem_mb=RES_DEFAULT["mem_mb"],
        runtime=RES_DEFAULT["runtime_min"]
    log: log_of("neighbornet")
    benchmark: bench_of("neighbornet")
    conda: "envs/r-env.yaml"
    shell:
        r"""
        mkdir -p {NET_DIR}
        Rscript scripts/R/neighborNet.R \
          --aln {input.msa} \
          --out {NET_DIR}/neighborNet \
          --metadata {input.meta} \
          --color-map {output.colors} \
          --max-tips {params.max_tips} \
          --per-species-cap {params.per_cap} \
          --seed {params.seed} \
          --pdf-width {params.pdf_w} \
          --pdf-height {params.pdf_h} \
          > {log} 2>&1
        """

# --- IQ-TREE ML phylogeny -----------------------------------------------------
rule iqtree_tree:
    input:  msa = f"{ALIGN_DIR}/dnaA.masked.fasta"
    output:
        treefile = f"{TREE_DIR}/dnaA.treefile",
        iqtree   = f"{TREE_DIR}/dnaA.iqtree",
        logf     = f"{TREE_DIR}/dnaA.log"
    params:
        seed = SEEDS["iqtree"],
        opts = IQTREE_OPTS            # ex: "-m GTR+G -bb 1000 -alrt 1000"
    threads: 8
    resources:
        mem_mb=RES_HEAVY["mem_mb"],
        runtime=RES_HEAVY["runtime_min"]
    log: log_of("iqtree_tree")
    benchmark: bench_of("iqtree_tree")
    conda: "envs/iqtree.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {TREE_DIR}
        # choisir binaire
        if command -v iqtree >/dev/null 2>&1; then
          IQ=iqtree
        elif command -v iqtree2 >/dev/null 2>&1; then
          IQ=iqtree2
        else
          echo "[ERROR] IQ-TREE introuvable dans l'environnement conda." >&2
          echo "Assure-toi que envs/iqtree.yaml installe 'iqtree' ou 'iqtree2'." >&2
          exit 127
        fi

        echo "[$($IQ -version | head -n1)]" > {log} 2>&1
        "$IQ" -s {input.msa} -seed {params.seed} -nt {threads} \
              -pre {TREE_DIR}/dnaA \
              -redo \
              {params.opts} >> {log} 2>&1
        """

# --- Final summary -------------------------------------------------------------
rule final_summary:
    input:
        consensus = f"{RECOMB_DIR}/consensus_breakpoints.tsv",
        calls     = f"{PARENT_DIR}/parent_calls.tsv"
    output:
        summary = f"{SUM_DIR}/recombination_summary.tsv"
    threads: 1
    resources:
        mem_mb=RES_DEFAULT["mem_mb"],
        runtime=RES_DEFAULT["runtime_min"]
    log: log_of("final_summary")
    benchmark: bench_of("final_summary")
    conda: "envs/base-python.yaml"
    shell:
        r"""
        mkdir -p {SUM_DIR}
        python - << 'PY' > {log} 2>&1
import pandas as pd
cons = pd.read_csv("{input.consensus}", sep="\t")
par  = pd.read_csv("{input.calls}", sep="\t")

if "breakpoint_interval" not in par.columns:
    par["breakpoint_interval"] = "NA"
if "event_id" not in cons.columns:
    cons["event_id"] = range(1, len(cons)+1)

join_cols = [c for c in ["breakpoint_interval","interval","range"] if c in cons.columns and c in par.columns]
if not join_cols:
    cons["key"] = cons.index.astype(str)
    par["key"] = "all"
    out = cons.merge(par, on="key", how="left")
else:
    key = join_cols[0]
    out = cons.merge(par, left_on=key, right_on=key, how="left")

out.to_csv("{output.summary}", sep="\t", index=False)
print("Wrote summary with", len(out), "rows.")
PY
        """

# --- Rule ordering to disambiguate wildcards ----------------------------------
ruleorder: hyphy_gard > gard_to_breaks > combine_events > consensus_breakpoints > parentage_blast > final_summary

################################################################################
# End of Snakefile
################################################################################

