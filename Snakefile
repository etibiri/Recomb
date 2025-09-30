################################################################################
# Recombination Analysis
#
# Reproducible, offline-ready Snakemake workflow.
# Only uses the existing Conda envs and the provided scripts.
#
# Envs (exist in envs/): base-python.yaml, blast.yaml, hyphy.yaml,
#                             iqtree.yaml, mafft.yaml, r-env.yaml
#
# Provided scripts (exist in scripts/):
#   neighborNet.R
#   combine_rdp_events.py
#   consensus_breakpoints.py
#   curate_sequences.py
#   gard_to_breaks.py
#   hyphy_batch_runner.py
#   mask_alignment.py
#
# Key outputs (see config and “rule all”):
#   - results/curation/dnaA.cleaned.fasta, results/curation/qc.tsv
#   - results/alignments/dnaA.aln.fasta, results/alignments/dnaA.masked.fasta
#   - results/recombination/gard_raw.json, results/recombination/gard_breakpoints.tsv
#   - results/recombination/events_all.tsv, results/recombination/consensus_breakpoints.tsv
#   - results/parentage/parent_calls.tsv
#   - results/phylogeny/iqtree/* (tree + logs)
#   - results/network/neighborNet.(pdf|png) (+ colors/annotations)
#   - results/summary/recombination_summary.tsv
#   - results/logs/* (per-rule logs) and results/benchmarks/* (per-rule metrics)
#
# Determinism: fixed seeds pulled from config. All I/O paths are stable.
# HPC-friendly: per-rule threads/resources/log/benchmark.
################################################################################

import os
import textwrap

configfile: "config.yaml"

# --- Shortcuts from config -----------------------------------------------------

DATA_FASTA         = config["inputs"]["dnaA_fasta"]
METADATA_TSV       = config["inputs"]["metadata"]
REFS_DIR           = config["inputs"]["refs_dir"]

SEEDS              = config["seeds"]
MAFFT_OPTS         = config["tools"]["mafft_opts"]
IQTREE_OPTS        = config["tools"]["iqtree_opts"]
HYPHY_OPTS         = config["tools"]["hyphy_opts"]
MASKING            = config["masking"]

RES_DEFAULT        = config["resources"]["default"]
RES_HEAVY          = config["resources"]["heavy"]

LOGDIR             = "results/logs"
BMDIR              = "results/benchmarks"

# ensure log/benchmark directories exist when dry-running
os.makedirs(LOGDIR, exist_ok=True)
os.makedirs(BMDIR, exist_ok=True)

# --- Helper: consistent logs/benchmarks paths ---------------------------------

def log_of(rule_name):
    return f"{LOGDIR}/{rule_name}.log"

def bench_of(rule_name):
    return f"{BMDIR}/{rule_name}.tsv"

# --- Rule: all ----------------------------------------------------------------
# Collects all final deliverables required by §3 of the spec.

rule all:
    message:
        "Aggregating all final outputs for the CMD DNA-A recombination analysis."
    input:
        # Curation
        "results/curation/dnaA.cleaned.fasta",
        "results/curation/qc.tsv",
        # Alignments
        "results/alignments/dnaA.aln.fasta",
        "results/alignments/dnaA.masked.fasta",
        # Recombination detection & consensus
        "results/recombination/gard_raw.json",
        "results/recombination/gard_breakpoints.tsv",
        "results/recombination/events_all.tsv",
        "results/recombination/consensus_breakpoints.tsv",
        # Parentage
        "results/parentage/parent_calls.tsv",
        # Phylogeny
        "results/phylogeny/iqtree/dnaA.treefile",
        "results/phylogeny/iqtree/dnaA.iqtree",
        "results/phylogeny/iqtree/dnaA.log",
        # Network
        "results/network/neighborNet.pdf",
        "results/network/neighborNet.png",
        "results/network/taxa_colors.csv",
        # Summary
        "results/summary/recombination_summary.tsv"
    threads: 1
    resources:
        mem_mb=RES_DEFAULT["mem_mb"],
        runtime=RES_DEFAULT["runtime_min"]
    log: log_of("all")
    benchmark: bench_of("all")
    conda:
        "envs/base-python.yaml"
    shell:
        r"""
        echo "[OK] Final targets exist." > {log}
        """

# --- Sequence curation & normalization ----------------------------------------
# Uses scripts/curate_sequences.py on data/dnaA.fasta.
# Produces cleaned FASTA and a QC report.

rule curate_sequences:
    input:
        fasta = DATA_FASTA
    output:
        cleaned = "results/curation/dnaA.cleaned.fasta",
        qc = "results/curation/qc.tsv"
    threads: 2
    resources:
        mem_mb=RES_DEFAULT["mem_mb"],
        runtime=RES_DEFAULT["runtime_min"]
    log: log_of("curate_sequences")
    benchmark: bench_of("curate_sequences")
    conda:
        "envs/base-python.yaml"
    shell:
        r"""
        mkdir -p results/curation
        python scripts/curate_sequences.py \
          --in {input.fasta} \
          --out {output.cleaned} \
          --qc {output.qc} \
          --seed {SEEDS[cleaning]} \
          > {log} 2>&1
        """

# --- Multiple sequence alignment (MAFFT) --------------------------------------
# Align cleaned sequences using MAFFT with options from config.
# Output: results/alignments/dnaA.aln.fasta

rule mafft_align:
    input:
        "results/curation/dnaA.cleaned.fasta"
    output:
        "results/alignments/dnaA.aln.fasta"
    threads: 8
    resources:
        mem_mb=RES_HEAVY["mem_mb"],
        runtime=RES_HEAVY["runtime_min"]
    log: log_of("mafft_align")
    benchmark: bench_of("mafft_align")
    conda:
        "envs/mafft.yaml"
    shell:
        r"""
        mkdir -p results/alignments
        mafft --thread {threads} --seed {SEEDS[mafft]} {MAFFT_OPTS} {input} > {output} 2> {log}
        """

# --- Mask alignment columns (gap-rich/low-informative) ------------------------
# Uses scripts/mask_alignment.py to mask the MAFFT alignment.
# Output: results/alignments/dnaA.masked.fasta

rule mask_alignment:
    input:
        aln = "results/alignments/dnaA.aln.fasta"
    output:
        masked = "results/alignments/dnaA.masked.fasta"
    threads: 2
    resources:
        mem_mb=RES_DEFAULT["mem_mb"],
        runtime=RES_DEFAULT["runtime_min"]
    log: log_of("mask_alignment")
    benchmark: bench_of("mask_alignment")
    conda:
        "envs/base-python.yaml"
    shell:
        r"""
        python scripts/mask_alignment.py \
          --in {input.aln} \
          --out {output.masked} \
          --gap-thresh {MASKING[gap_threshold]} \
          --treat-ambiguous-as-gap {MASKING[treat_ambiguous_as_gap]} \
          --seed {SEEDS[masking]} \
          > {log} 2>&1
        """

# --- Recombination detection (HyPhy GARD) -------------------------------------
# Runs provided hyphy_batch_runner.py on the masked MSA.
# Outputs raw GARD JSON.

rule hyphy_gard:
    input:
        msa = "results/alignments/dnaA.masked.fasta"
    output:
        json = "results/recombination/gard_raw.json"
    threads: 8
    resources:
        mem_mb=RES_HEAVY["mem_mb"],
        runtime=RES_HEAVY["runtime_min"]
    log: log_of("hyphy_gard")
    benchmark: bench_of("hyphy_gard")
    conda:
        "envs/hyphy.yaml"
    shell:
        r"""
        mkdir -p results/recombination
        python scripts/hyphy_batch_runner.py \
          --msa {input.msa} \
          --out {output.json} \
          --seed {SEEDS[hyphy]} \
          {HYPHY_OPTS} \
          > {log} 2>&1
        """

# --- Convert GARD JSON to tabular breakpoints ---------------------------------
# Uses scripts/gard_to_breaks.py to produce a standardized TSV of breakpoints.

rule gard_to_breaks:
    input:
        json = "results/recombination/gard_raw.json"
    output:
        tsv = "results/recombination/gard_breakpoints.tsv"
    threads: 1
    resources:
        mem_mb=RES_DEFAULT["mem_mb"],
        runtime=RES_DEFAULT["runtime_min"]
    log: log_of("gard_to_breaks")
    benchmark: bench_of("gard_to_breaks")
    conda:
        "envs/base-python.yaml"
    shell:
        r"""
        python scripts/gard_to_breaks.py \
          --in {input.json} \
          --out {output.tsv} \
          > {log} 2>&1
        """

# --- (Optional) Event consolidation scaffold ----------------------------------
# If other detectors are added later, combine them here with combine_rdp_events.py.
# For now, we simply pass-through GARD breaks as events_all.tsv for a stable interface.

rule combine_events:
    input:
        gard = "results/recombination/gard_breakpoints.tsv"
    output:
        combined = "results/recombination/events_all.tsv"
    threads: 1
    resources:
        mem_mb=RES_DEFAULT["mem_mb"],
        runtime=RES_DEFAULT["runtime_min"]
    log: log_of("combine_events")
    benchmark: bench_of("combine_events")
    conda:
        "envs/base-python.yaml"
    shell:
        r"""
        # If future detectors exist, replace this pass-through with:
        # python scripts/combine_rdp_events.py --inputs ... --out {output.combined}
        cp {input.gard} {output.combined}
        echo "Combined events currently equal to GARD-only breaks." >> {log}
        """

# --- Consensus breakpoints -----------------------------------------------------
# Generates a consensus set of breakpoints (robust across detectors).
# Here we run consensus_breakpoints.py over events_all.tsv (GARD-only by default).

rule consensus_breakpoints:
    input:
        combined = "results/recombination/events_all.tsv"
    output:
        consensus = "results/recombination/consensus_breakpoints.tsv"
    threads: 1
    resources:
        mem_mb=RES_DEFAULT["mem_mb"],
        runtime=RES_DEFAULT["runtime_min"]
    log: log_of("consensus_breakpoints")
    benchmark: bench_of("consensus_breakpoints")
    conda:
        "envs/base-python.yaml"
    shell:
        r"""
        python scripts/consensus_breakpoints.py \
          --in {input.combined} \
          --out {output.consensus} \
          --seed {SEEDS[consensus]} \
          > {log} 2>&1
        """

# --- Putative parentage assignment (BLAST on non-recombinant segments) --------
# Reads consensus breakpoints and the masked alignment, splits the alignment into
# non-recombinant segments (alignment positions), and assigns putative parents
# by BLASTing segment consensus against refs (if available).
# Output: results/parentage/parent_calls.tsv
#
# Notes:
#  - If data/refs/ is empty or missing, the rule still runs and returns NA calls.
#  - Pure-Python FASTA parsing used for segment slicing (no extra deps).
#  - All logic encapsulated here for deterministic single-output behavior.

rule parentage_blast:
    input:
        msa = "results/alignments/dnaA.masked.fasta",
        consensus = "results/recombination/consensus_breakpoints.tsv"
    output:
        calls = "results/parentage/parent_calls.tsv"
    threads: 8
    resources:
        mem_mb=RES_HEAVY["mem_mb"],
        runtime=RES_HEAVY["runtime_min"]
    log: log_of("parentage_blast")
    benchmark: bench_of("parentage_blast")
    conda:
        "envs/blast.yaml"
    shell:
        r"""
        mkdir -p results/parentage

        python - << 'PY' > {log} 2>&1
import os, sys, glob, subprocess, tempfile, shutil, textwrap, random

random.seed({SEEDS[blast]})

msa_path = "{input.msa}"
cons_path = "{input.consensus}"
refs_dir = "{REFS_DIR}"
out_tsv  = "{output.calls}"

def read_fasta(path):
    seqs = []
    with open(path) as fh:
        hdr = None
        buf = []
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if hdr is not None:
                    seqs.append((hdr, "".join(buf)))
                hdr = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if hdr is not None:
            seqs.append((hdr, "".join(buf)))
    return seqs  # list[(id, seq)]

def read_breaks(consensus_path):
    # Expecting a TSV with breakpoint positions (alignment-based).
    # We accept any header; we look for integer positions in columns named 'pos' or 'break' or 'position',
    # OR two-column '(start, end)' windows. To keep generic, we parse all ints per line.
    bps = []
    with open(consensus_path) as fh:
        header = fh.readline().rstrip("\n")
        for line in fh:
            if not line.strip():
                continue
            ints = []
            for tok in line.strip().replace(",", "\t").split():
                try:
                    ints.append(int(tok))
                except ValueError:
                    pass
            for x in ints:
                if x not in bps:
                    bps.append(x)
    bps = sorted(set([bp for bp in bps if bp > 0]))
    return bps

def slice_columns(msa_seqs, start, end):
    # start/end are 1-based inclusive alignment column coordinates
    s = start-1; e = end
    sliced = []
    for sid, seq in msa_seqs:
        sliced.append((sid, seq[s:e]))
    return sliced

def consensus_of_slice(sliced):
    # Very simple majority consensus (ACGTN-? kept as-is if tie)
    L = len(sliced[0][1]) if sliced and sliced[0][1] is not None else 0
    cons = []
    for i in range(L):
        col = [s[1][i] for s in sliced]
        freq = {}
        for c in col:
            freq[c] = freq.get(c, 0) + 1
        best = sorted(freq.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
        cons.append(best)
    return "".join(cons)

def build_blast_db_if_any(refs_dir, tmpdir):
    fa_list = sorted(glob.glob(os.path.join(refs_dir, "*.fa")) +
                     glob.glob(os.path.join(refs_dir, "*.fasta")) +
                     glob.glob(os.path.join(refs_dir, "*.fna")))
    if not fa_list:
        return None
    cat = os.path.join(tmpdir, "refs.cat.fasta")
    with open(cat, "w") as out:
        for fp in fa_list:
            with open(fp) as fh:
                out.write(fh.read())
    db = os.path.join(tmpdir, "refdb")
    subprocess.run(["makeblastdb", "-in", cat, "-dbtype", "nucl", "-out", db],
                   check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return db

msa = read_fasta(msa_path)
if not msa:
    raise SystemExit("Masked alignment is empty.")

bps = read_breaks(cons_path)
L = len(msa[0][1])

# define non-recombinant segments as: [1 .. bp1], (bp1+1 .. bp2], ... (bpK+1 .. L]
cut_points = [bp for bp in bps if 1 <= bp < L]
segments = []
prev = 1
for bp in cut_points + [L]:
    segments.append( (prev, bp) )
    prev = bp + 1

tmpdir = tempfile.mkdtemp(prefix="blast_parentage_")
db = None
try:
    have_refs = os.path.isdir(refs_dir)
    if have_refs:
        try:
            db = build_blast_db_if_any(refs_dir, tmpdir)
        except Exception as e:
            print(f"[WARN] Failed to build BLAST DB from {refs_dir}: {e}")

    rows = []
    seg_id = 0
    for (a,b) in segments:
        seg_id += 1
        sliced = slice_columns(msa, a, b)
        cons = consensus_of_slice(sliced)
        # write query fasta
        qfa = os.path.join(tmpdir, f"seg_{seg_id:03d}.fa")
        with open(qfa, "w") as out:
            out.write(f">seg_{seg_id:03d}_{a}_{b}\n{cons}\n")

        major = "NA"
        minor = "NA"
        score = "NA"

        if db is not None:
            # run blastn, take top hit as 'major_parent' (simple heuristic)
            cmd = [
                "blastn", "-query", qfa, "-db", db,
                "-outfmt", "6 qseqid sseqid pident length evalue bitscore qcovs",
                "-max_target_seqs", "5"
            ]
            try:
                cp = subprocess.run(cmd, check=True, text=True, capture_output=True)
                lines = [ln for ln in cp.stdout.splitlines() if ln.strip()]
                if lines:
                    best = lines[0].split("\t")
                    # best[1] = subject id; best[6] = qcovs
                    major = best[1]
                    score = best[6]  # percent query coverage as a proxy
                    # optional: set minor if a second distinct hit exists
                    if len(lines) > 1:
                        minor = lines[1].split("\t")[1]
            except Exception as e:
                print(f"[WARN] BLAST failed on segment {seg_id}: {e}")

        rows.append( (f"{a}-{b}", major, minor, score) )

    with open(out_tsv, "w") as out:
        out.write("breakpoint_interval\tmajor_parent\tminor_parent\tscore\n")
        for r in rows:
            out.write("\t".join(map(str, r)) + "\n")

finally:
    shutil.rmtree(tmpdir, ignore_errors=True)
PY
        """

# --- IQ-TREE ML phylogeny -----------------------------------------------------
# Runs IQ-TREE on masked alignment, with deterministic seed and config-specified options.
# Outputs with -pre results/phylogeny/iqtree/dnaA.*

rule iqtree_tree:
    input:
        msa = "results/alignments/dnaA.masked.fasta"
    output:
        treefile = "results/phylogeny/iqtree/dnaA.treefile",
        iqtree   = "results/phylogeny/iqtree/dnaA.iqtree",
        logf     = "results/phylogeny/iqtree/dnaA.log"
    threads: 8
    resources:
        mem_mb=RES_HEAVY["mem_mb"],
        runtime=RES_HEAVY["runtime_min"]
    log: log_of("iqtree_tree")
    benchmark: bench_of("iqtree_tree")
    conda:
        "envs/iqtree.yaml"
    shell:
        r"""
        mkdir -p results/phylogeny/iqtree
        iqtree2 -s {input.msa} -seed {SEEDS[iqtree]} -nt {threads} \
                -pre results/phylogeny/iqtree/dnaA \
                {IQTREE_OPTS} > {log} 2>&1
        """

# --- NeighborNet network (R) --------------------------------------------------
# Uses scripts/neighborNet.R with metadata and plot settings from config.
# Produces PDF, PNG, colors CSV; optionally an RDS if the script supports it.

rule neighbornet:
    input:
        msa = "results/alignments/dnaA.masked.fasta",
        meta = METADATA_TSV
    output:
        pdf   = "results/network/neighborNet.pdf",
        png   = "results/network/neighborNet.png",
        colors= "results/network/taxa_colors.csv"
    threads: 2
    resources:
        mem_mb=RES_DEFAULT["mem_mb"],
        runtime=RES_DEFAULT["runtime_min"]
    log: log_of("neighbornet")
    benchmark: bench_of("neighbornet")
    conda:
        "envs/r-env.yaml"
    shell:
        r"""
        mkdir -p results/network
        Rscript scripts/neighborNet.R \
          --aln {input.msa} \
          --out results/network/neighborNet \
          --metadata {input.meta} \
          --color-map {output.colors} \
          --max-tips {config[network][max_tips]} \
          --per-species-cap {config[network][per_species_cap]} \
          --seed {SEEDS[network]} \
          --pdf-width {config[network][pdf_w]} \
          --pdf-height {config[network][pdf_h]} \
          > {log} 2>&1
        """

# --- Final summary table -------------------------------------------------------
# Collates consensus breakpoints and parentage calls into a single TSV.

rule final_summary:
    input:
        consensus = "results/recombination/consensus_breakpoints.tsv",
        calls     = "results/parentage/parent_calls.tsv"
    output:
        summary = "results/summary/recombination_summary.tsv"
    threads: 1
    resources:
        mem_mb=RES_DEFAULT["mem_mb"],
        runtime=RES_DEFAULT["runtime_min"]
    log: log_of("final_summary")
    benchmark: bench_of("final_summary")
    conda:
        "envs/base-python.yaml"
    shell:
        r"""
        mkdir -p results/summary
        python - << 'PY' > {log} 2>&1
import pandas as pd

cons = pd.read_csv("{input.consensus}", sep="\t")
par  = pd.read_csv("{input.calls}", sep="\t")

# ensure expected columns exist; if not, fall back gracefully
if "breakpoint_interval" not in par.columns:
    par["breakpoint_interval"] = "NA"
if "event_id" not in cons.columns:
    cons["event_id"] = range(1, len(cons)+1)

# naive merge on nearest / shared interval-style key if present
join_cols = [c for c in ["breakpoint_interval","interval","range"] if c in cons.columns and c in par.columns]
if not join_cols:
    cons["key"] = cons.index.astype(str)
    par["key"] = "all"
    out = cons.merge(par, left_on="key", right_on="key", how="left")
else:
    key = join_cols[0]
    out = cons.merge(par, left_on=key, right_on=key, how="left")

out.to_csv("{output.summary}", sep="\t", index=False)
print("Wrote summary with", len(out), "rows.")
PY
        """

# --- Simple phony dependency graph to order steps ------------------------------

# Curation -> Align -> Mask
use rule curate_sequences
use rule mafft_align as mafft_align
use rule mask_alignment as mask_alignment

# Mask -> HyPhy -> Breaks -> Combine -> Consensus -> Parentage -> Summary
ruleorder: hyphy_gard > gard_to_breaks > combine_events > consensus_breakpoints > parentage_blast > final_summary

# Additional chaining (Snakemake infers most via inputs/outputs):
#   neighborNet and iqtree_tree both depend on masked alignment.


################################################################################
# End of Snakefile
################################################################################

