################################################################################
# Recombination Analysis — CMD DNA-A
# Author: Ezechiel B. TIBIRI
# Reproducible, offline-ready Snakemake workflow.
# Uses existing Conda envs (envs/) and provided scripts (scripts/).
################################################################################

import os

configfile: "config.yaml"
shell.executable("/bin/bash")
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

# à mettre une seule fois en haut du Snakefile si absent
#shell.executable("/bin/bash")

rule hyphy_gard:
    input:
        aln = f"{ALIGN_DIR}/dnaA.masked.fasta"
    output:
        json = f"{RECOMB_DIR}/dnaA.gard.json"
    params:
        extra = HYPHY_OPTS  # peut être vide
    threads: 8
    log: log_of("hyphy_gard")
    benchmark: bench_of("hyphy_gard")
    conda: "envs/hyphy.yaml"
    shell:
        r"""
        set -euo pipefail

        outdir="$(dirname "{output.json}")"
        mkdir -p "$outdir"

        TMPDIR="$(mktemp -d)"
        OUTJSON="$TMPDIR/GARD.json"

        # définir HYPHY_ANALYSES si absent, sans expansions fragiles
        HVAL="$(printenv HYPHY_ANALYSES 2>/dev/null || true)"
        if [ -z "$HVAL" ] && [ -d "external/hyphy-analyses" ]; then
          export HYPHY_ANALYSES=external/hyphy-analyses
        fi

        export TOLERATE_NUMERICAL_ERRORS=1

        if command -v hyphy >/dev/null 2>&1; then H=hyphy
        elif command -v hyphy-avx >/dev/null 2>&1; then H=hyphy-avx
        else
          echo "[ERROR] HyPhy introuvable dans l'env conda" > "{log}"
          rm -rf "$TMPDIR"; exit 127
        fi

        $H --version 2>/dev/null | head -n1 | sed 's/^/[which] /' > "{log}" 2>&1
        echo "[env] HYPHY_ANALYSES=$(printenv HYPHY_ANALYSES 2>/dev/null || echo '')" >> "{log}"
        echo "[env] TOLERATE_NUMERICAL_ERRORS=$TOLERATE_NUMERICAL_ERRORS" >> "{log}"
        echo "[cmd] $H gard --alignment {input.aln} --type nucleotide --output $OUTJSON --threads {threads} --model GTR {params.extra}" >> "{log}"

        set +e
        $H gard \
          --alignment "{input.aln}" \
          --type nucleotide \
          --output "$OUTJSON" \
          --threads {threads} \
          --model {params.extra} \
           ENV=TOLERATE_NUMERICAL_ERRORS=1; >> "{log}" 2>&1
        rc=$?
        set -e

        if [ -s "$OUTJSON" ]; then
          mv "$OUTJSON" "{output.json}"
          echo "[OK] GARD JSON: {output.json}" >> "{log}"
          rm -rf "$TMPDIR"
          exit 0
        fi

        cand="$(find "$TMPDIR" -maxdepth 2 -type f -name '*GARD*.json' | head -n1 || true)"
        if [ -n "$cand" ] && [ -s "$cand" ]; then
          mv "$cand" "{output.json}"
          echo "[OK] GARD JSON (fallback): {output.json}" >> "{log}"
          rm -rf "$TMPDIR"
          exit 0
        fi

        echo "[ERROR] GARD.json introuvable; rc=$rc" >> "{log}"
        rm -rf "$TMPDIR"
        exit 1
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
    log: log_of("gard_to_breaks")
    conda: "envs/base-python.yaml"   # ou ton env python
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.tsv}")"
        scripts/python/gard_to_breaks.py "{input.json}" "{output.tsv}" "dnaA" > "{log}" 2>&1
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
    input:
        combined = f"{RECOMB_DIR}/events_all.tsv"
    output:
        consensus = f"{RECOMB_DIR}/consensus_breakpoints.tsv"
    threads: 1
    resources:
        mem_mb = RES_DEFAULT["mem_mb"],
        runtime = RES_DEFAULT["runtime_min"]
    log: log_of("consensus_breakpoints")
    benchmark: bench_of("consensus_breakpoints")
    conda: "envs/base-python.yaml"
    shell:
        r"""
        set -euo pipefail

        python - << 'PY' > "{log}" 2>&1
import os, re, json

inp  = r"{input.combined}"
outp = r"{output.consensus}"

if not os.path.isfile(inp):
    raise SystemExit("[ERROR] Missing input: %s" % inp)

# 1) Essayer d'extraire des breaks depuis events_all.tsv
breaks = []
gard_json_path = None

with open(inp, "r") as fh:
    for line in fh:
        line = line.strip()
        if not line:
            continue
        if line.startswith("#SUMMARY"):
            cols = line.split("\t")
            # format attendu: #SUMMARY, alignment, gard_json, n_breakpoints, breaks, best_aicc
            if len(cols) >= 3 and gard_json_path is None:
                gard_json_path = cols[2].strip()
            if len(cols) >= 5:
                raw = cols[4].strip()
                if raw:
                    for tok in re.split(r"[,\s]+", raw):
                        if tok.isdigit():
                            breaks.append(int(tok))
        elif line.startswith("#DETAIL"):
            cols = line.split("\t")
            if len(cols) >= 4 and cols[-1].isdigit():
                breaks.append(int(cols[-1]))

# 2) Si pas de breaks trouvés, tenter directement le JSON GARD référencé en SUMMARY
if not breaks and gard_json_path and os.path.isfile(gard_json_path):
    try:
        with open(gard_json_path, "r") as jh:
            data = json.load(jh)
        # Schéma HyPhy récent: improvements -> dernière étape contient "breakpoints"
        imp = data.get("improvements")
        if isinstance(imp, dict) and imp:
            last_key = max(int(k) for k in imp.keys())
            bp = imp[str(last_key)].get("breakpoints", [])
            flat = []
            for x in bp:
                if isinstance(x, list):
                    flat.extend(x)
                else:
                    flat.append(x)
            for v in flat:
                try:
                    breaks.append(int(v))
                except Exception:
                    pass
        # Backups éventuels pour anciens schémas
        if not breaks:
            for alt in ("breakpointList", "breakpoints", "bestBreakPoints"):
                vals = data.get(alt)
                if isinstance(vals, list):
                    for v in vals:
                        try:
                            breaks.append(int(v))
                        except Exception:
                            pass
    except Exception as e:
        print("[WARN] fallback JSON parse failed: %s" % str(e))

# 3) Dédupliquer et trier
breaks = sorted(set(breaks))

# 4) Fusion optionnelle de points proches en mini intervalles
MERGE_WINDOW = 10
intervals = []
for b in breaks:
    if not intervals:
        intervals.append([b, b])
    else:
        if b - intervals[-1][1] <= MERGE_WINDOW:
            intervals[-1][1] = b
        else:
            intervals.append([b, b])

# 5) Écrire sortie; support=1 car un seul détecteur
os.makedirs(os.path.dirname(outp), exist_ok=True)
with open(outp, "w") as out:
    out.write("start\tend\tsupport\n")
    for s, e in intervals:
        out.write("%d\t%d\t1\n" % (s, e))

print("[INFO] Consensus clusters: %d -> %s" % (len(intervals), outp))
PY
        """

# --- Parentage assignment (BLAST on non-recombinant segments) -----------------
rule parentage_blast:
    input:
        msa = f"{ALIGN_DIR}/dnaA.masked.fasta",
        consensus = f"{RECOMB_DIR}/consensus_breakpoints.tsv"
    output:
        calls = "results/parentage/parent_calls.tsv"
    params:
        seed = SEEDS["blast"],
        refs_dir = REFS_DIR
    threads: 8
    resources:
        mem_mb = RES_HEAVY["mem_mb"],
        runtime = RES_HEAVY["runtime_min"]
    log: log_of("parentage_blast")
    benchmark: bench_of("parentage_blast")
    conda: "envs/blast.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.calls}")"
        scripts/python/parentage_blast.py \
          --msa "{input.msa}" \
          --consensus "{input.consensus}" \
          --refs "{params.refs_dir}" \
          --out "{output.calls}" \
          --seed {params.seed} \
          --threads {threads} \
          > "{log}" 2>&1
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
        mem_mb = RES_DEFAULT["mem_mb"],
        runtime = RES_DEFAULT["runtime_min"]
    log: log_of("neighbornet")
    benchmark: bench_of("neighbornet")
    conda: "envs/r-env.yaml"
    shell:
        r"""
        set -euo pipefail

        outdir="$(dirname "{output.pdf}")"
        mkdir -p "$outdir"
        outbase="$outdir/neighborNet"

        # lancer le script R (aucune accolade non-Snakemake dans ce bloc)
        if ! Rscript scripts/R/neighborNet.R \
          --aln "{input.msa}" \
          --out "$outbase" \
          --metadata "{input.meta}" \
          --color-map "{output.colors}" \
          --max-tips {params.max_tips} \
          --per-species-cap {params.per_cap} \
          --seed {params.seed} \
          --pdf-width {params.pdf_w} \
          --pdf-height {params.pdf_h} \
          > "{log}" 2>&1; then
          echo "[WARN] R a retourné un code non nul; tentative de récupération des fichiers." >> "{log}"
        fi

        # rattrapage: chercher des fichiers plausibles et normaliser les noms
        if [ ! -s "{output.pdf}" ]; then
          cand_pdf="$(find "$outdir" -maxdepth 3 -type f \( -iname '*neighbor*net*.pdf' -o -iname '*network*.pdf' \) | head -n1 || true)"
          [ -n "$cand_pdf" ] && mv -f "$cand_pdf" "{output.pdf}" || true
        fi

        if [ ! -s "{output.png}" ]; then
          cand_png="$(find "$outdir" -maxdepth 3 -type f \( -iname '*neighbor*net*.png' -o -iname '*network*.png' \) | head -n1 || true)"
          [ -n "$cand_png" ] && mv -f "$cand_png" "{output.png}" || true
        fi

        if [ ! -s "{output.colors}" ]; then
          cand_col="$(find "$outdir" -maxdepth 3 -type f \( -iname '*taxa*color*.csv' -o -iname '*color*.csv' \) | head -n1 || true)"
          [ -n "$cand_col" ] && mv -f "$cand_col" "{output.colors}" || true
        fi

        # courte attente pour latence FS
        tries=0
        while [ $tries -lt 5 ]; do
          if [ -s "{output.pdf}" ] && [ -s "{output.png}" ] && [ -s "{output.colors}" ]; then
            echo "[OK] neighborNet outputs présents." >> "{log}"
            exit 0
          fi
          tries=$((tries+1))
          sleep 2
        done

        echo "[ERROR] Fichiers neighborNet manquants. Voir {log}" >> "{log}"
        exit 1
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

