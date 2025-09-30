# 1. Recombination Analysis (CMD DNA-A demo)

A reproducible Snakemake workflow to **detect and synthesize recombination signal** in cassava mosaic Begomoviruses (CMD), demonstrated on **DNA-A**. It wires together only the **existing environments** and **provided scripts** to curate sequences, align, mask, detect breakpoints with HyPhy/GARD, consolidate events, assign putative parents (BLAST vs optional references), and summarize with both **phylogenetic trees (IQ-TREE)** and a **NeighborNet** network.

---

## 2. Scientific context

Cassava mosaic disease (CMD) is driven by Begomoviruses (family *Geminiviridae*), where recombination plays a key role in the emergence of novel variants. This workflow focuses on **DNA-A** as a demonstration and provides a clean orchestration for:

- quality-controlled **sequence curation**  
- **MAFFT** multiple sequence alignment + masking of gap-rich columns  
- **HyPhy GARD** recombination breakpoint inference  
- event **consolidation/consensus**  
- segment-wise **putative parentage** via BLAST (if references are provided)  
- **ML phylogeny** (IQ-TREE) and **NeighborNet** visualization (R)

The pipeline is deterministic (fixed seeds), offline-ready, and HPC-friendly.

---

## 3. Methods (high level)

1. **Curation**: `scripts/curate_sequences.py` cleans `data/dnaA.fasta` and emits a QC table.  
2. **Alignment**: **MAFFT** aligns curated sequences (`mafft_opts` in `config.yaml`).  
3. **Masking**: `scripts/mask_alignment.py` removes poorly informative / gap-rich columns.  
4. **Recombination**: `scripts/hyphy_batch_runner.py` runs **GARD**; `scripts/gard_to_breaks.py` standardizes breakpoints.  
5. **Consensus**: `scripts/consensus_breakpoints.py` synthesizes a robust set of breakpoints (GARD-only by default, expandable).  
6. **Parentage**: Non-recombinant segments (alignment coordinates) are BLASTed vs `data/refs/` (if present) for putative major/minor parents.  
7. **Phylogeny & Network**: **IQ-TREE** (ML tree) and `scripts/neighborNet.R` (network) on the masked alignment.  
8. **Summary**: Consolidated table links consensus events to parentage calls.

---

## Repository layout

