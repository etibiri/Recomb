#!/usr/bin/env python3
"""
Batch runner for HyPhy selection methods (SLAC/FEL/FUBAR/MEME) on non-recombinant fragments.

Expected fragment layout:
  <fragments_root>/
    dnaA/fragment_001.fasta
    dnaA/fragment_002.fasta
    dnaB/fragment_001.fasta
    ...

Tree handling:
  - If --tree-dir is provided, the script will look for a .nwk file with the same base name as the fragment:
        <tree_dir>/dnaA/fragment_001.nwk
    Fallback within the fragment's output folder if not found.
  - If --build-tree is set, IQ-TREE2 will be used to infer a maximum-likelihood tree for each fragment.
    (Model defaults are for nucleotide data; adjust with --iqtree-model for aa/codon if needed.)
  - If neither a tree is found nor --build-tree is set, methods requiring a tree will be skipped for that fragment.

Usage example:
  python scripts/python/hyphy_batch_runner.py \
    --fragments results/recombination/fragments \
    --outdir results/hyphy \
    --tree-dir results/trees \
    --build-tree \
    --methods SLAC FEL FUBAR MEME \
    --cpu 4 \
    --genetic-code 1 \
    --iqtree-model GTR+G
"""
import click
import glob
import os
import shutil
import subprocess
from pathlib import Path

# ---------- helpers ----------
def run(cmd, cwd=None, allow_fail=False):
    p = subprocess.run(cmd, cwd=cwd, text=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if p.returncode != 0 and not allow_fail:
        raise RuntimeError(
            f"Command failed ({' '.join(cmd)}):\n"
            f"--- STDOUT ---\n{p.stdout}\n"
            f"--- STDERR ---\n{p.stderr}\n"
        )
    return p

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def guess_tree_for_fragment(frag_path: Path, tree_dir: Path | None, outp: Path) -> Path | None:
    """
    Try a few sensible locations:
      1) <tree_dir>/<dnaX>/<fragment_xxx>.nwk
      2) <tree_dir>/<fragment_xxx>.nwk
      3) <outp>/<fragment_xxx>.treefile (if IQ-TREE was run here)
    """
    base = frag_path.stem  # fragment_001
    sub  = frag_path.parent.name  # dnaA / dnaB / etc.

    candidates = []
    if tree_dir:
        candidates += [
            tree_dir / sub / f"{base}.nwk",
            tree_dir / f"{base}.nwk",
        ]
    candidates += [
        outp / f"{base}.treefile",   # IQ-TREE default
        outp / f"{base}.nwk"
    ]
    for c in candidates:
        if c.exists() and c.stat().st_size > 0:
            return c
    return None

def build_tree_with_iqtree(frag_path: Path, outp: Path, cpu: int, iqtree_model: str):
    """
    Build a quick ML tree with IQ-TREE2. Output files go to outp.
    """
    ensure_dir(outp)
    base = frag_path.stem
    # IQ-TREE2 typical call for nucleotide alignment:
    #   iqtree2 -s <aln> -m GTR+G -nt <cpu> -pre <prefix>
    prefix = outp / base
    cmd = [
        "iqtree2",
        "-s", str(frag_path),
        "-m", iqtree_model,
        "-nt", str(cpu),
        "-pre", str(prefix),
        "-quiet"
    ]
    run(cmd)
    # Return tree path if created
    treefile = outp / f"{base}.treefile"
    return treefile if treefile.exists() else None

def hyphy_call(method: str, aln: Path, tree: Path, outp: Path, genetic_code: int, cpu: int):
    """
    Run a HyPhy method on a codon/nucleotide alignment with a provided tree.
    HyPhy CLI (modern analyses):
      hyphy <method> --alignment <aln> --tree <tree> --output <out_dir> [--code <gc>] [--threads N]
    """
    ensure_dir(outp)
    cmd = [
        "hyphy", method.lower(),            # slac / fel / fubar / meme
        "--alignment", str(aln),
        "--tree", str(tree),
        "--output", str(outp),
        "--code", str(genetic_code),
        "--threads", str(cpu)
    ]
    return run(cmd, allow_fail=False)

# ---------- CLI ----------
@click.command()
@click.option("--fragments", required=True,
              help="Racine contenant les fichiers dna*/fragment_*.fasta")
@click.option("--outdir", required=True,
              help="Répertoire racine des résultats HyPhy.")
@click.option("--tree-dir", default=None,
              help="Répertoire racine où chercher des arbres .nwk (même base que les fragments).")
@click.option("--build-tree", is_flag=True,
              help="Si aucun arbre n'est trouvé, construire un arbre avec IQ-TREE2.")
@click.option("--methods", multiple=True, default=["SLAC","FEL","FUBAR","MEME"],
              help="Méthodes HyPhy à exécuter (SLAC, FEL, FUBAR, MEME).")
@click.option("--cpu", default=1, show_default=True, type=int,
              help="Threads CPU pour HyPhy (et IQ-TREE2).")
@click.option("--genetic-code", default=1, show_default=True, type=int,
              help="Code génétique HyPhy (1=Standard).")
@click.option("--iqtree-model", default="GTR+G", show_default=True,
              help="Modèle pour IQ-TREE2 quand --build-tree est actif.")
def main(fragments, outdir, tree_dir, build_tree, methods, cpu, genetic_code, iqtree_model):
    fragments = Path(fragments)
    outdir = Path(outdir)
    tree_dir = Path(tree_dir) if tree_dir else None

    ensure_dir(outdir)

    frags = sorted(glob.glob(str(fragments / "dna*" / "fragment_*.fasta")))
    if not frags:
        raise SystemExit(f"[ERROR] Aucun fragment trouvé sous {fragments}/dna*/fragment_*.fasta")

    allowed = {"SLAC","FEL","FUBAR","MEME"}
    todo = [m.upper() for m in methods if m.upper() in allowed]
    if not todo:
        raise SystemExit("[ERROR] Aucune méthode valide. Choisir parmi: SLAC FEL FUBAR MEME")

    print(f"[INFO] N fragments: {len(frags)} | Méthodes: {', '.join(todo)}")

    for frag in frags:
        frag_path = Path(frag)
        sub = frag_path.parent.name           # dnaA / dnaB
        base = frag_path.stem                 # fragment_001
        outp = outdir / sub / base
        ensure_dir(outp)

        # 1) arbre
        tree = guess_tree_for_fragment(frag_path, tree_dir, outp)
        if not tree and build_tree:
            print(f"[INFO] Pas d'arbre pour {sub}/{base}, construction avec IQ-TREE2…")
            try:
                tree = build_tree_with_iqtree(frag_path, outp, cpu=cpu, iqtree_model=iqtree_model)
            except Exception as e:
                print(f"[WARN] IQ-TREE2 a échoué pour {frag_path}: {e}")
                tree = None

        if not tree or not tree.exists():
            print(f"[WARN] Aucun arbre pour {sub}/{base}. "
                  f"Les méthodes nécessitant un arbre seront sautées.")
            available_tree = False
        else:
            available_tree = True

        # 2) HyPhy methods
        for m in todo:
            m_out = outp / m.lower()
            ensure_dir(m_out)

            # SLAC/FEL/FUBAR/MEME nécessitent un arbre
            if not available_tree:
                with open(m_out / "SKIPPED.txt", "w") as fh:
                    fh.write("Skipped: no tree available for this fragment.\n")
                print(f"[SKIP] {m} sur {sub}/{base} (pas d'arbre)")
                continue

            try:
                print(f"[RUN ] {m:5s} → {sub}/{base}")
                hyphy_call(m, frag_path, tree, m_out, genetic_code, cpu)
                # HyPhy écrit typiquement <Method>.json dans m_out
                print(f"[OK  ] {m:5s} ✓ {sub}/{base}")
            except Exception as e:
                # On n'arrête pas tout le batch : on note l'échec pour ce fragment/méthode.
                with open(m_out / "ERROR.txt", "w") as fh:
                    fh.write(str(e) + "\n")
                print(f"[FAIL] {m:5s} ✗ {sub}/{base} -> voir {m_out}/ERROR.txt")

    print(f"[DONE] HyPhy batch terminé → {outdir}")

if __name__ == "__main__":
    main()

