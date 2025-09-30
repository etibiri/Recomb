#!/usr/bin/env python3
"""Run HyPhy selection methods on each fragment (non-recombinant).
This is a controller that would call hyphy SLAC/FEL/FUBAR/MEME; here we stub commands.
"""
import click, glob, os, subprocess

@click.command()
@click.option('--fragments', required=True)
@click.option('--outdir', required=True)
def main(fragments, outdir):
    os.makedirs(outdir, exist_ok=True)
    for frag in sorted(glob.glob(os.path.join(fragments, 'dna*', 'fragment_*.fasta'))):
        base = os.path.basename(frag).replace('.fasta','')
        outp = os.path.join(outdir, base)
        os.makedirs(outp, exist_ok=True)
        # Example: hyphy slac --alignment frag --tree tree (requires tree); left as TODO
        # subprocess.run(['hyphy','slac','--alignment',frag,'--tree',treefile,'--output',outp+'/slac.json'], check=False)
        open(os.path.join(outp,'README.txt'),'w').write('Run SLAC/FEL/FUBAR/MEME here with proper trees.')
    print(f"[INFO] HyPhy batch stub done → {outdir}")

if __name__ == '__main__':
    main()

