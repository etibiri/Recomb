#!/usr/bin/env python3
"""Parse HyPhy GARD JSON to a simple TSV: comp, breakpoint_start, breakpoint_end, support_metric.
Handles missing keys gracefully."""
import json, sys
import click

@click.command()
@click.option('--in', 'inf', required=True)
@click.option('--out', 'outf', required=True)
def main(inf, outf):
    rows = []
    try:
        data = json.load(open(inf))
        bps = data.get('breakpoint-data', []) or data.get('breakpoints', [])
        for bp in bps:
            start = bp.get('position', [None, None])[0]
            end = bp.get('position', [None, None])[1]
            aic = bp.get('AICc-improvement', None)
            rows.append((start, end, aic))
    except Exception as e:
        print(f"[WARN] GARD parse failed: {e}")
    with open(outf, 'w') as w:
        w.write("start\tend\tAICc_impr\n")
        for s,e,a in rows:
            w.write(f"{s}\t{e}\t{a}\n")
    print(f"[INFO] Parsed {len(rows)} breakpoints → {outf}")

if __name__ == '__main__':
    main()
