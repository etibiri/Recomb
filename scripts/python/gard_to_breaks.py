#!/usr/bin/env python3
"""
Parse HyPhy GARD output (GARD.json) and extract breakpoints to a compact TSV.

Output format:
  #SUMMARY\talignment\tgard_json\tn_breakpoints\tbreaks\tbest_aicc
  #DETAIL \talignment\tbreak_index\tposition

Usage:
  python scripts/python/gard_to_breaks.py --in results/recombination/dnaA.gard.json \
                                          --out results/recombination/dnaA.gard_breakpoints.tsv
"""
import json
from pathlib import Path
import sys
import click

def _load_json(fp: Path):
    try:
        with open(fp, "r") as fh:
            return json.load(fh)
    except FileNotFoundError:
        raise SystemExit(f"[ERROR] File not found: {fp}")
    except json.JSONDecodeError as e:
        raise SystemExit(f"[ERROR] Invalid JSON in {fp}: {e}")

def _parse_breakpoints(data: dict):
    """
    Try multiple known layouts used by HyPhy GARD.
    Returns (breaks_sorted_list, aicc_or_None)
    """
    breaks = []
    aicc = None

    # Common layout (hyphy-analyses/GARD)
    # e.g., {"breakpointData":{"breakpoints":[...]},"bestModelAICc":...}
    if isinstance(data, dict) and "breakpointData" in data:
        bpdata = data["breakpointData"]
        if isinstance(bpdata, dict) and "breakpoints" in bpdata:
            if isinstance(bpdata["breakpoints"], list):
                try:
                    breaks = [int(x) for x in bpdata["breakpoints"]]
                except Exception:
                    # fallback: coerce numbers from strings
                    cleaned = []
                    for x in bpdata["breakpoints"]:
                        try:
                            cleaned.append(int(str(x).strip()))
                        except Exception:
                            pass
                    breaks = cleaned
        aicc = data.get("bestModelAICc", None)

    # Alternate: partitions list with "from"/"to" boundaries; infer cuts at ends of partitions
    # e.g., {"data":{"partitions":[{"from":1,"to":312}, {"from":313,"to":1021}, ...]}}
    if not breaks and "data" in data and isinstance(data["data"], dict):
        parts = data["data"].get("partitions")
        if isinstance(parts, list) and len(parts) > 1:
            candidates = []
            for i in range(len(parts) - 1):
                left = parts[i]
                if isinstance(left, dict) and "to" in left:
                    try:
                        candidates.append(int(left["to"]))
                    except Exception:
                        pass
            breaks = candidates

    # Final cleanup
    breaks = sorted({b for b in breaks if isinstance(b, int) and b > 0})
    return breaks, aicc

@click.command()
@click.option("--in", "in_json", required=True, help="Path to HyPhy GARD.json")
@click.option("--out", "out_tsv", required=True, help="Path to TSV to write")
def main(in_json, out_tsv):
    in_path = Path(in_json)
    out_path = Path(out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    data = _load_json(in_path)
    breaks, aicc = _parse_breakpoints(data)

    # alignment name guessed from JSON filename stem (e.g., dnaA)
    # If your JSON is named "dnaA.gard.json", stem is "dnaA.gard"; take first token before ".gard"
    stem = in_path.stem  # e.g., "dnaA.gard" or "GARD"
    # Try to keep a clean alignment label
    if ".gard" in stem.lower():
        aln_name = stem.split(".")[0]
    else:
        aln_name = stem

    # Write TSV
    with open(out_path, "w") as out:
        out.write("#SUMMARY\talignment\tgard_json\tn_breakpoints\tbreaks\tbest_aicc\n")
        out.write(
            "#SUMMARY\t{aln}\t{json}\t{n}\t{lst}\t{aicc}\n".format(
                aln=aln_name,
                json=str(in_path),
                n=len(breaks),
                lst=";".join(map(str, breaks)) if breaks else "",
                aicc=aicc if aicc is not None else ""
            )
        )
        out.write("#DETAIL\talignment\tbreak_index\tposition\n")
        for i, bp in enumerate(breaks, start=1):
            out.write(f"#DETAIL\t{aln_name}\t{i}\t{bp}\n")

    print(f"[OK] Parsed {len(breaks)} breakpoint(s) → {out_path}")

if __name__ == "__main__":
    main()
