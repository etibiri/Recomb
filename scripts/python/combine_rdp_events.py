#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Combine RDP4 validated recombination events from two exports (dnaA & dnaB)
into a single harmonized table and plot a combined breakpoints histogram.

Inputs:  RDP4 CSV-like exports (e.g., "dnaA_masked.csv", "dnaB_masked.csv")
Outputs: combined CSV and a PNG figure

Validation rule (default): event is 'validated' if supported by >= MIN_METHODS
where support is: (i) numeric p-value < 0.05, or (ii) any non-empty token
different from 'NS' (to accommodate RDP exports that mark support with 'Yes', '*', etc.)

Dependencies: Python 3.8+, pandas, matplotlib
"""

import argparse
import csv
import re
from pathlib import Path
from typing import List, Dict, Optional, Tuple

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# -------------------------- Helpers --------------------------

def parse_intlike(s: Optional[str]) -> Optional[int]:
    """Extract first integer from a string such as '1011*' -> 1011; return None if none."""
    if s is None:
        return None
    m = re.search(r"\d+", str(s))
    return int(m.group(0)) if m else None


def sniff_delimiter(text_sample: str) -> Tuple[str, str]:
    """Infer CSV delimiter and quotechar from a text sample."""
    sniffer = csv.Sniffer()
    try:
        dialect = sniffer.sniff(text_sample, delimiters=[",", ";", "\t", "|"])
        delim = dialect.delimiter
        quotechar = getattr(dialect, "quotechar", '"') or '"'
    except Exception:
        delim = ","
        quotechar = '"'
    return delim, quotechar


def parse_rdp_events(csv_path: Path, min_methods: int = 3) -> pd.DataFrame:
    """
    Parse a RDP4 export (CSV-like) and return validated events
    according to min_methods threshold.
    """
    if not csv_path.exists():
        return pd.DataFrame()

    text = csv_path.read_text(encoding="utf-8", errors="ignore")
    sample = text[:200_000]
    delim, quotechar = sniff_delimiter(sample)

    events: List[Dict] = []
    header_methods: List[str] = []
    cols_idx: Dict[str, Optional[int]] = {}
    inside_table = False

    with csv_path.open(encoding="utf-8", errors="ignore", newline="") as f:
        reader = csv.reader(f, delimiter=delim, quotechar=quotechar)
        for row in reader:
            if not row:
                continue

            joined = delim.join((x or "").strip() for x in row)

            # Detect the header row (RDP variants)
            if re.search(r"Recombin(ation|ant)\s+Event\s+Number", joined, flags=re.IGNORECASE):
                # Normalize header cells (collapse spaces)
                hdr = [re.sub(r"\s+", " ", (x or "").strip()) for x in row]

                def find_idx(cands: List[str], default: Optional[int] = None) -> Optional[int]:
                    for i, k in enumerate(hdr):
                        lk = k.lower()
                        if any(c in lk for c in cands):
                            return i
                    return default

                idx_event = find_idx(["recombination event number", "recombinant event number", "event number"], 0)
                idx_begin = find_idx(["begin in alignment", "begin in alignment (nts)", "begin in alignment (bp)", "begin"], 2)
                idx_end   = find_idx(["end in alignment", "end in alignment (nts)", "end in alignment (bp)", "end"], 3)
                idx_recomb= find_idx(["recombinant"], 8)
                idx_minor = find_idx(["minor parent", "minor"], 9)
                idx_major = find_idx(["major parent", "major"], 10)

                # Methods usually start after the last core column
                last_core = max([i for i in [idx_major, idx_minor, idx_recomb, idx_end] if isinstance(i, int)] + [10])
                m_start = last_core + 1
                header_methods = [hdr[i] for i in range(m_start, len(hdr)) if hdr[i].strip()]

                cols_idx = dict(idx_event=idx_event, idx_begin=idx_begin, idx_end=idx_end,
                                idx_recomb=idx_recomb, idx_minor=idx_minor, idx_major=idx_major, m_start=m_start)
                inside_table = True
                continue

            # Data rows
            if inside_table:
                first_cell = row[0].strip() if len(row) > 0 and row[0] is not None else ""
                # RDP tables: data lines typically start with an integer event id in col0
                if re.fullmatch(r"\d+", first_cell):
                    def get(i: Optional[int]) -> str:
                        return row[i].strip() if i is not None and i < len(row) and row[i] is not None else ""

                    ev_id   = get(cols_idx.get("idx_event", 0))
                    begin_s = get(cols_idx.get("idx_begin", 2))
                    end_s   = get(cols_idx.get("idx_end", 3))
                    recomb  = get(cols_idx.get("idx_recomb", 8))
                    minor   = get(cols_idx.get("idx_minor", 9))
                    major   = get(cols_idx.get("idx_major", 10))

                    # Methods support: numeric p<0.05 or non-empty token != 'NS'
                    sig_methods: List[str] = []
                    for j, mname in enumerate(header_methods):
                        col = cols_idx["m_start"] + j
                        val = get(col)
                        if not val:
                            continue
                        v = val.strip().lower()
                        if v == "ns":
                            continue
                        supported = False
                        try:
                            pv = float(v.replace(",", "."))
                            supported = (pv < 0.05)
                        except Exception:
                            # e.g., "Yes", "*", "Sig"
                            supported = True
                        if supported:
                            sig_methods.append(mname)

                    events.append({
                        "event_id": parse_intlike(ev_id),
                        "begin": parse_intlike(begin_s),
                        "end": parse_intlike(end_s),
                        "recombinant": recomb,
                        "major_parent": major,
                        "minor_parent": minor,
                        "n_sig_methods": len(sig_methods),
                        "sig_methods": ", ".join(sig_methods),
                    })

    df = pd.DataFrame(events)

    if df.empty:
        return df

    # Validated: >= min_methods
    return df[df["n_sig_methods"] >= int(min_methods)].copy()


def ensure_begin_end(df: pd.DataFrame) -> pd.DataFrame:
    """If begin/end missing but a breakpoints text column exists, try to derive begin/end."""
    df = df.copy()
    if "begin" not in df.columns or "end" not in df.columns:
        if "breakpoints" in df.columns:
            begins, ends = [], []
            for s in df["breakpoints"].astype(str).tolist():
                nums = re.findall(r"\d+", s)
                if len(nums) >= 2:
                    begins.append(int(nums[0])); ends.append(int(nums[1]))
                elif len(nums) == 1:
                    begins.append(int(nums[0])); ends.append(np.nan)
                else:
                    begins.append(np.nan); ends.append(np.nan)
            df["begin"] = begins; df["end"] = ends
    for c in ("begin", "end"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def harmonize(df: pd.DataFrame, segment_name: str) -> pd.DataFrame:
    """Keep a fixed schema and add 'segment' column."""
    if df.empty:
        return df
    df = ensure_begin_end(df)
    needed = ["event_id", "recombinant", "major_parent", "minor_parent",
              "begin", "end", "n_sig_methods", "sig_methods"]
    for c in needed:
        if c not in df.columns:
            df[c] = np.nan
    out = df[needed].copy()
    out.insert(1, "segment", segment_name)
    return out


def make_histogram(ptsA: List[int], ptsB: List[int], out_png: Path) -> None:
    """Plot combined breakpoints histogram for dnaA and dnaB."""
    fig = plt.figure(figsize=(9, 3.2))
    if ptsA or ptsB:
        all_pts = (ptsA + ptsB) if (ptsA and ptsB) else (ptsA or ptsB or [])
        if all_pts:
            rng = (min(all_pts), max(all_pts))
            bins = 50
            if ptsA:
                plt.hist(ptsA, bins=bins, range=rng, alpha=0.5, label="dnaA")
            if ptsB:
                plt.hist(ptsB, bins=bins, range=rng, alpha=0.5, label="dnaB")
            plt.xlabel("Genomic position (nt)")
            plt.ylabel("Breakpoint frequency (start + end)")
            plt.title("Combined distribution of recombination breakpoints: DNA-A + DNA-B")
            if ptsA and ptsB:
                plt.legend()
        else:
            plt.text(0.5, 0.5, "No breakpoints detected", ha="center", va="center")
            plt.axis("off")
    else:
        plt.text(0.5, 0.5, "No breakpoints detected", ha="center", va="center")
        plt.axis("off")
    plt.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=180)
    plt.close(fig)


# -------------------------- CLI --------------------------

def main():
    ap = argparse.ArgumentParser(description="Combiner des événements RDP4 validés (dnaA+dnaB) en une table harmonisée et tracer les cassures.")
    ap.add_argument("--dnaA", type=Path, required=True, help="Export RDP4 pour dnaA (CSV)")
    ap.add_argument("--dnaB", type=Path, required=True, help="Export RDP4 pour dnaB (CSV)")
    ap.add_argument("--out-csv", type=Path, default=Path("rdp4_validated_events_dnaA_dnaB.csv"),
                    help="Chemin de sortie CSV (table harmonisée)")
    ap.add_argument("--out-png", type=Path, default=Path("recomb_breakpoints_dnaA_dnaB.png"),
                    help="Chemin de sortie PNG (figure combinée)")
    ap.add_argument("--min-methods", type=int, default=3, help="Nombre minimal de méthodes pour valider un événement (défaut: 3)")
    args = ap.parse_args()

    # Parse each segment with identical rule
    A = parse_rdp_events(args.dnaA, min_methods=args.min_methods)
    A["segment"] = "dnaA"
    B = parse_rdp_events(args.dnaB, min_methods=args.min_methods)
    B["segment"] = "dnaB"

    # Harmonize and combine
    A_h = harmonize(A, "dnaA") if not A.empty else pd.DataFrame(columns=["event_id","segment","recombinant","major_parent","minor_parent","begin","end","n_sig_methods","sig_methods"])
    B_h = harmonize(B, "dnaB") if not B.empty else pd.DataFrame(columns=A_h.columns)

    combined = pd.concat([A_h, B_h], ignore_index=True)
    if not combined.empty:
        combined = combined.sort_values(["segment", "event_id", "begin", "end"], na_position="last")

    # Save CSV
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(args.out_csv, index=False)

    # Prepare points and plot
    def collect_pts(df, seg):
        m = df[df["segment"] == seg]
        pts = []
        if "begin" in m.columns:
            pts += [int(x) for x in m["begin"].dropna().astype(int).tolist()]
        if "end" in m.columns:
            pts += [int(x) for x in m["end"].dropna().astype(int).tolist()]
        return pts

    ptsA = collect_pts(combined, "dnaA")
    ptsB = collect_pts(combined, "dnaB")
    make_histogram(ptsA, ptsB, args.out_png)

    # Console summary
    print(f"[OK] Événements validés (dnaA + dnaB) : {combined.shape[0]}")
    print(f"[OK] CSV  → {args.out_csv}")
    print(f"[OK] PNG  → {args.out_png}")


if __name__ == "__main__":
    main()
