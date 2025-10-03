#!/usr/bin/env python3
import json, sys, math, pathlib

if len(sys.argv) < 3:
    print(f"Usage: {sys.argv[0]} <gard.json> <out.tsv> [alignment_name]", file=sys.stderr)
    sys.exit(2)

jin, jout = sys.argv[1], sys.argv[2]
aln_name = sys.argv[3] if len(sys.argv) > 3 else pathlib.Path(jin).stem

with open(jin, "r") as fh:
    data = json.load(fh)

# 1) Récupérer les breakpoints finaux (nouveau schéma HyPhy : 'improvements')
breaks = []
best_aicc = data.get("bestModelAICc")
if isinstance(data.get("improvements"), dict) and data["improvements"]:
    # prendre l’étape au plus grand index (finale)
    last_key = max((int(k) for k in data["improvements"].keys()))
    step = data["improvements"][str(last_key)]
    bp = step.get("breakpoints", [])
    # flattener (car souvent [[530],[1005],...])
    flat = []
    for x in bp:
        if isinstance(x, list):
            flat.extend(x)
        else:
            flat.append(x)
    breaks = sorted(set(int(v) for v in flat))

# 2) Rétro-compat: si pas trouvé, essayer anciens champs possibles
# (Exemples d’anciens schémas; on ne devine pas, on reste conservateur)
for alt in ("breakpointList", "breakpoints", "bestBreakPoints"):
    if not breaks and isinstance(data.get(alt), list):
        breaks = sorted(set(int(v) for v in data[alt]))

# 3) À défaut, on renonce proprement
if best_aicc is None and "bestModelAICc" in data:
    best_aicc = data["bestModelAICc"]
if not isinstance(best_aicc, (int, float)):
    # fallback : utiliser singleTreeAICc si présent pour éviter 'None' dans la sortie
    best_aicc = data.get("singleTreeAICc", float("nan"))

# 4) Écrire le TSV
with open(jout, "w") as out:
    out.write("#SUMMARY\talignment\tgard_json\tn_breakpoints\tbreaks\tbest_aicc\n")
    out.write(f"#SUMMARY\t{aln_name}\t{jin}\t{len(breaks)}\t{','.join(map(str, breaks))}\t{best_aicc}\n")
    out.write("#DETAIL\talignment\tbreak_index\tposition\n")
    for i, b in enumerate(breaks, start=1):
        out.write(f"#DETAIL\t{aln_name}\t{i}\t{b}\n")
