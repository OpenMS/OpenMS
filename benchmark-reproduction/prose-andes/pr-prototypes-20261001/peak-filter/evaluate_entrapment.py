"""Entrapment FDP for arms searched against the shuffled paired-mass entrapment databases (r = 1).

A PSM is an entrapment hit when every protein it maps to is ENTRAP_*. FDP estimates (Wen et al. 2025):
lower bound N_E / N, combined 2 N_E / N. Per seed at Percolator q <= 0.01, averaged over seeds; plus the
PSMs an arm gains over base_E (by scan, pooled over seeds) and the combined FDP among them."""
import csv, json, statistics, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
ARMS = sys.argv[1:] or ["base_E", "wt100_E", "pf_E"]
FILES = ["astral_A2", "astral_B1", "astral_B3", "tims_plasma_30min", "tims_plasma_50min"]


def accepted(path):
    out = {}
    for line in open(path).read().splitlines()[1:]:
        f = line.split("\t")
        if float(f[2]) <= 0.01:
            scan = f[0].rsplit("_", 2)[-2]
            out[scan] = (f[4], all(p.startswith(PREFIX) for p in f[5:] if p))
    return out


rows = []
for d in FILES:
  for arm in ARMS:
    suffix = arm.rsplit("_", 1)[1]
    PREFIX = {"E": "ENTRAP_", "F": "FOREIGN_"}[suffix]
    base = {s: accepted(ROOT / "results" / d / f"base_{suffix}" / f"s{s}" / "target.tsv") for s in (1, 42, 137)
            if (ROOT / "results" / d / f"base_{suffix}" / f"s{s}" / "target.tsv").exists()}
    if True:
        folder = ROOT / "results" / d / arm
        if not (folder / "summary.json").exists(): continue
        n, ne, gained, gained_e = [], [], 0, 0
        for s in (1, 42, 137):
            acc = accepted(folder / f"s{s}" / "target.tsv")
            n.append(len(acc)); ne.append(sum(e for _, e in acc.values()))
            if not arm.startswith("base_") and s in base:
                g = [e for scan, (_, e) in acc.items() if scan not in base[s]]
                gained += len(g); gained_e += sum(g)
        row = dict(dataset=d, arm=arm, psms=round(statistics.mean(n), 1), entrapment=round(statistics.mean(ne), 1),
                   fdp_lower_pct=round(100 * sum(ne) / sum(n), 2), fdp_combined_pct=round(200 * sum(ne) / sum(n), 2),
                   gained_psms=gained, gained_entrapment=gained_e,
                   gained_fdp_combined_pct=round(200 * gained_e / gained, 1) if gained else "")
        rows.append(row); print("\t".join(map(str, row.values())))
if rows:
    with open(ROOT / "entrapment" / "entrapment.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
