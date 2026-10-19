#!/usr/bin/env python3
"""andes reference for PR #10635: andes 72ee0a5 (no rule) vs 3dded7b (bigbio/andes#102 rule) with --refine, and a
no-refine control, on the Astral files; Percolator seeds 1/42/137, target PSMs at q <= 0.01. Same residual metric as
evaluate.py (precursor error peptide - spectrum at the PSM's isotope offset, err = -dm * z).

usage: evaluate_andes.py <runs/andes> <out tsv>"""
import csv, statistics, sys
from pathlib import Path

runs, out = Path(sys.argv[1]), Path(sys.argv[2])
DEAM = "+0.98402"
lines = ["arm\tdataset\tpsms_s1/s42/s137\tmean\tpass2_mean\tdeam_mean\tdeam_near0_pct\tdeam_mispick_pct\t"
         "deam_mimic_rule_pct\tdropped_as_mimics"]
for arm in ["rule_norefine", "norule", "rule"]:
    for d in sorted((runs / arm).glob("*")):
        rows = {}
        with open(d / "input.pin") as f:
            r = csv.reader(f, delimiter="\t")
            h = next(r)
            ix = {n: i for i, n in enumerate(h)}
            zc = [(int(n[6:]), ix[n]) for n in h if n.startswith("charge") and n[6:].isdigit()]
            for row in r:
                z = next(c for c, i in zc if row[i] == "1")
                rows[(row[ix["SpecId"]], row[ix["Peptide"]])] = (-float(row[ix["dm"]]) * z,
                                                                 int(row[ix["IsRefinement"]]))
        n, p2, dm, near, mis, rule = [], [], [], [], [], []
        for s in ("1", "42", "137"):
            acc = [x for x in csv.DictReader(open(d / f"s{s}" / "target.tsv"), delimiter="\t")
                   if float(x["q-value"]) <= 0.01]
            n.append(len(acc))
            p2.append(sum(rows[(x["PSMId"], x["peptide"])][1] for x in acc))
            de = [rows[(x["PSMId"], x["peptide"])][0] * 1000 for x in acc if DEAM in x["peptide"]]
            k = [x["peptide"].count(DEAM) for x in acc if DEAM in x["peptide"]]
            dm.append(len(de))
            near.append(sum(abs(e) <= 5 for e in de))
            mis.append(sum(-26 <= e <= -12 for e in de))
            gap = 19.3394
            rule.append(sum(abs(e + gap * c) < abs(e) for e, c in zip(de, k)))
        dropped = ""
        log = (d / "andes.log").read_text()
        if "dropped as isotope mimics" in log:
            dropped = log.split("winners; ")[1].split(" dropped")[0]
        tot = sum(dm) or 1
        lines.append(f"{arm}\t{d.name}\t{'/'.join(map(str, n))}\t{statistics.mean(n):.1f}\t{statistics.mean(p2):.1f}\t"
                     f"{statistics.mean(dm):.1f}\t{100 * sum(near) / tot:.1f}\t{100 * sum(mis) / tot:.1f}\t"
                     f"{100 * sum(rule) / tot:.1f}\t{dropped}")
out.write_text("\n".join(lines) + "\n")
print(out.read_text())
