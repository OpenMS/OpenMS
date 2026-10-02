"""Break the Sage-only / ProSE-only / shared spectra of each group down by precursor charge, peptide length and
Sage's matched-peak count, from the frozen spectrum provenance and Sage's PIN."""
import collections, csv, json, sys
from pathlib import Path
from compare import accepted, candidates, classify, SAGE, PROSE, MANIFEST, B, group_of
arm = sys.argv[1] if len(sys.argv) > 1 else "base"
groups = sys.argv[2:]
stat = collections.defaultdict(lambda: collections.defaultdict(collections.Counter))
for r in MANIFEST:
    d = r["id"]; g = group_of(d)
    if groups and g not in groups: continue
    if not (SAGE / d / "summary.json").exists() or not (PROSE / d / arm / "summary.json").exists(): continue
    prov = {s["scan"]: s for s in json.loads((B / f"prose-andes-reproduction/broad_benchmark/provenance/{d}.json").read_text())["scans"]}
    sa, pa = accepted(SAGE / d), accepted(PROSE / d / arm)
    pc = candidates(PROSE / d / arm / "native.pin", "score")
    for scan in set(sa) | set(pa):
        if scan in sa and pa.get(scan) == sa[scan]: cls = "both"
        elif scan in sa: cls = "sage_only:" + classify(scan, sa[scan], pa, pc)
        else: cls = "prose_only"
        z = prov[scan]["charge"]; seq = sa.get(scan) or pa[scan]
        stat[g]["charge"][(cls, min(z, 4))] += 1
        stat[g]["length"][(cls, "<=12" if len(seq) <= 12 else "13-20" if len(seq) <= 20 else ">20")] += 1
for g, s in stat.items():
    print(f"\n== {g} ({arm})")
    for dim in ["charge", "length"]:
        cats = sorted({k[1] for k in s[dim]}, key=str)
        clss = sorted({k[0] for k in s[dim]})
        print(f"  {dim:8s}" + "".join(str(c).rjust(9) for c in cats) + "   (share of row)")
        for cl in clss:
            n = sum(s[dim][(cl, c)] for c in cats)
            print(f"  {cl[:34]:34s}" if False else f"  {cl:34s}"[:36] + "".join(f"{100*s[dim][(cl, c)]/n:8.0f}%" for c in cats) + f"  n={n}")
