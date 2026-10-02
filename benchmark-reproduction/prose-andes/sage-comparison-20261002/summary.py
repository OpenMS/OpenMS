"""Group means of Sage and ProSE arms on the 20-file suite, with the share of the develop-to-Sage gap each arm closes."""
import json, statistics, sys
from pathlib import Path
from compare import SAGE, PROSE, MANIFEST, group_of
ARMS = sys.argv[1:] or ["base", "pr10335", "ds1", "fz2", "ds1_fz2", "ds1_mm4", "ds1_mm3"]
def mean(folder):
    p = folder / "summary.json"
    return json.loads(p.read_text())["mean_psms"] if p.exists() else None
groups = list(dict.fromkeys(group_of(r["id"]) for r in MANIFEST))
print(f"{'group':22s}{'Sage':>9s}" + "".join(f"{a:>17s}" for a in ARMS))
for g in groups:
    ds = [r["id"] for r in MANIFEST if group_of(r["id"]) == g]
    sage = statistics.mean(mean(SAGE / d) for d in ds)
    base = statistics.mean(mean(PROSE / d / "base") for d in ds)
    cells = []
    for a in ARMS:
        v = [mean(PROSE / d / a) for d in ds]
        if None in v: cells.append(f"{'—':>17s}"); continue
        m = statistics.mean(v)
        gap = f"{100*(m-base)/(sage-base):+.0f}%gap" if a != "base" and abs(sage - base) > 1 else ""
        cells.append(f"{m:8.1f} {gap:>8s}")
    print(f"{g:22s}{sage:9.1f}" + "".join(cells))
