"""Rank-1 decoys vs entrapment hits: share with a small rt_abs_error (below the median of confident targets)."""
import sys, statistics
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from decoy_vs_entrap import load, category  # noqa
for spec in sys.argv[1:]:
    dataset, arm = spec.split(":")
    peps = load(f"results/{dataset}/{arm}/search.idXML.gz")
    rows = []
    for pi in peps:
        hs = pi.getHits()
        if hs: rows.append((hs[0].getScore(), category(hs[0], dataset), float(hs[0].getMetaValue("rt_abs_error"))))
    rows.sort(key=lambda r: -r[0])
    good = statistics.median(r[2] for r in rows[:len(rows) // 4] if r[1] == "target")
    out = []
    for cat in ("decoy", "entrapment"):
        sel = [r for r in rows if r[1] == cat]
        top = sorted(sel, key=lambda r: -r[0])[: len(sel) // 3]  # best-scoring third of this class
        out.append(f"{cat} n={len(sel)} small-RT-error share all {sum(r[2] < good for r in sel) / len(sel):.3f}"
                   f" top-third {sum(r[2] < good for r in top) / len(top):.3f}")
    print(dataset, arm, f"threshold {good:.0f} s |", " | ".join(out))
