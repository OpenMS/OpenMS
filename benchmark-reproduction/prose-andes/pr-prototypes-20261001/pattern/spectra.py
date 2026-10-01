"""Raw MS2 spectrum properties per file (the 8000 benchmark spectra, before ProSE preprocessing).

peaks: median peak count; per100: median peaks per 100 Da of occupied m/z range; zero: share of zero-intensity
peaks; top20_peaks / top20_tic: share of peaks and of TIC that survive a top-20-per-100-Da window filter
(ProSE peaks:window_top 20, sliding windows approximated by fixed 100 Da bins); dyn: median log10(max / median
intensity); low_frac: share of TIC below m/z 200.
"""
import csv, json, math, statistics, sys
from pathlib import Path
import pyopenms as oms

ROOT = Path(__file__).resolve().parent.parent
BENCH = ROOT.parent / "prose-andes-reproduction/broad_benchmark"
MANIFEST = json.loads((BENCH / "manifest.json").read_text())
rows = []
for r in MANIFEST:
    d = r["id"]
    exp = oms.MSExperiment(); oms.MzMLFile().load(str(BENCH / "data" / f"{d}.mzML"), exp)
    peaks, per100, zero, keep_p, keep_t, dyn, low = [], [], [], [], [], [], []
    for s in exp:
        if s.getMSLevel() != 2: continue
        mz, it = s.get_peaks()
        n = len(mz)
        if n < 5: continue
        peaks.append(n)
        span = max(mz) - min(mz)
        per100.append(100 * n / max(span, 100))
        zero.append(sum(1 for x in it if x <= 0) / n)
        bins = {}
        for m, i in zip(mz, it): bins.setdefault(int(m // 100), []).append(i)
        kept = [sorted(v, reverse=True)[:20] for v in bins.values()]
        tic = sum(it)
        keep_p.append(sum(map(len, kept)) / n)
        keep_t.append(sum(map(sum, kept)) / tic if tic > 0 else 1)
        pos = sorted(x for x in it if x > 0)
        if pos: dyn.append(math.log10(pos[-1] / pos[len(pos) // 2]))
        low.append(sum(i for m, i in zip(mz, it) if m < 200) / tic if tic > 0 else 0)
    inst = r.get("instrument") or r.get("group") or ""
    rows.append(dict(dataset=d, spectra=len(peaks), peaks=int(statistics.median(peaks)),
                     per100=round(statistics.median(per100), 1), zero_pct=round(100 * statistics.mean(zero), 2),
                     top20_peaks_pct=round(100 * statistics.mean(keep_p), 1),
                     top20_tic_pct=round(100 * statistics.mean(keep_t), 1),
                     dyn_log10=round(statistics.median(dyn), 2), tic_below200_pct=round(100 * statistics.mean(low), 1)))
    print("\t".join(str(v) for v in rows[-1].values()), flush=True)
with open(ROOT / "pattern" / "spectra.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
