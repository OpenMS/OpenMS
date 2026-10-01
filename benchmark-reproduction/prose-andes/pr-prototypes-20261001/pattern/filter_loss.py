"""Share of fragment-region (m/z >= 200) peaks and intensity removed by a top-K-per-100-Da filter, K = 20 / 40 / 100,
and the share of MS2 spectra with any 100 Da bin above 20 peaks."""
import csv, json, statistics
from pathlib import Path
import pyopenms as oms
ROOT = Path(__file__).resolve().parent.parent
BENCH = ROOT.parent / "prose-andes-reproduction/broad_benchmark"
rows = []
for r in json.loads((BENCH / "manifest.json").read_text()):
    exp = oms.MSExperiment(); oms.MzMLFile().load(str(BENCH / "data" / f"{r['id']}.mzML"), exp)
    lost = {20: [], 40: [], 100: []}; lostp = {20: [], 40: [], 100: []}; crowded = []
    for s in exp:
        if s.getMSLevel() != 2: continue
        mz, it = s.get_peaks()
        bins = {}
        for m, i in zip(mz, it):
            if m >= 200: bins.setdefault(int(m // 100), []).append(i)
        tic = sum(map(sum, bins.values())); n = sum(map(len, bins.values()))
        if tic <= 0: continue
        crowded.append(any(len(v) > 20 for v in bins.values()))
        for k in lost:
            kept = [sorted(v, reverse=True)[:k] for v in bins.values()]
            lost[k].append(1 - sum(map(sum, kept)) / tic); lostp[k].append(1 - sum(map(len, kept)) / n)
    rows.append(dict(dataset=r["id"], crowded_spectra_pct=round(100 * statistics.mean(crowded), 1),
                     **{f"tic_lost_top{k}_pct": round(100 * statistics.mean(v), 1) for k, v in lost.items()},
                     **{f"peaks_lost_top{k}_pct": round(100 * statistics.mean(v), 1) for k, v in lostp.items()}))
    print("\t".join(str(v) for v in rows[-1].values()), flush=True)
with open(ROOT / "pattern" / "filter_loss.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
