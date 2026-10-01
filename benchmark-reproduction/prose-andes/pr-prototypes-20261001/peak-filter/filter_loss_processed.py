"""Per-spectrum share of intensity that ProSE's local window filter (jump_full, 100 Da, top K) removes, after ProSE's
own preprocessing up to that step: drop zero intensities, normalize to max, deisotope (high-res only, same parameters
as ProSEAlgorithm::preprocessSpectra_). Reports the mean and the share of spectra above 10/20/30% for K = 20."""
import csv, json, statistics, sys
from pathlib import Path
import pyopenms as oms
ROOT = Path(__file__).resolve().parent.parent
BENCH = ROOT.parent / "prose-andes-reproduction/broad_benchmark"


def local_top(mz, it, k):
    keep, b = [], 0
    while b < len(mz):
        e = b + 1
        while e < len(mz) and mz[e] - mz[b] < 100.0: e += 1
        keep += sorted(range(b, e), key=lambda i: (-it[i], i))[:k]
        b = e
    return keep


rows = []
for r in json.loads((BENCH / "manifest.json").read_text()):
    ppm = r["fragment_unit"] == "ppm"; tol = r["fragment_tolerance"]
    exp = oms.MSExperiment(); oms.MzMLFile().load(str(BENCH / "data" / f"{r['id']}.mzML"), exp)
    lost, n_after = [], []
    for s in exp:
        if s.getMSLevel() != 2: continue
        mz, it = s.get_peaks()
        pk = [(m, i) for m, i in zip(mz, it) if i > 0]
        if len(pk) < 5: continue
        top = max(i for _, i in pk)
        t = oms.MSSpectrum(); t.set_peaks(([m for m, _ in pk], [i / top for _, i in pk])); t.sortByPosition()
        if ppm:
            oms.Deisotoper.deisotopeAndSingleCharge(t, tol, True, 1, 3, False, 3, 10, True)
            t.sortByPosition()
        mz, it = t.get_peaks(); tic = float(sum(it))
        if tic <= 0: continue
        kept = local_top(list(mz), list(it), 20)
        lost.append(1 - sum(it[i] for i in kept) / tic); n_after.append(len(mz))
    row = dict(dataset=r["id"], peaks_after_deisotope=int(statistics.median(n_after)), lost_mean_pct=round(100 * statistics.mean(lost), 1),
               **{f"spectra_lost_gt{x}_pct": round(100 * sum(v > x / 100 for v in lost) / len(lost), 1) for x in (10, 20, 30)})
    rows.append(row); print("\t".join(map(str, row.values())), flush=True)
with open(ROOT / "pattern" / "filter_loss_processed.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
