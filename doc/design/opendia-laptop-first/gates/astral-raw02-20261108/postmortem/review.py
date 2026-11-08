"""Extraction-quality review from selection.xipm (spent raw 02, post-mortem only).
Per precursor, around its CNN winner apex (+-4 s): fragments with signal, fragment co-elution,
library-intensity correlation, fragment and MS1 ppm error, MS1 monoisotope and isotope evidence."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyarrow.parquet as pq, sqlite3
from numpress import decode_row
t = pq.read_table("selection.xipm").to_pandas()
sel = pd.read_csv("xipm_selection.tsv", sep="\t")
aw = pd.read_csv("../astral-entrapment-raw02-20261108/all_windows.tsv", sep="\t", usecols=["sample_id", "rt_apex"])
cf = pd.read_csv("../astral-entrapment-raw02-20261108/crossfit/crossfit.tsv", sep="\t", usecols=["sample_id", "precursor_id", "precursor_best"])
apex = cf[cf.precursor_best == 1].merge(aw, on="sample_id").set_index("precursor_id").rt_apex
con = sqlite3.connect("../astral-entrapment-20261107/lib/library.pqp")
lib = pd.read_sql("select t.ID as tid, t.LIBRARY_INTENSITY as li from TRANSITION t", con).set_index("tid").li
rows = []
for pid, g in t.groupby("PRECURSOR_ID"):
    a = apex.get(pid, np.nan)
    if not np.isfinite(a): continue
    frag = g[g.MS_LEVEL == 2]; ms1 = g[g.MS_LEVEL == 1]
    traces = {}; apex_int = {}; ppm = []
    for _, r in frag.iterrows():
        mz, rt, im, it = decode_row(r)
        m = np.abs(rt - a) <= 4.0
        if m.any():
            s = pd.Series(it[m]).groupby(np.round(rt[m], 3)).sum()
            traces[int(r.TRANSITION_ID)] = s
            apex_int[int(r.TRANSITION_ID)] = s.max()
            w = it[m]; ppm.append(np.average((mz[m] - r.TARGET_MZ) / r.TARGET_MZ * 1e6, weights=w))
    n_frag = len(frag); n_sig = len(traces)
    # co-elution: mean pairwise Pearson of fragment traces on the shared RT grid
    co = np.nan
    if n_sig >= 2:
        df = pd.DataFrame(traces).fillna(0.0)
        if len(df) >= 3:
            c = df.corr().values; iu = np.triu_indices_from(c, 1); co = np.nanmean(c[iu])
    libc = np.nan
    if n_sig >= 3:
        ids = list(apex_int); x = np.array([apex_int[i] for i in ids]); y = lib.reindex(ids).values
        if np.isfinite(y).all(): libc = np.corrcoef(np.sqrt(x), np.sqrt(y))[0, 1]
    mono = np.nan; iso_ok = np.nan; ms1_ppm = np.nan
    ms1d = {}
    for _, r in ms1.iterrows():
        mz, rt, im, it = decode_row(r); m = np.abs(rt - a) <= 4.0
        k = int(str(r.ANNOTATION).split("_i")[-1]) if "_i" in str(r.ANNOTATION) else 0
        ms1d[k] = it[m].sum() if m.any() else 0.0
        if k == 0 and m.any(): ms1_ppm = np.average((mz[m] - r.TARGET_MZ) / r.TARGET_MZ * 1e6, weights=it[m])
    if ms1d: mono = ms1d.get(0, 0.0); iso_ok = float(ms1d.get(0, 0) > 0 and ms1d.get(1, 0) > 0)
    rows.append(dict(precursor_id=pid, n_frag=n_frag, n_sig=n_sig, frac_sig=n_sig / max(n_frag, 1), coelution=co,
                     lib_corr=libc, frag_ppm=np.median(ppm) if ppm else np.nan, ms1_mono=mono, ms1_iso12=iso_ok, ms1_ppm=ms1_ppm))
q = pd.DataFrame(rows).merge(sel, on="precursor_id")
q.to_csv("review_per_precursor.tsv", sep="\t", index=False)
order = ["accepted_human", "accepted_entrap", "pair_of_accepted_entrap", "high_decoy", "missed_human_z6_10", "entrap_z6_10",
         "pair_of_entrap_z6_10", "pair_of_missed_human_z6_10", "human_low", "pair_of_accepted_human"]
agg = q.groupby("group").agg(n=("precursor_id", "size"), frac_sig=("frac_sig", "median"), coelution=("coelution", "median"),
                             lib_corr=("lib_corr", "median"), abs_frag_ppm=("frag_ppm", lambda s: s.abs().median()),
                             ms1_present=("ms1_mono", lambda s: (s > 0).mean()), ms1_iso12=("ms1_iso12", "mean"),
                             abs_ms1_ppm=("ms1_ppm", lambda s: s.abs().median()))
print(agg.reindex([o for o in order if o in agg.index]).round(3).to_string())
print("\naccepted entrapment hits:")
print(q[q.group == "accepted_entrap"][["precursor_id", "key", "n_sig", "n_frag", "coelution", "lib_corr", "frag_ppm", "ms1_mono", "ms1_iso12", "ms1_ppm"]].round(3).to_string(index=False))
