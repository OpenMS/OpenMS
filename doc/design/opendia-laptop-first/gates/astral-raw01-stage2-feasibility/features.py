"""Extraction features per precursor from an .xipm, at the CNN winner apex (+-4 s).
Usage: features.py <xipm> <selection.tsv with precursor_id,sample_id> <all_windows.tsv> <library.pqp> <out.tsv>"""
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyarrow.parquet as pq, sqlite3
from numpress import decode_row
xipm, selp, awp, libp, outp = sys.argv[1:6]
t = pq.read_table(xipm).to_pandas()
sel = pd.read_csv(selp, sep="\t")
aw = pd.read_csv(awp, sep="\t", usecols=["sample_id", "rt_apex"]).set_index("sample_id").rt_apex
apex = sel.set_index("precursor_id").sample_id.map(aw)
lib = pd.read_sql("select ID as tid, LIBRARY_INTENSITY as li from TRANSITION", sqlite3.connect(libp)).set_index("tid").li
rows = []
for pid, g in t.groupby("PRECURSOR_ID"):
    a = apex.get(pid, np.nan)
    if not np.isfinite(a): continue
    traces, apex_int, ppm = {}, {}, []
    frag = g[g.MS_LEVEL == 2]
    for _, r in frag.iterrows():
        mz, rt, im, it = decode_row(r); m = np.abs(rt - a) <= 4.0
        if m.any():
            s = pd.Series(it[m]).groupby(np.round(rt[m], 3)).sum()
            traces[int(r.TRANSITION_ID)] = s; apex_int[int(r.TRANSITION_ID)] = s.max()
            ppm.append(np.average((mz[m] - r.TARGET_MZ) / r.TARGET_MZ * 1e6, weights=it[m]))
    co = libc = np.nan
    if len(traces) >= 2:
        df = pd.DataFrame(traces).fillna(0.0)
        if len(df) >= 3:
            c = df.corr().values; co = np.nanmean(c[np.triu_indices_from(c, 1)])
    if len(traces) >= 3:
        ids = list(apex_int); y = lib.reindex(ids).values
        if np.isfinite(y).all(): libc = np.corrcoef(np.sqrt([apex_int[i] for i in ids]), np.sqrt(y))[0, 1]
    ms1 = {}; ms1_ppm = np.nan; ms1_tr = None
    for _, r in g[g.MS_LEVEL == 1].iterrows():
        mz, rt, im, it = decode_row(r); m = np.abs(rt - a) <= 4.0
        k = int(str(r.ANNOTATION).rsplit("_i", 1)[-1]) if "_i" in str(r.ANNOTATION) else 0
        ms1[k] = it[m].sum() if m.any() else 0.0
        if k == 0 and m.any():
            ms1_ppm = np.average((mz[m] - r.TARGET_MZ) / r.TARGET_MZ * 1e6, weights=it[m])
            ms1_tr = pd.Series(it[m]).groupby(np.round(rt[m], 3)).sum()
    # MS1-to-fragment co-elution: correlate MS1 mono trace with summed fragment trace (nearest-RT match)
    ms1_co = np.nan
    if ms1_tr is not None and len(ms1_tr) >= 3 and traces:
        fs = pd.DataFrame(traces).fillna(0.0).sum(axis=1).sort_index()
        idx = np.searchsorted(fs.index.values, ms1_tr.index.values).clip(0, len(fs) - 1)
        ms1_co = np.corrcoef(ms1_tr.values, fs.values[idx])[0, 1] if len(fs) >= 3 else np.nan
    rows.append(dict(precursor_id=pid, n_frag=len(frag), frac_sig=len(traces) / max(len(frag), 1), coelution=co, lib_corr=libc,
                     abs_frag_ppm=abs(np.median(ppm)) if ppm else np.nan, ms1_present=float(ms1.get(0, 0) > 0),
                     ms1_iso=float(ms1.get(0, 0) > 0 and ms1.get(1, 0) > 0), abs_ms1_ppm=abs(ms1_ppm), ms1_coelution=ms1_co))
pd.DataFrame(rows).merge(sel, on="precursor_id").to_csv(outp, sep="\t", index=False)
print("wrote", outp, len(rows))
