"""Our raw-scan extraction (OpenDIACalibrate -out_anchor_traces: 6 fragments, calibrated m/z, our Da tolerance) vs the
OSW XIC (same precursor, same fragments by library order, ~34 ppm) on the 512 raw 01 anchors. Interpolate OSW onto our
scan RTs. Report per-fragment Pearson, intensity ratio, and the fraction of scans where one side is zero and the other
is not. Large disagreement would mean an extraction problem upstream of any scoring."""
import os, numpy as np, pandas as pd
D = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-design-raw01-20261107'
S = os.environ['JCODE_SCRATCH_DIR'] + '/xic01'
frag = np.load(S + '/frag.npy', mmap_mode='r'); rt = np.load(S + '/rt.npy'); npts = np.load(S + '/npts.npy'); meta = pd.read_parquet(S + '/meta.parquet')
row = pd.Series(np.arange(len(meta)), index=meta.precursor_id)
t = pd.read_csv(D + '/anchor-peaks/anchor_traces.tsv', sep='\t')
res = []
for pid, g in t.groupby('precursor_id'):
    i = row.get(pid)
    if i is None: continue
    n = npts[i]; ort = rt[i, :n]
    for f in range(6):
        ours = g[f'i{f}'].values.astype(float)
        osw = np.interp(g.rt.values, ort, np.asarray(frag[i, f, :n], float), left=np.nan, right=np.nan)
        m = np.isfinite(osw)
        if m.sum() < 10: continue
        a, b = ours[m], osw[m]
        pr = np.corrcoef(a, b)[0, 1] if a.std() > 0 and b.std() > 0 else np.nan
        res.append(dict(pid=pid, f=f, pearson=pr, ratio=np.log2((a.sum() + 1) / (b.sum() + 1)),
                        ours_zero_osw_not=np.mean((a == 0) & (b > 0)), osw_zero_ours_not=np.mean((b == 0) & (a > 0)),
                        ours_nonzero=np.mean(a > 0), osw_nonzero=np.mean(b > 0)))
r = pd.DataFrame(res)
print('anchors compared', r.pid.nunique(), 'fragment traces', len(r))
print(r.drop(columns=['pid']).groupby('f').median().round(3))
print('overall median pearson', round(r.pearson.median(), 3), '| q10', round(r.pearson.quantile(.1), 3), '| median log2(ours/osw)', round(r.ratio.median(), 2))
