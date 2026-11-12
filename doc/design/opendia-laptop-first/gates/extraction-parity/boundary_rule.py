"""How DIA-NN places its apex and boundaries, measured on its own --xic traces (raw01, q <= 1%), and the same
quantities for our picker on our traces. Evaluation only.

For each precursor, on the scan grid of DIA-NN's export:
  S_dn  = sum of DIA-NN's exported fragments, S_us = sum of our 6 fragments (same scans)
  - apex: is RT the raw max of S_dn / of the best fragment / of a [1 2 1]-smoothed S_dn, within RT.Start..RT.Stop?
  - boundaries: S_dn(start) / S_dn(apex), S_dn(stop) / S_dn(apex); scans apex-start vs stop-apex;
    width in scans vs FWHM in scans; is the boundary at a local minimum?
  - our boundary heights on S_us for our matched peak (from peak_parity_rows_raw01.tsv)
Usage: boundary_rule.py <run dir>
"""
import sys, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; X = B + '/extraction-parity-20261112'
RUN = sys.argv[1] if len(sys.argv) > 1 else 'rtmono'
out = open(X + '/' + RUN + '/BOUNDARY_RULE_raw01.txt', 'w')
def say(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')

rows = pd.read_csv(X + '/' + RUN + '/peak_parity_rows_raw01.tsv', sep='\t')
rows = rows[(rows.grp == 'q<=1%') & rows.contains]
idm = pd.read_csv(X + '/idmap.tsv', sep='\t')
rows = rows.merge(idm[['precursor_id', 'Precursor.Id']], on='precursor_id')
dn = pd.read_parquet(X + '/diann/report_raw01_xic/raw01.xic.parquet')
dn = dn[dn.pr.isin(set(rows['Precursor.Id'])) & ~dn.feature.isin(['index', 'ms1'])].copy()
dn['rt_s'] = dn.rt.astype(float) * 60
us = pd.read_csv(X + '/' + RUN + '/traces.tsv', sep='\t', usecols=['precursor_id', 'rt', 'slot', 'intensity'])
us = us[us.precursor_id.isin(set(rows.precursor_id))]
usS = us.groupby(['precursor_id', 'rt']).intensity.sum()
dnS = dn.groupby(['pr', 'rt_s']).value.sum()
dnF = dn.pivot_table(index=['pr', 'rt_s'], columns='feature', values='value', aggfunc='sum')

def at(rt_arr, val, x):
    i = int(np.argmin(np.abs(rt_arr - x)));
    return i, (abs(rt_arr[i] - x) < 0.01)

res = []
for r in rows.itertuples():
    try:
        s = dnS.loc[r._21]
    except KeyError:
        continue
    t = s.index.values; v = s.values
    ia, oka = at(t, v, r.RT); i0, ok0 = at(t, v, r._5); i1, ok1 = at(t, v, r._6)
    if not (oka and ok0 and ok1) or v[ia] <= 0:
        continue
    sm = np.convolve(np.pad(v, 1, mode='edge'), [0.25, 0.5, 0.25], mode='valid')
    F = dnF.loc[r._21].reindex(t).fillna(0).values  # scans x fragments
    best = F[:, np.argmax(F[i0:i1 + 1].sum(0))]
    seg = slice(i0, i1 + 1)
    # our consensus on the same scans
    try:
        u = usS.loc[r.precursor_id]
    except KeyError:
        continue
    ut = u.index.values; uv = u.values
    def uat(x):
        j = int(np.argmin(np.abs(ut - x))); return uv[j] if abs(ut[j] - x) < 0.01 else np.nan
    uapex = uat(r.rt_apex)
    res.append(dict(precursor_id=r.precursor_id,
        apex_is_rawmax=int(np.argmax(v[seg])) + i0 == ia, apex_is_smmax=int(np.argmax(sm[seg])) + i0 == ia,
        apex_is_bestfrag_max=int(np.argmax(best[seg])) + i0 == ia,
        left_sc=ia - i0, right_sc=i1 - ia, w_sc=i1 - i0, fwhm_sc=r.FWHM / 1.646,
        h_start=v[i0] / v[ia], h_stop=v[i1] / v[ia],
        h_start_best=best[i0] / max(best[ia], 1e-9), h_stop_best=best[i1] / max(best[ia], 1e-9),
        start_is_min=(i0 == 0 or v[i0] <= v[i0 - 1]) and v[i0] <= v[i0 + 1],
        stop_is_min=(i1 == len(v) - 1 or v[i1] <= v[i1 + 1]) and v[i1] <= v[i1 - 1],
        # height of DIA-NN's boundaries on our consensus, and of our boundaries on our consensus
        us_h_dn_start=uat(r._5) / uat(r.RT) if uat(r.RT) else np.nan, us_h_dn_stop=uat(r._6) / uat(r.RT) if uat(r.RT) else np.nan,
        us_h_start=uat(r.rt_start) / uapex if uapex else np.nan, us_h_stop=uat(r.rt_end) / uapex if uapex else np.nan,
        us_left_sc=round((r.rt_apex - r.rt_start) / 1.646), us_right_sc=round((r.rt_end - r.rt_apex) / 1.646),
        # beyond DIA-NN's stop: does our consensus keep falling or rise (tail vs. next peak)?
        us_after_dn_stop=uat(r._6 + 1.646) / uat(r._6) if uat(r._6) else np.nan))
R = pd.DataFrame(res)
def q(s): return {k: round(float(v), 3) for k, v in s.quantile([.1, .25, .5, .75, .9]).items()}
say(f'q<=1% precursors with DIA-NN apex inside our peak and DIA-NN apex/start/stop on export scans: {len(R)}')
say(f'DIA-NN apex = raw max of S_dn in [start,stop]: {R.apex_is_rawmax.mean():.4f}; = [1 2 1] smoothed max: {R.apex_is_smmax.mean():.4f}; = best-fragment max: {R.apex_is_bestfrag_max.mean():.4f}')
say('DIA-NN scans apex-start', q(R.left_sc), ' stop-apex', q(R.right_sc), ' width', q(R.w_sc))
say('DIA-NN width share by scans:', R.w_sc.value_counts().sort_index().to_dict())
say('DIA-NN (stop-apex) - (apex-start):', (R.right_sc - R.left_sc).value_counts().sort_index().to_dict())
say('DIA-NN FWHM in scans', q(R.fwhm_sc), ' width / FWHM', q(R.w_sc / R.fwhm_sc))
say('DIA-NN boundary height / apex on S_dn: start', q(R.h_start), ' stop', q(R.h_stop))
say('   on best fragment: start', q(R.h_start_best), ' stop', q(R.h_stop_best))
say(f'DIA-NN boundary is a local minimum of S_dn: start {R.start_is_min.mean():.4f}, stop {R.stop_is_min.mean():.4f}')
say('our scans apex-start', q(R.us_left_sc), ' end-apex', q(R.us_right_sc))
say('on our consensus: our boundary height / apex: start', q(R.us_h_start), ' end', q(R.us_h_stop))
say('                  DIA-NN boundary height / DIA-NN apex: start', q(R.us_h_dn_start), ' stop', q(R.us_h_dn_stop))
say('our consensus one scan after DIA-NN stop / at DIA-NN stop', q(R.us_after_dn_stop))
R.to_csv(X + '/' + RUN + '/boundary_rule_rows_raw01.tsv', sep='\t', index=False)
