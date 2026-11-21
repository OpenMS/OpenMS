"""Why is H's RT calibration worse than D's? (diagnostic)
Each run's calibration.json holds the isotonic RT map fitted on that run's own anchors. Compare it with DIA-NN's RT for the
DIA-NN IDs of the same run: residual by gradient position, map knots and slopes, anchor coverage, and the RT range of the
DIA-NN IDs against the map's knot range. Also: how does DIA-NN's own RT relate to library RT (linear? curved?)."""
import json, glob, sqlite3, numpy as np, pandas as pd
C = '/nfs/research/juan/singj/opendia'
LIB = f'{C}/entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp'
c = sqlite3.connect(f'file:{LIB}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID pid, p.DECOY d, p.CHARGE z, p.LIBRARY_RT lrt, pe.MODIFIED_SEQUENCE ms,
  min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c).set_index('pid')
tg = pq[pq.d == 0]; sk = (tg.ms + '/' + tg.z.astype(str)).reset_index().drop_duplicates(0).set_index(0).pid
print(f'library RT range {tg.lrt.min():.3f} .. {tg.lrt.max():.3f} (median {tg.lrt.median():.3f})')

def apply_map(cal, lrt):
    r = cal['rt']; kx = np.asarray(r['knots_library_rt']); ky = np.asarray(r['knots_run_rt_s'])
    y = np.interp(lrt, kx, ky)
    lo, hi = lrt < kx[0], lrt > kx[-1]
    y[lo] = ky[0] + (lrt[lo] - kx[0]) * r['low_slope_s']; y[hi] = ky[-1] + (lrt[hi] - kx[-1]) * r['high_slope_s']
    return y

for run, rd in (('D', f'{C}/entlib-mt-20261120/runs/D_mt_s0'), ('H', f'{C}/entlib-mt-20261120/runs/H_mt_s0')):
    cal = json.load(open(glob.glob(f'{rd}/*.calibration.json')[0])); r = cal['rt']
    d = pd.read_parquet(f'{C}/entlib-diann5-20261119/{run}/report.parquet', columns=['Modified.Sequence', 'Precursor.Charge', 'Q.Value', 'Decoy', 'RT', 'Precursor.Quantity'])
    d = d[(d['Q.Value'] <= 0.01) & (d.Decoy == 0)].copy()
    d['pid'] = sk.reindex(d['Modified.Sequence'] + '/' + d['Precursor.Charge'].astype(str)).values
    d = d.dropna(subset=['pid']); d['pid'] = d.pid.astype(np.int64)
    d = d[pq.e.reindex(d.pid).values == 0].drop_duplicates('pid').set_index('pid')
    d['rt_s'] = d.RT * 60; d['lrt'] = tg.lrt.reindex(d.index).values
    d['pred'] = apply_map(cal, d.lrt.values); d['res'] = d.rt_s - d.pred
    kx = np.asarray(r['knots_library_rt']); ky = np.asarray(r['knots_run_rt_s'])
    print(f'\n===== {run}: map {r["map"]}, {len(kx)} knots, anchors {r["anchors_used"]}, tolerance {r["tolerance_half_width_s"]:.1f} s, sigma {r["sigma_s"]:.1f} s')
    print(f'  knots library RT {kx.min():.3f}..{kx.max():.3f} -> run {ky.min():.0f}..{ky.max():.0f} s; tail slopes low {r["low_slope_s"]:.0f} high {r["high_slope_s"]:.0f} s/unit')
    print(f'  DIA-NN IDs {len(d)}: run RT {d.rt_s.min():.0f}..{d.rt_s.max():.0f} s; library RT {d.lrt.min():.3f}..{d.lrt.max():.3f}')
    inside = (d.lrt >= kx.min()) & (d.lrt <= kx.max())
    print(f'  share of DIA-NN IDs inside the knot range: {inside.mean():.3f}')
    a = d.res.abs()
    print(f'  |DIA-NN RT - map| p50 {a.median():.1f} p90 {a.quantile(.9):.1f} p95 {a.quantile(.95):.1f} s; inside knots p95 {a[inside].quantile(.95):.1f}, outside p95 {a[~inside].quantile(.95) if (~inside).any() else float("nan"):.1f}')
    # residual by gradient decile of DIA-NN RT
    d['dec'] = pd.qcut(d.rt_s, 10, labels=False)
    g = d.groupby('dec').agg(rt_lo=('rt_s', 'min'), rt_hi=('rt_s', 'max'), n=('res', 'size'), bias=('res', 'median'),
                            p95=('res', lambda x: x.abs().quantile(.95)))
    print('  residual by DIA-NN RT decile:\n' + g.round(1).to_string())
    # how well can ANY monotone map do? isotonic fit of DIA-NN RT on library RT (upper bound for this library)
    from sklearn.isotonic import IsotonicRegression
    o = np.argsort(d.lrt.values); iso = IsotonicRegression(out_of_bounds='clip').fit(d.lrt.values[o], d.rt_s.values[o])
    best = (d.rt_s - iso.predict(d.lrt.values)).abs()
    print(f'  best monotone map fitted on DIA-NN IDs themselves: |res| p50 {best.median():.1f} p95 {best.quantile(.95):.1f} s  (library RT noise floor)')
    # linear fit for reference
    A = np.vstack([d.lrt.values, np.ones(len(d))]).T; coef = np.linalg.lstsq(A, d.rt_s.values, rcond=None)[0]
    lin = (d.rt_s - A @ coef).abs()
    print(f'  linear map on DIA-NN IDs: slope {coef[0]:.0f} s/unit, |res| p50 {lin.median():.1f} p95 {lin.quantile(.95):.1f} s')
    # residual scaled by gradient length: in units of the run's elution span
    span = d.rt_s.quantile(.99) - d.rt_s.quantile(.01)
    print(f'  elution span p1..p99 {span:.0f} s; map p95 residual = {a.quantile(.95) / span * 100:.1f} % of span; noise floor {best.quantile(.95) / span * 100:.1f} % of span')
