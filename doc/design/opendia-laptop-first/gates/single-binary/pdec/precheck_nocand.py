"""PRECHECK_NOCAND (diagnostic): why do 17.8 % of DIA-NN IDs at r~7 D have no OpenDIA candidate near DIA-NN's RT?
Data: r~7 D s0 default (07131138b8, entlib-feat36-20261119), calibration.json of that run, library
hela_r5.predicted_decoys.pqp. DIA-NN 2.2 RT/IM for evaluation only (non-entrapment q<=0.01).
Sorts each missed ID into:
  rt_out   DIA-NN RT outside the calibrated search window (map(lib RT) +- rt tolerance)
  im_out   inside the RT window, DIA-NN 1/K0 outside the calibrated mobility window
  cap      inside both, and the precursor already has max_peaks (5) candidates elsewhere
  picker   inside both, fewer than 5 candidates (the picker found no peak at DIA-NN's RT)"""
import glob, json, sqlite3, numpy as np, pandas as pd
G = '/home/singjc/.jcode/scratch/dia2/gap7'
L = '/home/singjc/.jcode/scratch/dia2/entlib10/cut/hela_r5.predicted_decoys.pqp'
DN = '/home/singjc/.jcode/scratch/dia2/entlib12/D/report.parquet'
cal = json.load(open(f'{G}/D36.calibration.json'))
kx, ky = np.array(cal['rt']['knots_library_rt']), np.array(cal['rt']['knots_run_rt_s'])
lo, hi, tol = cal['rt']['low_slope_s'], cal['rt']['high_slope_s'], cal['rt']['tolerance_half_width_s']
im_a, im_b, im_tol = cal['im']['intercept_vssc'], cal['im']['slope'], cal['im']['tolerance_half_width_vssc']
def rtmap(x):
    y = np.interp(x, kx, ky)
    y = np.where(x < kx[0], ky[0] + lo * (x - kx[0]), y)
    return np.where(x > kx[-1], ky[-1] + hi * (x - kx[-1]), y)
F = pd.read_parquet(f'{G}/D36_feat.parquet', columns=['precursor_id', 'is_decoy', 'rt_apex', 'score', 'transition_support'])
c = sqlite3.connect(f'file:{L}?mode=ro', uri=True)
tg = pd.read_sql("""select p.ID pid, p.CHARGE z, p.LIBRARY_RT lrt, p.LIBRARY_DRIFT_TIME lim, pe.MODIFIED_SEQUENCE ms,
  min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID where p.DECOY=0 group by p.ID""", c)
tg = tg[tg.e == 0].assign(sk=tg.ms + '/' + tg.z.astype(str)).drop_duplicates('sk').set_index('sk')
d = pd.read_parquet(DN, columns=['Modified.Sequence', 'Precursor.Charge', 'Q.Value', 'RT', 'IM', 'Precursor.Quantity'])
d = d[d['Q.Value'] <= 0.01]
k = d['Modified.Sequence'] + '/' + d['Precursor.Charge'].astype(str)
d = d[k.isin(tg.index).values].copy(); k = k[k.isin(tg.index)]
for col in ('pid', 'lrt', 'lim'): d[col] = tg[col].reindex(k).values
d = d.drop_duplicates('pid').set_index('pid')
d['rt_s'] = d.RT * 60; d['pred'] = rtmap(d.lrt.values); d['rt_err'] = d.rt_s - d.pred
d['im_pred'] = im_a + im_b * d.lim; d['im_err'] = d.IM - d.im_pred
T = F[(F.is_decoy == 0) & F.precursor_id.isin(d.index)].copy()
T['drt'] = np.abs(T.rt_apex.values - d.rt_s.reindex(T.precursor_id).values)
g = T.groupby('precursor_id')
d['n'] = g.size().reindex(d.index).fillna(0).astype(int)
d['right'] = g.drt.min().reindex(d.index) <= 6
d['qq'] = pd.qcut(np.log10(d['Precursor.Quantity'].clip(lower=1)), 4, labels=['Q1', 'Q2', 'Q3', 'Q4'])
print(f'DIA-NN IDs {len(d)}, no candidate within 6 s {(~d.right).sum()} ({(~d.right).mean():.3f})')
print(f'RT tolerance +-{tol:.1f} s, IM tolerance +-{im_tol:.4f}, map anchors {cal["rt"]["anchors_used"]}')
print(f'All DIA-NN IDs: |RT error| p50 {d.rt_err.abs().median():.1f} p95 {d.rt_err.abs().quantile(.95):.1f} s; '
      f'outside RT window {(d.rt_err.abs() > tol).mean():.3f}; |IM error| p95 {d.im_err.abs().quantile(.95):.4f}, outside {(d.im_err.abs() > im_tol).mean():.3f}')
M = d[~d.right].copy()
rt_out = M.rt_err.abs() > tol; im_out = ~rt_out & (M.im_err.abs() > im_tol)
cap = ~rt_out & ~im_out & (M.n >= 5)
M['cat'] = np.select([rt_out, im_out, cap], ['rt_out', 'im_out', 'cap (5 candidates elsewhere)'], 'picker (<5 candidates, none at DIA-NN RT)')
tab = pd.crosstab(M.cat, M.qq); tab['all'] = M.cat.value_counts(); tab['share'] = (tab['all'] / len(M)).round(3)
print('\nMissed IDs by cause and DIA-NN quantity quartile:'); print(tab.to_string())
print(f'\nCandidates per precursor: missed p50 {M.n.median():.0f}, found p50 {d[d.right].n.median():.0f}; '
      f'share with 5: missed {(M.n >= 5).mean():.3f}, found {(d[d.right].n >= 5).mean():.3f}')
for name, s in (('rt_out', M[rt_out]), ('cap', M[cap])):
    if len(s): print(f'{name}: signed RT error p5/p50/p95 {s.rt_err.quantile(.05):.1f}/{s.rt_err.median():.1f}/{s.rt_err.quantile(.95):.1f} s; '
                     f'library RT p5/p50/p95 {s.lrt.quantile(.05):.3f}/{s.lrt.median():.3f}/{s.lrt.quantile(.95):.3f}')
P = M[cap | (M.cat.str.startswith('picker'))]
TT = T[T.precursor_id.isin(P.index)]
print(f'In-window misses: nearest candidate to DIA-NN RT p50 {TT.groupby("precursor_id").drt.min().median():.1f} s; '
      f'their candidates span RT +-{np.median(TT.groupby("precursor_id").rt_apex.agg(lambda x: x.max() - x.min())) / 2:.1f} s (p50 half-span)')
# Gradient position of misses vs found
d['pos'] = pd.cut(d.rt_s, [0, 120, 180, 240, 300, 1e9], labels=['<2 min', '2-3', '3-4', '4-5', '>5'])
print('\nMiss share by DIA-NN RT:', ' '.join(f'{p} {(~d.right[d.pos == p]).mean():.3f} (n {int((d.pos == p).sum())})' for p in d.pos.cat.categories))
