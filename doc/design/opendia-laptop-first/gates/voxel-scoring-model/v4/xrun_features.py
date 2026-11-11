"""PLAN_V7 cross-run features (label-free, model-free). Usage: xrun_features.py <raw01|raw02>.
For each population window w of precursor p in run a, the other run b's window of the same precursor p whose apex,
mapped to run a, is nearest to w's apex and within 10 s (ties: lower sample_id). Map b -> a: binned medians (20
equal-count bins, linear between knots, constant ends) of run_rt of the calibration anchors used (anchors.tsv used=1)
that are shared by both runs. These come from OpenDIACalibrate, so no model score and no label is used.
Output <raw>/featX.npy [n, 19] float32: the 11 heuristics and the 6 feat12 values of that other-run window, its
|mapped apex - w apex| in s, and a presence flag. When no such window exists: presence 0, |drt| 10, and the 17 feature
columns are set to their median over the present rows (so they carry no signal after standardisation)."""
import sys, numpy as np, pandas as pd
A = sys.argv[1]; Bn = {'raw01': 'raw02', 'raw02': 'raw01'}[A]
R = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
O = R + '/astral-voxel-pilot-20261110/v2'
CAL = {'raw01': R + '/astral-design-raw01-20261107/window300k', 'raw02': R + '/astral-entrapment-raw02-20261108/cal'}
TOL = 10.0
an = {r: pd.read_csv(CAL[r] + '/anchors.tsv', sep='\t').query('used == 1').drop_duplicates('precursor_id').set_index('precursor_id').run_rt
      for r in (A, Bn)}
sh = an[A].to_frame('a').join(an[Bn].rename('b'), how='inner')
o = np.argsort(sh.b.values); xb, ya = sh.b.values[o], sh.a.values[o]
e = np.array_split(np.arange(len(xb)), 20)
kx = np.array([np.median(xb[i]) for i in e]); ky = np.maximum.accumulate(np.array([np.median(ya[i]) for i in e]))
res = np.abs(ya - np.interp(xb, kx, ky))
print(f'{A}: shared anchors {len(sh)}; map {Bn} -> {A} residual median {np.median(res):.2f} s, q95 {np.quantile(res, 0.95):.2f} s')
pa = pd.read_csv(f'{O}/{A}/population.tsv', sep='\t', usecols=['sample_id', 'precursor_id', 'rt_apex'])
pb = pd.read_csv(f'{O}/{Bn}/population.tsv', sep='\t', usecols=['sample_id', 'precursor_id', 'rt_apex'])
Fb = np.c_[np.load(f'{O}/{Bn}/heuristics.npy'), np.load(f'{O}/{Bn}/feat12.npy')].astype(np.float64)
pa['row'] = np.arange(len(pa)); pb['rowb'] = np.arange(len(pb)); pb['rt_a'] = np.interp(pb.rt_apex.values, kx, ky)
m = pa.merge(pb[['precursor_id', 'rowb', 'sample_id', 'rt_a']].rename(columns={'sample_id': 'sid_b'}), on='precursor_id', how='inner')
m['d'] = (m.rt_a - m.rt_apex).abs(); m = m[m.d <= TOL].sort_values(['row', 'd', 'sid_b']).drop_duplicates('row')
out = np.empty((len(pa), 19), np.float64)
present = np.zeros(len(pa), bool); present[m.row.values] = True
med = np.median(Fb[m.rowb.values], axis=0)
out[:, :17] = med; out[m.row.values, :17] = Fb[m.rowb.values]
out[:, 17] = TOL; out[m.row.values, 17] = m.d.values
out[:, 18] = present
np.save(f'{O}/{A}/featX.npy', out.astype(np.float32))
pop = pd.read_csv(f'{O}/{A}/population.tsv', sep='\t', usecols=['is_decoy'])
print(f'{A}: windows {len(pa)}, with other-run window {present.mean():.3f} (targets {present[pop.is_decoy.values == 0].mean():.3f}, '
      f'decoys {present[pop.is_decoy.values == 1].mean():.3f}); |drt| median of present {np.median(m.d):.2f} s')
