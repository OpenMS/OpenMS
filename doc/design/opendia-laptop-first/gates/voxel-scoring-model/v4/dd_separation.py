"""Stage-D deep dive, step 2 (descriptive; labels only define the groups). Uses dd_feats_<raw>.tsv (dd_evidence.py).
(1) Separation per raw-evidence measure: AUC of D against T1 (top-1 % decoys) and of E against T1. A measure that
    separates D from T1 well but is not used (or used weakly) by the model is a candidate cause; one with AUC ~0.5 shows
    what D shares with decoys.
(2) Per-miss 'decoy-like' profile: for each D window, which measures fall inside the T1 inter-quartile range (decoy-like)
    versus the E inter-quartile range (target-like). Count the patterns, so each miss gets its own reason.
(3) Theoretical isotope check: iso_ratio (M+1/M+0) against the averagine expectation for the precursor mass; decoys share
    the target's composition, so a mismatch marks wrong-precursor MS1 signal (co-isolated other ion)."""
import sys, sqlite3, numpy as np, pandas as pd
R = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; O = R + '/astral-voxel-pilot-20261110/v2'
M = ['n_sig', 'n_peak', 'log_area', 'frac_core', 'corr_med', 'top1_share', 'lib_cos', 'rank_agree', 'ppm_abs', 'ms1_m0', 'ms1_corr',
     'iso_dev', 'n_sig12', 'extra_sig']
HIGH_GOOD = {'n_sig': 1, 'n_peak': 1, 'log_area': 1, 'frac_core': 1, 'corr_med': 1, 'top1_share': -1, 'lib_cos': 1, 'rank_agree': 1,
             'ppm_abs': -1, 'ms1_m0': 1, 'ms1_corr': 1, 'iso_dev': -1, 'n_sig12': 1, 'extra_sig': 1}
c = sqlite3.connect('file:' + R + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
pre = pd.read_sql('select ID precursor_id, PRECURSOR_MZ mz, CHARGE z from PRECURSOR', c).set_index('precursor_id')


def auc(a, b):
  a, b = np.asarray(a, float), np.asarray(b, float); a, b = a[np.isfinite(a)], b[np.isfinite(b)]
  r = pd.Series(np.r_[a, b]).rank().values; return (r[:len(a)].sum() - len(a) * (len(a) + 1) / 2) / (len(a) * len(b))


for raw in sys.argv[1:] or ['raw01', 'raw02']:
  f = pd.read_csv(f'{O}/v4/dd_feats_{raw}.tsv', sep='\t').set_index(['group', 'row'])
  pop = pd.read_csv(f'{O}/{raw}/population.tsv', sep='\t', usecols=['precursor_id'])
  pid = pop.precursor_id.values[f.index.get_level_values('row')]
  mass = (pre.mz.reindex(pid).values - 1.007276) * pre.z.reindex(pid).values
  exp = 0.000594 * mass - 0.03091  # averagine-like linear fit of the M+1/M+0 ratio (about 0.55 at 1 kDa, 1.15 at 2 kDa)
  f['iso_dev'] = np.abs(f.iso_ratio.values - exp)
  g = {k: f.xs(k, level='group') for k in f.index.get_level_values('group').unique()}
  D, E, T = g['D below cut'], g['E accepted'], g['T1 top-1 % decoys']
  print(f'\n== {raw} (1) AUC oriented so 1 = target-like (D vs T1, E vs T1, D vs E)')
  rows = []
  for m in M:
    s = HIGH_GOOD[m]; rows.append((m, auc(s * D[m], s * T[m]), auc(s * E[m], s * T[m]), auc(s * D[m], s * E[m])))
  print(pd.DataFrame(rows, columns=['measure', 'D_vs_T1', 'E_vs_T1', 'D_vs_E']).round(3).to_string(index=False))
  # (2) per-miss decoy-like counts: a measure is decoy-like when the D value is on the decoy side of the midpoint
  #     between the E and T1 medians.
  mid = {m: (E[m].median() + T[m].median()) / 2 for m in M}
  like = pd.DataFrame({m: (HIGH_GOOD[m] * (D[m] - mid[m]) < 0) for m in M})
  print(f'{raw} (2) share of D misses that are decoy-side on each measure:')
  print(like.mean().round(3).sort_values(ascending=False).to_string())
  nlike = like.sum(1)
  print(f'{raw} decoy-side measures per D miss (of {len(M)}): quantiles', nlike.quantile([0.1, 0.25, 0.5, 0.75, 0.9]).tolist())
  for name, grp in (('E', E), ('T1', T)):
    l2 = pd.DataFrame({m: (HIGH_GOOD[m] * (grp[m] - mid[m]) < 0) for m in M}).sum(1)
    print(f'{raw} same count for {name}: median {l2.median()}, q25 {l2.quantile(0.25)}, q75 {l2.quantile(0.75)}')
