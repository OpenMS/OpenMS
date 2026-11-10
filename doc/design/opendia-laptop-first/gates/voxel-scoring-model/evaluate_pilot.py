"""Voxel pilot step 4: PLAN_PILOT.md reading of the 8 cells. Only this script reads eval_labels.tsv.
Per cell on held-out folds, at the base winning window of each sample precursor:
  AUC present over decoy, AUC entrapment over decoy (groups from make_sample.py),
  N at FDP_c 1 % and 2 % over the sample target windows (present + entrap + low_target), FDP_c = N_e (1 + 1/r) / N.
Pass (frozen): |AUC_entrap - 0.5| <= 0.01, AUC_present >= (a, T0), N@2 % > (a, T0) with the paired bootstrap 90 %
interval of the difference (1000 precursor resamples, seed 5) above 0. A proxy: the sample is not the full run."""
import os, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
O = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110'
R = 148983 / 151239
lab = pd.read_csv(O + '/eval_labels.tsv', sep='\t')
def auc(a, b): return mannwhitneyu(a, b).statistic / (len(a) * len(b))
def n_at(s, e, t):
    o = np.argsort(-s, kind='stable'); e = e[o]
    fdp = np.cumsum(e) * (1 + 1 / R) / np.arange(1, len(e) + 1); ok = np.where(fdp <= t)[0]
    return int(ok.max() + 1) if len(ok) else 0
cells = [f'{o}_{t}' for t in ('T0', 'T1') for o in 'abcd']
S = {}
for c in cells:
    p = f'{O}/scores_{c}.tsv'
    if os.path.exists(p): S[c] = pd.read_csv(p, sep='\t').merge(lab, on=['sample_id', 'precursor_id'])
ref = S.get('a_T0'); rng = np.random.default_rng(5)
res = {}
print('cell   AUC_present  AUC_entrap  N@1%  N@2%  dN@2% vs a_T0 [90% CI]  verdict')
for c, d in S.items():
    pres, ent, dec = (d[d.group == g].score.values for g in ('present', 'entrap', 'decoy'))
    tg = d[d.group != 'decoy']
    n1, n2 = n_at(tg.score.values, tg.e.values, 0.01), n_at(tg.score.values, tg.e.values, 0.02)
    ap, ae = auc(pres, dec), auc(ent, dec)
    line = f'{c}  {ap:.4f}  {ae:.4f}  {n1}  {n2}'
    if ref is not None and c != 'a_T0':
        m = tg[['precursor_id', 'score', 'e']].merge(ref[ref.group != 'decoy'][['precursor_id', 'score']].rename(columns={'score': 'r'}), on='precursor_id')
        s, r, e = m.score.values, m.r.values, m.e.values; diffs = []
        for _ in range(1000):
            i = rng.integers(0, len(m), len(m)); diffs.append(n_at(s[i], e[i], 0.02) - n_at(r[i], e[i], 0.02))
        lo, hi = np.percentile(diffs, [5, 95]); rp = auc(ref[ref.group == 'present'].score.values, ref[ref.group == 'decoy'].score.values)
        ok = abs(ae - 0.5) <= 0.01 and ap >= rp and n2 > n_at(r, e, 0.02) and lo > 0
        line += f'  {n2 - n_at(r, e, 0.02):+d} [{lo:+.0f}, {hi:+.0f}]  {"PASS" if ok else "fail"}'
    print(line)
