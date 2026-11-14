"""Diagnostic ceiling (uses DIA-NN labels, evaluation only, never a method): is the evidence needed to accept the
626 right-peak-but-rejected DIA-NN IDs already in OpenDIA's per-window arrays, or is it missing?
Positives: OpenDIA winners at DIA-NN's peak for human DIA-NN q<=0.01 IDs (2200). Negatives: decoy winners.
Classifier: sklearn HistGradientBoosting, 5-fold CV by precursor, on feature sets
  S  the 20 scorer inputs
  SX S + log1p of the 6 fragment XICs and 4 MS1 isotope traces over the central 12 bins (+-5.6 s) and the 6 library intensities
  X  the raw arrays alone (no engineered features)
Reported: TPR on the rejected 626 and on all positives at the decoy-winner rate that the current cut admits,
and the same for the current D11FDWT score. FPR is estimated on the held-out decoy winners."""
import numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
D = B + '/gap-20261114/armE_raw01'
exec(open(B + '/gap-20261114/funnel.py').read().split("# Why stage 1 loses")[0].split('print(\'\\nFunnel')[0])
wr = w.set_index('precursor_id')
pos = np.array(sorted(near_is_winner)); rej = set(near_is_winner - accepted)
dw = w[w.is_decoy == 1]
fpr_cut = float((dw.s >= cut).mean())
print(f'decoy winners {len(dw)}, above cut {int((dw.s >= cut).sum())}, FPR at cut {fpr_cut:.2e}')
rng = np.random.default_rng(0)
neg_rows = dw.row.values
pos_rows = wr.loc[pos, 'row'].values
rows = np.r_[pos_rows, neg_rows]; y = np.r_[np.ones(len(pos_rows)), np.zeros(len(neg_rows))]
is_rej = np.r_[np.isin(pos, list(rej)), np.zeros(len(neg_rows), bool)]
order = np.argsort(rows); inv = np.empty_like(order); inv[order] = np.arange(len(order))
XIC = np.load(f'{D}/xic.npy', mmap_mode='r'); MS1 = np.load(f'{D}/ms1.npy', mmap_mode='r'); LIB = np.load(f'{D}/lib.npy', mmap_mode='r')
sr = rows[order]
xic = np.log1p(np.asarray(XIC[sr][:, :, 10:22], np.float32))[inv].reshape(len(rows), -1)
ms1 = np.log1p(np.asarray(MS1[sr][:, :, 10:22], np.float32))[inv].reshape(len(rows), -1)
lib = np.asarray(LIB[sr], np.float32)[inv]
S = X[rows]
sets = {'S': S, 'SX': np.c_[S, xic, ms1, lib], 'X': np.c_[xic, ms1, lib]}
cur = w.set_index('row').loc[rows, 's'].values
def report(name, p):
    pn = p[y == 0]
    thr = np.quantile(pn, 1 - fpr_cut)
    print(f'  {name:6s} TPR all {np.mean(p[y == 1] > thr):.3f}  rejected {np.mean(p[is_rej] > thr):.3f}  '
          f'(at 4x FPR: all {np.mean(p[y == 1] > np.quantile(pn, 1 - 4 * fpr_cut)):.3f} rejected {np.mean(p[is_rej] > np.quantile(pn, 1 - 4 * fpr_cut)):.3f})')
report('D11FDWT', cur)
skf = StratifiedKFold(5, shuffle=True, random_state=0)
for name, F in sets.items():
    p = np.zeros(len(y))
    for tr, te in skf.split(F, y):
        m = HistGradientBoostingClassifier(max_iter=400, learning_rate=0.05, class_weight='balanced', random_state=0)
        m.fit(F[tr], y[tr]); p[te] = m.decision_function(F[te])
    report(name, p)
