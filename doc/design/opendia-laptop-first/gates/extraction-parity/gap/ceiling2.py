"""Ceiling as ceiling.py (labels used only to train a diagnostic, never a method), adding the evidence.py features
(SE = S + 9 evidence features; E alone). Same positives, negatives, folds and FPR points."""
import numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
D = B + '/gap-20261114/armE_raw01'
exec(open(B + '/gap-20261114/funnel.py').read().split("# Why stage 1 loses")[0].split('print(\'\\nFunnel')[0])
wr = w.set_index('precursor_id')
pos = np.array(sorted(near_is_winner)); rej = set(near_is_winner - accepted)
dw = w[w.is_decoy == 1]; fpr_cut = float((dw.s >= cut).mean())
rows = np.r_[wr.loc[pos, 'row'].values, dw.row.values]; y = np.r_[np.ones(len(pos)), np.zeros(len(dw))]
is_rej = np.r_[np.isin(pos, list(rej)), np.zeros(len(dw), bool)]
XIC = np.load(f'{D}/xic.npy', mmap_mode='r'); MS1 = np.load(f'{D}/ms1.npy', mmap_mode='r'); PPM = np.load(f'{D}/ppm.npy', mmap_mode='r')
def corr(a, b):
    a = a - a.mean(-1, keepdims=True); b = b - b.mean(-1, keepdims=True)
    return (a * b).sum(-1) / np.sqrt((a * a).sum(-1) * (b * b).sum(-1) + 1e-12)
order = np.argsort(rows); sr = rows[order]; E = np.zeros((len(rows), 9))
import warnings; warnings.filterwarnings('ignore')
for b in range(0, len(sr), 50000):
    r = sr[b:b + 50000]
    x = np.asarray(XIC[r], np.float64)[:, :, 10:22]; m = np.asarray(MS1[r], np.float64)[:, :, 10:22]; p = np.asarray(PPM[r], np.float64)[:, :, 10:22]
    tot = x.sum(1); area = x.sum(2); has = area > 0
    fppm = np.where(has, p.sum(2) / np.maximum(area, 1e-12), np.nan)
    wmean = np.nansum(fppm * area, 1) / np.maximum(area.sum(1), 1e-12)
    E[order[b:b + 50000]] = np.c_[corr(m[:, 0], tot), corr(m[:, 1], m[:, 0]), np.log10(1 + m[:, 0].sum(1)), np.abs(wmean),
        np.where(has.sum(1) >= 2, np.nanstd(fppm, 1), 20.0), np.log10(1 + tot.sum(1)), ((x > 0).sum(2) >= 3).sum(1),
        area.max(1) / np.maximum(area.sum(1), 1e-12), np.nan_to_num(np.nanmean(np.where(has, corr(x, tot[:, None, :] - x), np.nan), 1))]
S = X[rows]; cur = w.set_index('row').loc[rows, 's'].values
def report(name, p):
    pn = p[y == 0]; out = []
    for k in (1, 4):
        thr = np.quantile(pn, 1 - k * fpr_cut); out.append(f'{k}x FPR: all {np.mean(p[y == 1] > thr):.3f} rejected {np.mean(p[is_rej] > thr):.3f}')
    print(f'  {name:8s} ' + ' | '.join(out))
report('D11FDWT', cur)
skf = StratifiedKFold(5, shuffle=True, random_state=0)
for name, F in {'S': S, 'SE': np.c_[S, E], 'E': E, 'S+D11': np.c_[S, cur], 'SE+D11': np.c_[S, E, cur]}.items():
    p = np.zeros(len(y))
    for tr, te in skf.split(F, y):
        mdl = HistGradientBoostingClassifier(max_iter=400, learning_rate=0.05, class_weight='balanced', random_state=0)
        mdl.fit(F[tr], y[tr]); p[te] = mdl.decision_function(F[te])
    report(name, p)
