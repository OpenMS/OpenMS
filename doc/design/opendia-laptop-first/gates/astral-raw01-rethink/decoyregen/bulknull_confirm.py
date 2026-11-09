"""Confirmation of the bulk-null density-ratio TDC (k6b, declared before running).
bulknull.py passed its rule on both runs with 2 of 3 seeds each (marginal). Fixed configuration for confirmation:
score = mean of 5 MLP seeds (0,10,20,30,40; the calib_mitig ensemble), bulk = winners below the median decoy-winner
score, weights from histogram gradient boosting (bulk targets vs bulk decoys, pair-fold cross-fitted), clip [0.2, 5],
renormalised to mean 1 over bulk decoys. Primary rule: pass on a run if weighted q<=0.01 FDP_c <= 0.015 AND weighted
q<=0.05 FDP_c <= 0.06; confirmed if it passes on both runs.
Sensitivity variants (reported, no rule): bulk threshold at the 30th / 70th decoy percentile; clip [0.1, 10]; weights
from a logistic regression instead of boosting. Also: N at q<=0.01, entrapment-matched N at FDP_c 1%, the lower-bound
FDP (N_e/N) and the mean decoy weight in the tail.
Usage: OPENDIA_RUN=R02."""
import os, sys, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/../gap2')
import oswlike_score as S
_sel = S.select
def select(idx, score, fdr):
    pos, neg = _sel(idx, score, fdr)
    if len(pos) < 200:
        w = S.winners(idx, score); tw = w[~S.dec[w]]; pos = tw[np.argsort(-score[tw], kind='stable')[:200]]
    return pos, neg
S.select = select
meta = S.m; fam = meta.family.values; keep = (fam == 'T') | (fam == 'O')
for k_, v in {k: getattr(S, k).copy() for k in ('X0', 'prec', 'dec', 'fold', 'pair', 'entr', 'right')}.items(): setattr(S, k_, v[keep])
X = np.nan_to_num(S.X0); Xz = (X - X.mean(0)) / (X.std(0) + 1e-9)
def tdc_w(wn, score, wdec):
    d = pd.DataFrame({'i': wn, 'p': S.pair[wn], 's': score[wn], 'd': S.dec[wn]})
    d = d.sort_values(['s', 'd'], ascending=[False, False]).drop_duplicates('p').sort_values('s', ascending=False)
    dd = d.d.values; T = np.cumsum(~dd); Dn = np.cumsum(np.where(dd, wdec[d.i.values], 0.0))
    return d.i.values, np.minimum.accumulate(np.minimum(1, (Dn + 1) / np.maximum(T, 1))[::-1])[::-1]
def weights(sc, wn, pct=50, clip=(0.2, 5), kind='hgb'):
    dw = wn[S.dec[wn]]; tw = wn[~S.dec[wn]]; thr = np.percentile(sc[dw], pct)
    bt, bd = tw[sc[tw] < thr], dw[sc[dw] < thr]; idx = np.r_[bt, bd]; y = np.r_[np.ones(len(bt)), np.zeros(len(bd))]
    pall = np.zeros(len(sc)); poof = np.zeros(len(idx))
    for k in range(3):
        trm = S.fold[idx] != k
        m = (HistGradientBoostingClassifier(max_iter=200, learning_rate=0.05, max_leaf_nodes=31, random_state=0) if kind == 'hgb'
             else LogisticRegression(max_iter=500, C=1.0)).fit((X if kind == 'hgb' else Xz)[idx[trm]], y[trm])
        te = dw[S.fold[dw] == k]; pall[te] = m.predict_proba((X if kind == 'hgb' else Xz)[te])[:, 1]
        poof[~trm] = m.predict_proba((X if kind == 'hgb' else Xz)[idx[~trm]])[:, 1]
    r = np.clip(pall[dw] / np.maximum(1 - pall[dw], 1e-9) * len(bd) / len(bt), *clip)
    w = np.ones(len(sc)); w[dw] = r / r[sc[dw] < thr].mean()
    return w, roc_auc_score(y, poof)
def report(name, sc, wn, w):
    kept, q = tdc_w(wn, sc, w); rs = []
    for qq in (0.01, 0.05):
        a = kept[(~S.dec[kept]) & (q <= qq)]; ne = S.entr[a].sum()
        rs.append(f'q<={qq}: {len(a)} FDP_c {ne * (1 + 1 / S.R_ENT) / max(1, len(a)):.3f} (lower {ne / max(1, len(a)):.3f})')
    return f'{name:30s} ' + ' | '.join(rs)
sc = np.mean([S.crossfit(S.ALL, kind='mlp', seed=s) for s in (0, 10, 20, 30, 40)], axis=0)
wn = S.winners(np.arange(len(sc)), sc)
tg = wn[~S.dec[wn]]; o = tg[np.argsort(-sc[tg], kind='stable')]; fd = np.cumsum(S.entr[o]) * (1 + 1 / S.R_ENT) / np.arange(1, len(o) + 1)
print(f'run {os.environ.get("OPENDIA_RUN")}: 5-seed ensemble, entrapment-matched N at FDP_c 1%: {int(np.where(fd <= 0.01)[0].max() + 1)}', flush=True)
print(report('plain TDC', sc, wn, np.ones(len(sc))), flush=True)
w, auc = weights(sc, wn)
kept, q0 = S.tdc_q(wn, sc); t5 = sc[kept[(~S.dec[kept]) & (q0 <= 0.05)]].min(); dw = wn[S.dec[wn]]
print(report(f'PRIMARY bulk-null (AUC {auc:.3f})', sc, wn, w) + f' | tail decoy weight {w[dw][sc[dw] >= t5].mean():.2f}', flush=True)
for nm, kw in (('bulk p30', dict(pct=30)), ('bulk p70', dict(pct=70)), ('clip [0.1,10]', dict(clip=(0.1, 10))), ('logistic', dict(kind='lr'))):
    w2, a2 = weights(sc, wn, **kw)
    print(report(f'{nm} (AUC {a2:.3f})', sc, wn, w2), flush=True)
