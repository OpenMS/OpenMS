"""Contamination check for the bulk-null weights (k6c, declared before running).
Concern: the bulk (winners below the median decoy score) also holds some TRUE targets; the classifier could then learn
'looks like a real peptide' and up-weight decoys in the tail for the wrong reason (over-correction).
Check (diagnostic, uses DIA-NN identity only to drop rows, never as a training label): refit the primary weights with
every target whose precursor is a DIA-NN 1% ID (dn_rt present) removed from the bulk training set. Compare with the
primary weights: Spearman correlation of decoy weights, tail mean weight, and the weighted q FDP_c.
Reading: contamination is not driving the correction if q<=0.05 FDP_c changes by < 0.01 and the decoy-weight Spearman
correlation is >= 0.8 on both runs. Also reported: the share of bulk targets that are DIA-NN IDs.
Same configuration as bulknull_confirm.py (5-seed MLP ensemble, median threshold, HGB, clip [0.2, 5])."""
import os, sys, numpy as np, pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingClassifier
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
dnid = meta.dn_rt.notna().values[keep]
for k_, v in {k: getattr(S, k).copy() for k in ('X0', 'prec', 'dec', 'fold', 'pair', 'entr', 'right')}.items(): setattr(S, k_, v[keep])
X = np.nan_to_num(S.X0)
def tdc_w(wn, score, wdec):
    d = pd.DataFrame({'i': wn, 'p': S.pair[wn], 's': score[wn], 'd': S.dec[wn]})
    d = d.sort_values(['s', 'd'], ascending=[False, False]).drop_duplicates('p').sort_values('s', ascending=False)
    dd = d.d.values; T = np.cumsum(~dd); Dn = np.cumsum(np.where(dd, wdec[d.i.values], 0.0))
    return d.i.values, np.minimum.accumulate(np.minimum(1, (Dn + 1) / np.maximum(T, 1))[::-1])[::-1]
def weights(sc, wn, drop_ids):
    dw = wn[S.dec[wn]]; tw = wn[~S.dec[wn]]; thr = np.median(sc[dw])
    bt = tw[sc[tw] < thr]
    if drop_ids: bt = bt[~dnid[bt]]
    bd = dw[sc[dw] < thr]; idx = np.r_[bt, bd]; y = np.r_[np.ones(len(bt)), np.zeros(len(bd))]
    pall = np.zeros(len(sc))
    for k in range(3):
        trm = S.fold[idx] != k
        m = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.05, max_leaf_nodes=31, random_state=0).fit(X[idx[trm]], y[trm])
        te = dw[S.fold[dw] == k]; pall[te] = m.predict_proba(X[te])[:, 1]
    r = np.clip(pall[dw] / np.maximum(1 - pall[dw], 1e-9) * len(bd) / len(bt), 0.2, 5)
    w = np.ones(len(sc)); w[dw] = r / r[sc[dw] < thr].mean()
    return w, dnid[tw[sc[tw] < thr]].mean()
def rep(sc, wn, w):
    kept, q = tdc_w(wn, sc, w); out = []
    for qq in (0.01, 0.05):
        a = kept[(~S.dec[kept]) & (q <= qq)]
        out.append(f'q<={qq}: {len(a)} FDP_c {S.entr[a].sum() * (1 + 1 / S.R_ENT) / max(1, len(a)):.3f}')
    return ', '.join(out)
sc = np.mean([S.crossfit(S.ALL, kind='mlp', seed=s) for s in (0, 10, 20, 30, 40)], axis=0)
wn = S.winners(np.arange(len(sc)), sc); dw = wn[S.dec[wn]]
kept, q0 = S.tdc_q(wn, sc); t5 = sc[kept[(~S.dec[kept]) & (q0 <= 0.05)]].min(); tail = sc[dw] >= t5
w1, share = weights(sc, wn, False); w2, _ = weights(sc, wn, True)
print(f'run {os.environ.get("OPENDIA_RUN")}: bulk targets that are DIA-NN IDs {share:.3f}', flush=True)
print(f'  primary      {rep(sc, wn, w1)} | tail decoy weight {w1[dw][tail].mean():.2f}', flush=True)
print(f'  IDs dropped  {rep(sc, wn, w2)} | tail decoy weight {w2[dw][tail].mean():.2f} | Spearman(decoy weights) {spearmanr(w1[dw], w2[dw]).correlation:.3f}', flush=True)
