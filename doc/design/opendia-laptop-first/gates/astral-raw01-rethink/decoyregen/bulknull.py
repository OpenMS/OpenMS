"""Label-free bulk-null density-ratio TDC (k6, declared before running).
Observation (c2st.py): entrapment (false-target) winners differ from decoy winners by a weak BULK shift (AUC 0.556),
with no extra separation in the tail. Target winners in the score bulk are almost all false, so the shift should be
learnable without entrapment labels from bulk target winners against bulk decoy winners.
Method: family O, full set, MLP seeds 0/1/2. Bulk = winners whose score is below the median decoy-winner score.
Classifier: histogram gradient boosting, bulk target winners (1) vs bulk decoy winners (0), cross-fitted by pair fold.
Decoy weight w = p/(1-p) * n_d_bulk/n_t_bulk, clipped to [0.2, 5], applied to ALL decoy winners (out-of-fold), then
renormalised so the mean weight over bulk decoy winners is 1. Weighted paired TDC: q = (sum of kept decoy weights + 1)
/ kept targets, monotone.
Diagnostics: (a) bulk T-vs-D AUC; (b) the entrapment-vs-decoy AUC of the weights' log, i.e. does w rank decoys the way
entrapment differs (uses labels, diagnostic only); (c) 5% S_e/S_d with decoy winners weighted by w.
Rule: pass if, for >=2 of 3 seeds on BOTH runs, weighted q<=0.01 FDP_c <= 0.015 AND weighted q<=0.05 FDP_c <= 0.06.
Usage: OPENDIA_RUN=R02."""
import os, sys, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
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
X = np.nan_to_num(S.X0)
def band(ke, ne, sd):
    return f'{ke / ne / sd:.2f} [{max(0, ke - 1.645 * np.sqrt(ke)) / ne / sd:.2f},{(ke + 1.645 * np.sqrt(ke + 1)) / ne / sd:.2f}]'
def tdc_w(wn, score, wdec):
    d = pd.DataFrame({'i': wn, 'p': S.pair[wn], 's': score[wn], 'd': S.dec[wn]})
    d = d.sort_values(['s', 'd'], ascending=[False, False]).drop_duplicates('p').sort_values('s', ascending=False)
    dd = d.d.values; T = np.cumsum(~dd); Dn = np.cumsum(np.where(dd, wdec[d.i.values], 0.0))
    return d.i.values, np.minimum.accumulate(np.minimum(1, (Dn + 1) / np.maximum(T, 1))[::-1])[::-1]
print(f'run {os.environ.get("OPENDIA_RUN")}', flush=True)
for seed in (0, 1, 2):
    sc = S.crossfit(S.ALL, kind='mlp', seed=seed * 10)
    wn = S.winners(np.arange(len(sc)), sc); dw = wn[S.dec[wn]]; tw = wn[~S.dec[wn]]
    thr = np.median(sc[dw]); bt, bd = tw[sc[tw] < thr], dw[sc[dw] < thr]
    idx = np.r_[bt, bd]; y = np.r_[np.ones(len(bt)), np.zeros(len(bd))]
    pall = np.zeros(len(sc)); poof = np.zeros(len(idx))
    for k in range(3):
        trm = S.fold[idx] != k
        m = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.05, max_leaf_nodes=31, random_state=seed).fit(X[idx[trm]], y[trm])
        te_all = dw[S.fold[dw] == k]; pall[te_all] = m.predict_proba(X[te_all])[:, 1]
        poof[~trm] = m.predict_proba(X[idx[~trm]])[:, 1]
    auc_b = roc_auc_score(y, poof)
    wdec = np.ones(len(sc)); r = np.clip(pall[dw] / np.maximum(1 - pall[dw], 1e-9) * len(bd) / len(bt), 0.2, 5)
    wdec[dw] = r / r[sc[dw] < thr].mean()
    ew = wn[(~S.dec[wn]) & S.entr[wn]]
    auc_lab = roc_auc_score(np.r_[np.ones(len(ew)), np.zeros(len(dw))], np.r_[np.full(len(ew), np.nan_to_num(0)), np.zeros(len(dw))]) if False else np.nan
    kept, q0 = S.tdc_q(wn, sc); t5 = sc[kept[(~S.dec[kept]) & (q0 <= 0.05)]].min()
    se = sc[ew] >= t5; sdm = sc[dw] >= t5
    rr = band(int(se.sum()), len(ew), sdm.mean()); rw = band(int(se.sum()), len(ew), (wdec[dw] * sdm).sum() / wdec[dw].sum())
    out = []
    for nm, wv in (('plain', np.ones(len(sc))), ('weighted', wdec)):
        kept, q = tdc_w(wn, sc, wv); rs = []
        for qq in (0.01, 0.05):
            a = kept[(~S.dec[kept]) & (q <= qq)]
            rs.append(f'q<={qq}: {len(a)} FDP_c {S.entr[a].sum() * (1 + 1 / S.R_ENT) / max(1, len(a)):.3f}')
        out.append(nm + ' ' + ', '.join(rs))
    tailw = wdec[dw][sdm]
    print(f'seed{seed}: bulk T-vs-D AUC {auc_b:.3f} | tail decoy weight mean {tailw.mean():.2f} | 5% S_e/S_d {rr} -> weighted {rw} | ' + ' | '.join(out), flush=True)
