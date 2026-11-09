"""Classifier two-sample test (C2ST): do entrapment winners and decoy winners differ in the measured subscores?
Family O, full set. Winners by MLP seed 0 (seed fallback as score2_seed.py). Classifier: histogram gradient boosting,
entrapment winner (1) vs decoy winner (0), 3-fold cross-fitted by pair-hash fold, on the 20 subscores. Reported:
  (1) out-of-fold AUC (all winners), and per subscore the AUC drop when that column is permuted;
  (2) 5% S_e/S_d after reweighting decoy winners by the out-of-fold density ratio p/(1-p) * (n_d/n_e), clipped at 20.
Reading, declared before running: the mismatch is 'in the measured features' on a run if (2) brings the 5% band to
include 1 (by construction this should hold if the classifier captures the difference, so the main information is
(1): AUC and which features carry it). A C2ST AUC below 0.52 means the winners are not separable by these features,
so the tail excess would come from a difference too small or too localised for the classifier.
Diagnostic only: entrapment labels are never used for a production model."""
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
S.m = meta[keep].reset_index(drop=True)
sc = S.crossfit(S.ALL, kind='mlp', seed=0)
w = S.winners(np.arange(len(sc)), sc); kept, q = S.tdc_q(w, sc)
t = sc[kept[(~S.dec[kept]) & (q <= 0.05)]].min()
dw = w[S.dec[w]]; ew = w[(~S.dec[w]) & S.entr[w]]
idx = np.r_[ew, dw]; y = np.r_[np.ones(len(ew)), np.zeros(len(dw))]
X = S.X0[idx]; fo = S.fold[idx]; p = np.zeros(len(idx)); models = []
for k in range(3):
    tr, te = fo != k, fo == k
    m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=31, random_state=0).fit(X[tr], y[tr])
    p[te] = m.predict_proba(X[te])[:, 1]; models.append((m, te))
auc = roc_auc_score(y, p)
rng = np.random.default_rng(0); drops = []
for j, n in enumerate(S.NAMES):
    pp = np.zeros(len(idx))
    for m, te in models:
        Xt = X[te].copy(); Xt[:, j] = rng.permutation(Xt[:, j]); pp[te] = m.predict_proba(Xt)[:, 1]
    drops.append((auc - roc_auc_score(y, pp), n))
print(f'run {os.environ.get("OPENDIA_RUN")}: C2ST AUC entrapment vs decoy winners {auc:.4f} (n e {len(ew)}, d {len(dw)})', flush=True)
print('  permutation AUC drop: ' + ', '.join(f'{n} {d:.4f}' for d, n in sorted(drops, reverse=True)[:8]), flush=True)
pd_ = p[len(ew):]; wts = np.clip(pd_ / np.maximum(1 - pd_, 1e-9) * len(dw) / len(ew), 0, 20)
se = (sc[ew] >= t); sdr = (sc[dw] >= t).mean(); sdw = (wts * (sc[dw] >= t)).sum() / wts.sum()
def band(ke, ne, sd):
    return f'{ke / ne / sd:.2f} [{max(0, ke - 1.645 * np.sqrt(ke)) / ne / sd:.2f},{(ke + 1.645 * np.sqrt(ke + 1)) / ne / sd:.2f}]'
print(f'  5% S_e/S_d raw {band(int(se.sum()), len(ew), sdr)} -> density-ratio reweighted {band(int(se.sum()), len(ew), sdw)} | weight max {wts.max():.1f}, ESS {wts.sum() ** 2 / (wts ** 2).sum():.0f}', flush=True)
# tail-local view: classifier score among winners above a looser cut (q<=0.20)
t2 = sc[kept[(~S.dec[kept]) & (q <= 0.20)]].min(); loc = np.r_[sc[ew] >= t2, sc[dw] >= t2]
if loc.sum() > 50 and y[loc].min() < y[loc].max():
    print(f'  C2ST AUC within winners above the q<=0.20 cut: {roc_auc_score(y[loc], p[loc]):.3f} (n e {int(y[loc].sum())}, d {int((1 - y[loc]).sum())})', flush=True)
