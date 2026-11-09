"""Ranking head test (h7, declared before running): does a boosted-tree final head rank better than the MLP head?
Same semi-supervised loop as oswlike_score.crossfit (LDA rounds, then 2 final rounds), but the 2 final rounds fit a
HistGradientBoostingClassifier (positives = target winners at q<=0.02, negatives = decoy winners; class-balanced via
sample weights). Seeds 0/10/20, averaged. Compared with the 5-seed MLP ensemble on the same extraction (family O, full
set). Metric: entrapment-matched N at FDP_c 0.5/1/1.5/2% and the bulk-null weighted q (bulknull_confirm primary).
Rule: the boosted head is better if its N at FDP_c 1% exceeds the MLP ensemble's by >= 3% on BOTH runs.
Usage: OPENDIA_RUN=R02."""
import os, sys, numpy as np, pandas as pd
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
for k_, v in {k: getattr(S, k).copy() for k in ('X0', 'prec', 'dec', 'fold', 'pair', 'entr', 'right')}.items(): setattr(S, k_, v[keep])
X = np.nan_to_num(S.X0)
def crossfit_hgb(seed):
    Z = np.nan_to_num((S.X0 - np.nanmean(S.X0, 0)) / (np.nanstd(S.X0, 0) + 1e-9)); out = np.zeros(len(Z))
    for k in range(3):
        tr = np.where(S.fold != k)[0]; sc = Z[:, S.NAMES.index('xcorr_shape_w')].copy()
        for fdr in [0.15, 0.02, 0.02, 0.02, 0.02]:
            pos, neg = S.select(tr, sc, fdr); sc = Z @ S.lda_fit(Z, pos, neg)
        for it in range(2):
            pos, neg = S.select(tr, sc, 0.02)
            idx = np.r_[pos, neg]; y = np.r_[np.ones(len(pos)), np.zeros(len(neg))]
            sw = np.where(y == 1, len(neg) / max(1, len(pos)), 1.0)
            m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=31, l2_regularization=1.0,
                                               random_state=seed + it).fit(X[idx], y, sample_weight=sw)
            p = np.clip(m.predict_proba(X)[:, 1], 1e-7, 1 - 1e-7); sc = np.log(p / (1 - p))
        te = np.where(S.fold == k)[0]; wd = S.winners(te, sc); wd = wd[S.dec[wd]]
        med = np.median(sc[wd]); mad = 1.4826 * np.median(np.abs(sc[wd] - med)) + 1e-12
        out[te] = (sc[te] - med) / mad
    return out
def matched(sc):
    wn = S.winners(np.arange(len(sc)), sc); tg = wn[~S.dec[wn]]; o = tg[np.argsort(-sc[tg], kind='stable')]
    fd = np.cumsum(S.entr[o]) * (1 + 1 / S.R_ENT) / np.arange(1, len(o) + 1)
    return [int(np.where(fd <= t)[0].max() + 1) if (fd <= t).any() else 0 for t in (0.005, 0.01, 0.015, 0.02)]
mlp = np.mean([S.crossfit(S.ALL, kind='mlp', seed=s) for s in (0, 10, 20, 30, 40)], axis=0)
hgb = np.mean([crossfit_hgb(s) for s in (0, 10, 20)], axis=0)
print(f'run {os.environ.get("OPENDIA_RUN")}', flush=True)
for nm, sc in (('MLP 5-seed', mlp), ('HGB 3-seed', hgb), ('mean of both (rank)', (pd.Series(mlp).rank().values + pd.Series(hgb).rank().values))):
    S.report(nm, sc); print(f'{"":44s} matched N at FDP_c 0.5/1/1.5/2%: {matched(sc)}', flush=True)
