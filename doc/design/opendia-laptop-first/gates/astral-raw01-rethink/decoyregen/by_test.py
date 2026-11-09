"""Test of our own b/y features (h11, declared before running). Add by_feats.npy (by_nb, by_ny, by_run, by_int) to the
20 subscores; 5-seed MLP ensemble, family O, full set. Report entrapment-matched N at FDP_c 0.5/1/1.5/2% for ours and
ours + by, and plain / bulk-null weighted q FDP_c at q<=0.01 / 0.05 for ours + by.
Rule: pass if ours + by raises N at FDP_c 1% by >= 5% on BOTH runs, AND bulk-null q<=0.05 FDP_c <= 0.06 on both.
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
BYSET = os.environ.get('BY_SET', 'feats')  # 'feats' (v1 fractions) or 'counts' (v2 integer counts, declared after by_vs_osw02.txt)
BY = np.load(os.environ['JCODE_SCRATCH_DIR'] + '/xic' + os.environ['OPENDIA_RUN'] + ('/by_feats.npy' if BYSET == 'feats' else '/by_counts.npy'))
meta = S.m; fam = meta.family.values; keep = (fam == 'T') | (fam == 'O')
for k_, v in {k: getattr(S, k).copy() for k in ('X0', 'prec', 'dec', 'fold', 'pair', 'entr', 'right')}.items(): setattr(S, k_, v[keep])
BY = BY[keep]; X0, N0 = S.X0.copy(), list(S.NAMES)
def ens(X, names):
    S.X0 = X; S.NAMES = names
    return np.mean([S.crossfit(names, kind='mlp', seed=s) for s in (0, 10, 20, 30, 40)], axis=0)
def matched(score):
    w = S.winners(np.arange(len(score)), score); tg = w[~S.dec[w]]; o = tg[np.argsort(-score[tg], kind='stable')]
    fd = np.cumsum(S.entr[o]) * (1 + 1 / S.R_ENT) / np.arange(1, len(o) + 1)
    return [int(np.where(fd <= t)[0].max() + 1) if (fd <= t).any() else 0 for t in (0.005, 0.01, 0.015, 0.02)]
def qreport(sc, X):
    wn = S.winners(np.arange(len(sc)), sc); dw = wn[S.dec[wn]]; tw = wn[~S.dec[wn]]; thr = np.median(sc[dw])
    bt, bd = tw[sc[tw] < thr], dw[sc[dw] < thr]; idx = np.r_[bt, bd]; y = np.r_[np.ones(len(bt)), np.zeros(len(bd))]
    Xn = np.nan_to_num(X); pall = np.zeros(len(sc))
    for k in range(3):
        trm = S.fold[idx] != k
        mdl = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.05, max_leaf_nodes=31, random_state=0).fit(Xn[idx[trm]], y[trm])
        te = dw[S.fold[dw] == k]; pall[te] = mdl.predict_proba(Xn[te])[:, 1]
    r = np.clip(pall[dw] / np.maximum(1 - pall[dw], 1e-9) * len(bd) / len(bt), 0.2, 5); w = np.ones(len(sc)); w[dw] = r / r[sc[dw] < thr].mean()
    out = []
    for nm, wv in (('plain', np.ones(len(sc))), ('bulk-null', w)):
        d = pd.DataFrame({'i': wn, 'p': S.pair[wn], 's': sc[wn], 'd': S.dec[wn]}).sort_values(['s', 'd'], ascending=[False, False]).drop_duplicates('p').sort_values('s', ascending=False)
        dd = d.d.values; T = np.cumsum(~dd); Dn = np.cumsum(np.where(dd, wv[d.i.values], 0.0)); q = np.minimum.accumulate(np.minimum(1, (Dn + 1) / np.maximum(T, 1))[::-1])[::-1]
        kept = d.i.values; rs = []
        for qq in (0.01, 0.05):
            a = kept[(~S.dec[kept]) & (q <= qq)]; rs.append(f'q<={qq}: {len(a)} FDP_c {S.entr[a].sum() * (1 + 1 / S.R_ENT) / max(1, len(a)):.3f}')
        out.append(f'{nm} ' + ', '.join(rs))
    return ' | '.join(out)
print(f'run {os.environ["OPENDIA_RUN"]}', flush=True)
s0 = ens(X0, N0); print(f'ours            matched N {matched(s0)}', flush=True)
BN = ['by_nb', 'by_ny', 'by_run', 'by_int'] if BYSET == 'feats' else ['by_count_b', 'by_count_y']
s1 = ens(np.c_[X0, BY], N0 + BN); print(f'ours + own b/y  matched N {matched(s1)}', flush=True)
print(f'ours + own b/y  {qreport(s1, np.c_[X0, BY])}', flush=True)
