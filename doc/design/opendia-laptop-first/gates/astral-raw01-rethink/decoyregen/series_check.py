"""Follow-up to osw_info.py (h9b, declared before running). 'ours + series' raised N at FDP_c 1% from 1083 to 1281 on
raw 02. Two confounds: (a) the osw_hit flag (OSW reported a peak group within 5 s) was added with every group, and is
itself information from OSW's peak picker; (b) only 47% of candidates match an OSW feature.
Variants, 5-seed MLP ensemble, entrapment-matched N at FDP_c 0.5/1/1.5/2%:
  ours + hit flag only; ours + series without the flag (missing -> decoy median); ours + series + flag;
  ours + series, restricted to candidates that match an OSW feature (both arms restricted the same way: also 'ours'
  on the restricted set).
Then for 'ours + series + flag': plain and bulk-null weighted q (bulknull_confirm primary) FDP_c at q<=0.01/0.05.
Reading: series information 'carries the gap' if 'ours + series without flag' gains >= 5% over ours AND on the
restricted set ours + series gains >= 5% over ours restricted. Raw 01 replication follows when raw01_feat.osw exists.
Usage: OPENDIA_RUN=R02 (or R01)."""
import os, sys, sqlite3, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/../gap2')
import oswlike_score as S
HERE = os.path.dirname(os.path.abspath(__file__)); RUN = os.environ['OPENDIA_RUN']
_sel = S.select
def select(idx, score, fdr):
    pos, neg = _sel(idx, score, fdr)
    if len(pos) < 200:
        w = S.winners(idx, score); tw = w[~S.dec[w]]; pos = tw[np.argsort(-score[tw], kind='stable')[:200]]
    return pos, neg
S.select = select
meta = S.m; fam = meta.family.values; keep = (fam == 'T') | (fam == 'O')
FULL = {k: getattr(S, k).copy()[keep] for k in ('X0', 'prec', 'dec', 'fold', 'pair', 'entr', 'right')}
for k_, v in FULL.items(): setattr(S, k_, v)
apex = meta[keep].apex.values
c = sqlite3.connect('file:' + HERE + f'/raw{RUN[-2:]}_feat.osw?mode=ro', uri=True)
f = pd.read_sql('select f.PRECURSOR_ID pid, f.EXP_RT rt, m2.VAR_BSERIES_SCORE b, m2.VAR_YSERIES_SCORE y from FEATURE f join FEATURE_MS2 m2 on m2.FEATURE_ID=f.ID', c)
f = f.sort_values(['pid', 'rt']); fg = {p: (g.rt.values, g[['b', 'y']].values) for p, g in f.groupby('pid')}
O = np.full((len(apex), 2), np.nan); hit = np.zeros(len(apex), bool)
for i, (p, a) in enumerate(zip(S.prec, apex)):
    g = fg.get(p)
    if g is None: continue
    j = np.argmin(np.abs(g[0] - a))
    if abs(g[0][j] - a) <= 5: O[i] = g[1][j]; hit[i] = True
O = np.where(np.isnan(O), np.nanmedian(O[S.dec], axis=0), O)
print(f'run {RUN}: candidates matched {hit.mean():.3f}', flush=True)
X0, N0 = FULL['X0'], list(S.NAMES)
def ens(X, names):
    S.X0 = X; S.NAMES = names
    return np.mean([S.crossfit(names, kind='mlp', seed=s) for s in (0, 10, 20, 30, 40)], axis=0)
def matched(score):
    w = S.winners(np.arange(len(score)), score); tg = w[~S.dec[w]]; o = tg[np.argsort(-score[tg], kind='stable')]
    fd = np.cumsum(S.entr[o]) * (1 + 1 / S.R_ENT) / np.arange(1, len(o) + 1)
    return [int(np.where(fd <= t)[0].max() + 1) if (fd <= t).any() else 0 for t in (0.005, 0.01, 0.015, 0.02)]
SER = ['bseries', 'yseries']
runs = {'ours': (X0, N0), 'ours + hit flag': (np.c_[X0, hit], N0 + ['osw_hit']),
        'ours + series (no flag)': (np.c_[X0, O], N0 + SER), 'ours + series + flag': (np.c_[X0, O, hit], N0 + SER + ['osw_hit'])}
scores = {}
for nm, (X, names) in runs.items():
    scores[nm] = ens(X, names); print(f'{nm:28s} matched N {matched(scores[nm])}', flush=True)
# restricted to candidates with an OSW match
for k_, v in FULL.items(): setattr(S, k_, v[hit])
for nm, (X, names) in (('restricted: ours', (X0[hit], N0)), ('restricted: ours + series', (np.c_[X0, O][hit], N0 + SER))):
    s = ens(X, names); print(f'{nm:28s} matched N {matched(s)}', flush=True)
for k_, v in FULL.items(): setattr(S, k_, v)
# calibration of ours + series + flag
sc = scores['ours + series + flag']; S.X0 = np.c_[X0, O, hit]
wn = S.winners(np.arange(len(sc)), sc); dw = wn[S.dec[wn]]; tw = wn[~S.dec[wn]]; thr = np.median(sc[dw])
bt, bd = tw[sc[tw] < thr], dw[sc[dw] < thr]; idx = np.r_[bt, bd]; y = np.r_[np.ones(len(bt)), np.zeros(len(bd))]
X = np.nan_to_num(S.X0); pall = np.zeros(len(sc))
for k in range(3):
    trm = S.fold[idx] != k
    m = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.05, max_leaf_nodes=31, random_state=0).fit(X[idx[trm]], y[trm])
    te = dw[S.fold[dw] == k]; pall[te] = m.predict_proba(X[te])[:, 1]
r = np.clip(pall[dw] / np.maximum(1 - pall[dw], 1e-9) * len(bd) / len(bt), 0.2, 5); w = np.ones(len(sc)); w[dw] = r / r[sc[dw] < thr].mean()
for nm, wv in (('plain', np.ones(len(sc))), ('bulk-null', w)):
    d = pd.DataFrame({'i': wn, 'p': S.pair[wn], 's': sc[wn], 'd': S.dec[wn]}).sort_values(['s', 'd'], ascending=[False, False]).drop_duplicates('p').sort_values('s', ascending=False)
    dd = d.d.values; T = np.cumsum(~dd); Dn = np.cumsum(np.where(dd, wv[d.i.values], 0.0)); q = np.minimum.accumulate(np.minimum(1, (Dn + 1) / np.maximum(T, 1))[::-1])[::-1]
    kept = d.i.values; out = []
    for qq in (0.01, 0.05):
        a = kept[(~S.dec[kept]) & (q <= qq)]; out.append(f'q<={qq}: {len(a)} FDP_c {S.entr[a].sum() * (1 + 1 / S.R_ENT) / max(1, len(a)):.3f}')
    print(f'ours + series + flag, {nm} TDC: ' + ', '.join(out), flush=True)
