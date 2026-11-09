"""Why do OSW's b/y series scores add information and ours not? (h11b, diagnostic, raw 02.)
On candidates matched to an OSW feature (same precursor, apex within 5 s): Spearman correlation of our by_nb/by_ny with
OSW var_bseries/var_yseries; and the conditional AUC (target vs decoy winners, within 20 bins of our 5-seed MLP score,
top 5% of decoy winners and above) of each, computed on the SAME matched winners. Also the conditional AUC of
'OSW-matched' itself (whether OSW reported a peak group there) within those bins, on all winners.
Usage: OPENDIA_RUN=R02."""
import os, sys, sqlite3, numpy as np, pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/../gap2')
import oswlike_score as S
HERE = os.path.dirname(os.path.abspath(__file__))
_sel = S.select
def select(idx, score, fdr):
    pos, neg = _sel(idx, score, fdr)
    if len(pos) < 200:
        w = S.winners(idx, score); tw = w[~S.dec[w]]; pos = tw[np.argsort(-score[tw], kind='stable')[:200]]
    return pos, neg
S.select = select
BY = np.load(os.environ['JCODE_SCRATCH_DIR'] + '/xicR02/by_feats.npy')
meta = S.m; fam = meta.family.values; keep = (fam == 'T') | (fam == 'O')
for k_, v in {k: getattr(S, k).copy() for k in ('X0', 'prec', 'dec', 'fold', 'pair', 'entr', 'right')}.items(): setattr(S, k_, v[keep])
BY = BY[keep]; apex = meta[keep].apex.values
c = sqlite3.connect('file:' + HERE + '/raw02_feat.osw?mode=ro', uri=True)
f = pd.read_sql('select f.PRECURSOR_ID pid, f.EXP_RT rt, m2.VAR_BSERIES_SCORE b, m2.VAR_YSERIES_SCORE y from FEATURE f join FEATURE_MS2 m2 on m2.FEATURE_ID=f.ID', c)
f = f.sort_values(['pid', 'rt']); fg = {p: (g.rt.values, g[['b', 'y']].values) for p, g in f.groupby('pid')}
O = np.full((len(apex), 2), np.nan); hit = np.zeros(len(apex), bool)
for i, (p, a) in enumerate(zip(S.prec, apex)):
    g = fg.get(p)
    if g is None: continue
    j = np.argmin(np.abs(g[0] - a))
    if abs(g[0][j] - a) <= 5: O[i] = g[1][j]; hit[i] = True
print(f'OSW series value range: b {np.nanmin(O[:, 0])}..{np.nanmax(O[:, 0])}, integer share {np.mean(np.nan_to_num(O[hit, 0]) == np.round(np.nan_to_num(O[hit, 0]))):.3f}', flush=True)
print(f'Spearman (matched candidates): ours by_nb vs OSW b {spearmanr(BY[hit, 0], O[hit, 0]).correlation:.3f}, ours by_ny vs OSW y {spearmanr(BY[hit, 1], O[hit, 1]).correlation:.3f}', flush=True)
sc = np.mean([S.crossfit(S.ALL, kind='mlp', seed=s) for s in (0, 10, 20, 30, 40)], axis=0)
wn = S.winners(np.arange(len(sc)), sc); dw = wn[S.dec[wn]]
lo = np.quantile(sc[dw], 0.95)
def cauc(sel, v):
    edges = np.quantile(sc[sel], np.linspace(0, 1, 21)); b = np.clip(np.searchsorted(edges, sc[sel], 'right') - 1, 0, 19); num = den = 0.0
    for k in range(20):
        s = sel[b == k]; y = ~S.dec[s]
        if y.sum() >= 10 and (~y).sum() >= 10:
            a = roc_auc_score(y, v[s]); num += len(s) * a; den += len(s)
    return num / den if den else np.nan
selm = wn[(sc[wn] >= lo) & hit[wn]]; sela = wn[sc[wn] >= lo]
print(f'winners above the cut: {len(sela)}, OSW-matched {len(selm)} ({len(selm) / len(sela):.3f}); matched share targets {hit[sela][~S.dec[sela]].mean():.3f} decoys {hit[sela][S.dec[sela]].mean():.3f}', flush=True)
for nm, v in (('ours by_nb', BY[:, 0]), ('ours by_ny', BY[:, 1]), ('ours by_run', BY[:, 2]), ('ours by_int', BY[:, 3]), ('OSW bseries', O[:, 0]), ('OSW yseries', O[:, 1])):
    print(f'  conditional AUC on matched winners: {nm:12s} {cauc(selm, np.nan_to_num(v)):.3f}', flush=True)
print(f'  conditional AUC on all winners: OSW-matched flag {cauc(sela, hit.astype(float)):.3f}, ours by_ny {cauc(sela, BY[:, 1]):.3f}', flush=True)
