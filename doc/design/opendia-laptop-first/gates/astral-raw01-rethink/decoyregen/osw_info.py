"""Which information do our 20 subscores lack? (h9, diagnostic; OSW subscores are never used in a shipped model.)
Input: raw02_feat.osw (OSW scores of the same extraction, up to 5 peak groups per precursor) and xicR02.
Each of our candidates is matched to the OSW feature of the same precursor nearest in RT (within 5 s); unmatched
candidates get the column's decoy median plus a missing flag.
(1) Conditional separation: for each OSW subscore, target vs decoy winners (our 5-seed MLP winners), the AUC within
    20 bins of our score over the top 5% of decoy winners and above, averaged with bin weights ('cAUC'). A value near
    0.5 means no information beyond our score.
(2) Augmentation (declared before running): add OSW subscore groups to our 20 features and rerun the 5-seed MLP
    ensemble; report entrapment-matched N at FDP_c 0.5/1/1.5/2%. Groups: massdev (MS2 + MS1), mi (MS2 MI scores + MS1
    MI), series (b/y series), elution (elution model fit), iso (isotope correlation/overlap MS2 + MS1), all.
    Reading: a group 'carries the gap' if it raises N at FDP_c 1% by >= 5% (>= 54) over our 20 features.
Usage: OPENDIA_RUN=R02 (only R02 has an OSW feature table)."""
import os, sys, sqlite3, numpy as np, pandas as pd
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
meta = S.m; fam = meta.family.values; keep = (fam == 'T') | (fam == 'O')
for k_, v in {k: getattr(S, k).copy() for k in ('X0', 'prec', 'dec', 'fold', 'pair', 'entr', 'right')}.items(): setattr(S, k_, v[keep])
m = meta[keep].reset_index(drop=True); apex = m.apex.values
c = sqlite3.connect('file:' + HERE + '/raw02_feat.osw?mode=ro', uri=True)
f = pd.read_sql('select f.PRECURSOR_ID pid, f.EXP_RT rt, m2.*, m1.VAR_MASSDEV_SCORE ms1_massdev, m1.VAR_MI_SCORE ms1_mi, m1.VAR_MI_COMBINED_SCORE ms1_mi_comb, '
                'm1.VAR_ISOTOPE_CORRELATION_SCORE ms1_isocorr, m1.VAR_ISOTOPE_OVERLAP_SCORE ms1_isoovl, m1.VAR_XCORR_SHAPE_COMBINED ms1_xshape_comb, m1.VAR_XCORR_COELUTION_COMBINED ms1_xcoel_comb '
                'from FEATURE f join FEATURE_MS2 m2 on m2.FEATURE_ID=f.ID left join FEATURE_MS1 m1 on m1.FEATURE_ID=f.ID', c)
f.columns = [x.lower() for x in f.columns]
cols = [x for x in f.columns if x.startswith('var_') or x.startswith('ms1_')]
cols = [x for x in cols if f[x].notna().mean() > 0.5 and f[x].std() > 0]
print(f'OSW features {len(f)}, precursors {f.pid.nunique()}, usable columns {len(cols)}', flush=True)
# nearest OSW feature within 5 s, per candidate
f = f.sort_values(['pid', 'rt']); fg = {p: (g.rt.values, g[cols].values) for p, g in f.groupby('pid')}
O = np.full((len(apex), len(cols)), np.nan); hit = np.zeros(len(apex), bool)
for i, (p, a) in enumerate(zip(S.prec, apex)):
    g = fg.get(p)
    if g is None: continue
    j = np.argmin(np.abs(g[0] - a))
    if abs(g[0][j] - a) <= 5: O[i] = g[1][j]; hit[i] = True
print(f'candidates matched to an OSW feature: {hit.mean():.3f} (targets {hit[~S.dec].mean():.3f}, decoys {hit[S.dec].mean():.3f})', flush=True)
dmed = np.nanmedian(O[S.dec], axis=0); O = np.where(np.isnan(O), dmed, O)
sc = np.mean([S.crossfit(S.ALL, kind='mlp', seed=s) for s in (0, 10, 20, 30, 40)], axis=0)
wn = S.winners(np.arange(len(sc)), sc); dw = wn[S.dec[wn]]; tw = wn[~S.dec[wn]]
lo = np.quantile(sc[dw], 0.95); sel = wn[sc[wn] >= lo]
edges = np.quantile(sc[sel], np.linspace(0, 1, 21)); b = np.clip(np.searchsorted(edges, sc[sel], 'right') - 1, 0, 19)
res = []
for j, n in enumerate(cols):
    num = den = 0.0
    for k in range(20):
        s = sel[b == k]; y = ~S.dec[s]
        if y.sum() >= 10 and (~y).sum() >= 10:
            a = roc_auc_score(y, O[s, j]); num += len(s) * max(a, 1 - a); den += len(s)
    res.append((num / den if den else np.nan, n))
print('conditional AUC beyond our score (top 5% of decoy winners and above; max(a, 1-a)):', flush=True)
for a, n in sorted(res, reverse=True)[:15]: print(f'  {n:34s} {a:.3f}', flush=True)
def matched(score):
    w = S.winners(np.arange(len(score)), score); tg = w[~S.dec[w]]; o = tg[np.argsort(-score[tg], kind='stable')]
    fd = np.cumsum(S.entr[o]) * (1 + 1 / S.R_ENT) / np.arange(1, len(o) + 1)
    return [int(np.where(fd <= t)[0].max() + 1) if (fd <= t).any() else 0 for t in (0.005, 0.01, 0.015, 0.02)]
G = {'massdev': [x for x in cols if 'massdev' in x], 'mi': [x for x in cols if '_mi' in x or 'mi_' in x.replace('var_', '', 1)[:3]],
     'series': [x for x in cols if 'series' in x], 'elution': [x for x in cols if 'elution_model' in x],
     'iso': [x for x in cols if 'isotope' in x or 'iso' in x.split('_')[-1] or 'isocorr' in x or 'isoovl' in x]}
G['all'] = cols
print('groups: ' + '; '.join(f'{k}: {v}' for k, v in G.items() if k != 'all'), flush=True)
X0, N0 = S.X0.copy(), list(S.NAMES)
print(f'ours (20 subscores): matched N {matched(sc)}', flush=True)
for gname, gc in G.items():
    if not gc: continue
    idx = [cols.index(x) for x in gc]
    S.X0 = np.c_[X0, O[:, idx], hit.astype(float)]; S.NAMES = N0 + gc + ['osw_hit']
    s2 = np.mean([S.crossfit(S.NAMES, kind='mlp', seed=s) for s in (0, 10, 20, 30, 40)], axis=0)
    print(f'ours + {gname} ({len(gc)} cols): matched N {matched(s2)}', flush=True)
S.X0, S.NAMES = X0, N0
