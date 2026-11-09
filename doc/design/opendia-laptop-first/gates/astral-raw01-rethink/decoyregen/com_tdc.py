"""Label-free commonness-weighted TDC (mitigation test, declared before running).
commonness.py found that entrapment top-6 fragments overlap human-library fragments in the same window more than decoy
fragments do (AUC 0.56 on both runs). Commonness is a library property and needs no entrapment labels. Mitigation:
weight each decoy in the TDC count by w(bin) = P_target(bin) / P_decoy(bin), where bins are 10 commonness quantiles and
both P are computed over ALL library target and decoy precursors (no run data, no entrapment labels). The commonness of
a target is computed against all other library targets (human and entrapment alike, because a real pipeline cannot
tell them apart); this differs from commonness.py, which used human targets only.
q = (sum of kept-decoy weights + 1) / kept targets, monotone.
Rule: the mitigation passes if, for >=2 of 3 MLP seeds on BOTH runs, q<=0.01 FDP_c <= 0.015 and q<=0.05 FDP_c <= 0.06.
Also reported: the unweighted q for the same scores, and N at q<=0.01.
Family O, full set. Usage: OPENDIA_RUN=R02."""
import os, sys, sqlite3, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/../gap2')
import oswlike_score as S
_sel = S.select
def select(idx, score, fdr):
    pos, neg = _sel(idx, score, fdr)
    if len(pos) < 200:
        w = S.winners(idx, score); tw = w[~S.dec[w]]; pos = tw[np.argsort(-score[tw], kind='stable')[:200]]
    return pos, neg
S.select = select
HERE = os.path.dirname(os.path.abspath(__file__)); CACHE = HERE + '/commonness_all.parquet'
if not os.path.exists(CACHE):
    c = sqlite3.connect('file:' + HERE + '/merged.pqp?mode=ro', uri=True)
    tr = pd.read_sql('select m.PRECURSOR_ID pid, t.PRODUCT_MZ f, t.LIBRARY_INTENSITY li, t.ID t from TRANSITION t join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID', c)
    tr = tr.sort_values(['pid', 'li', 't'], ascending=[True, False, True]).groupby('pid').head(6)
    pr = pd.read_sql('select ID pid, PRECURSOR_MZ mz, DECOY d, TRAML_ID tid from PRECURSOR', c)
    pr['fam'] = np.where(pr.d == 0, 'T', np.where(pr.tid.str.startswith('DECOYP_'), 'P', np.where(pr.tid.str.startswith('DECOYS_'), 'S', 'O')))
    pr['win'] = np.floor((pr.mz - 380) / 2).astype(int)
    tr = tr.merge(pr[['pid', 'win', 'fam']], on='pid')
    tg = tr[tr.fam == 'T']
    ref = {w: (np.sort(g.f.values), g.pid.values[np.argsort(g.f.values)]) for w, g in tg.groupby('win')}
    out = []
    for pid, g in tr.groupby('pid'):
        a, ids = ref.get(g.win.values[0], (np.array([]), np.array([])))
        cnt = [int((ids[np.searchsorted(a, f * (1 - 17e-6)):np.searchsorted(a, f * (1 + 17e-6))] != pid).sum()) for f in g.f.values]
        out.append((pid, float(np.log1p(np.array(cnt)).mean()), g.fam.values[0]))
    pd.DataFrame(out, columns=['pid', 'com', 'fam']).to_parquet(CACHE)
cm = pd.read_parquet(CACHE); COM = dict(zip(cm.pid, cm.com))
ct, co = cm.com[cm.fam == 'T'].values, cm.com[cm.fam == 'O'].values
edges = np.unique(np.quantile(np.r_[ct, co], np.linspace(0, 1, 11)))
binof = lambda x: np.clip(np.searchsorted(edges, x, 'right') - 1, 0, len(edges) - 2)
pt = np.bincount(binof(ct), minlength=len(edges) - 1) / len(ct); po = np.bincount(binof(co), minlength=len(edges) - 1) / len(co)
WB = pt / np.maximum(po, 1e-12)
print(f'run {os.environ.get("OPENDIA_RUN")}: library commonness mean target {ct.mean():.3f} decoy O {co.mean():.3f} | bin weights {np.round(WB, 2).tolist()}', flush=True)
meta = S.m; fam = meta.family.values; keep = (fam == 'T') | (fam == 'O')
for k_, v in {k: getattr(S, k).copy() for k in ('X0', 'prec', 'dec', 'fold', 'pair', 'entr', 'right')}.items(): setattr(S, k_, v[keep])
wdec = WB[binof(np.array([COM.get(p, 0.0) for p in S.prec]))]
def tdc_w(wn, score, weighted):
    d = pd.DataFrame({'i': wn, 'p': S.pair[wn], 's': score[wn], 'd': S.dec[wn]})
    d = d.sort_values(['s', 'd'], ascending=[False, False]).drop_duplicates('p').sort_values('s', ascending=False)
    dd = d.d.values; T = np.cumsum(~dd); Dn = np.cumsum(np.where(dd, wdec[d.i.values] if weighted else 1.0, 0.0))
    q = np.minimum.accumulate(np.minimum(1, (Dn + 1) / np.maximum(T, 1))[::-1])[::-1]
    return d.i.values, q
for seed in (0, 1, 2):
    sc = S.crossfit(S.ALL, kind='mlp', seed=seed * 10)
    wn = S.winners(np.arange(len(sc)), sc); out = []
    for weighted in (False, True):
        kept, q = tdc_w(wn, sc, weighted); r = []
        for qq in (0.01, 0.05):
            a = kept[(~S.dec[kept]) & (q <= qq)]
            r.append(f'q<={qq}: {len(a)} FDP_c {S.entr[a].sum() * (1 + 1 / S.R_ENT) / max(1, len(a)):.3f}')
        out.append(('weighted ' if weighted else 'plain ') + ', '.join(r))
    print(f'seed{seed}: ' + ' | '.join(out), flush=True)
