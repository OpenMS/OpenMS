"""Does family S calibrate because its decoys borrow their own target's signal? Shuffled decoys sit at the target's
precursor m/z (97.8%) and share more top-6 fragments with it (family_overlap.txt). Score S with the MLP (seed fallback as
score2_seed.py), then compare the share of decoys sharing >=2 / >=3 fragments with their target among the decoy winners
above each q cut against all decoy winners. Repeat for P and O. Post-hoc diagnostic, no pass rule."""
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
# shared top-6 fragments between each decoy and its target (merged PQP IDs)
c = sqlite3.connect('file:merged.pqp?mode=ro', uri=True)
tr = pd.read_sql('select m.PRECURSOR_ID pid, t.PRODUCT_MZ f, t.LIBRARY_INTENSITY li, t.ID t from TRANSITION t join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID', c)
tr = tr.sort_values(['pid', 'li', 't'], ascending=[True, False, True]).groupby('pid').head(6)
F = tr.groupby('pid').f.apply(np.array).to_dict()
meta = S.m; dpid = meta.precursor_id.values; dpair = meta.pair.values
u = pd.DataFrame({'p': dpid, 'q': dpair, 'd': S.dec}).drop_duplicates('p'); u = u[u.d]
sh = {}
for p, q in zip(u.p, u.q):
    a, b = F.get(p), F.get(q)
    sh[p] = -1 if a is None or b is None else int((np.abs(a[:, None] - b[None, :]) <= 34e-6 * a[:, None]).any(1).sum())
full = {k: getattr(S, k).copy() for k in ('X0', 'prec', 'dec', 'fold', 'pair', 'entr', 'right')}
fam = meta.family.values
for Fm in ('S', 'P', 'O'):
    keep = (fam == 'T') | (fam == Fm)
    for k, v in full.items(): setattr(S, k, v[keep])
    S.m = meta[keep].reset_index(drop=True)
    sc = S.crossfit(S.ALL, kind='mlp')
    w = S.winners(np.arange(len(sc)), sc); kept, q = S.tdc_q(w, sc); dw = w[S.dec[w]]
    shw = np.array([sh[p] for p in S.prec[dw]])
    line = [f'all decoy winners n {len(dw)}: >=2 {(shw >= 2).mean():.3f}, >=3 {(shw >= 3).mean():.3f}']
    for qq in (0.01, 0.02, 0.05, 0.10):
        sel = kept[(~S.dec[kept]) & (q <= qq)]; t = sc[sel].min(); a = shw[sc[dw] >= t]
        line.append(f'above q<={qq} cut n {len(a)}: >=2 {(a >= 2).mean():.3f}, >=3 {(a >= 3).mean():.3f}')
    print(Fm, ' | '.join(line), flush=True)
