"""Clean-family check (post-hoc follow-up to score2, declared before running). Shuffled decoys (S) that share >=2 of their
top-6 fragments (34 ppm) with their own target borrow that target's signal (s_tail02.txt). Drop every pair whose S decoy
shares >=2 fragments, for ALL families (same targets everywhere), then score each family with LDA and the MLP at seeds 0,1,2
(seed fallback as score2_seed.py). Also: P restricted the same way.
Reading, declared before running: a family calibrates if on both runs its median-over-seeds MLP q<=0.01 FDP_c <= 0.015 and
its 5%-cut S_e/S_d band includes 1, while O on the same reduced set keeps the excess (band above 1 at 5% for >=2 of 3 seeds).
Usage: OPENDIA_RUN=R02."""
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
c = sqlite3.connect('file:merged.pqp?mode=ro', uri=True)
tr = pd.read_sql('select m.PRECURSOR_ID pid, t.PRODUCT_MZ f, t.LIBRARY_INTENSITY li, t.ID t from TRANSITION t join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID', c)
tr = tr.sort_values(['pid', 'li', 't'], ascending=[True, False, True]).groupby('pid').head(6)
F = tr.groupby('pid').f.apply(np.array).to_dict()
meta = S.m; fam = meta.family.values
u = meta[fam == 'S'][['precursor_id', 'pair']].drop_duplicates('precursor_id')
bad = set(q for p, q in zip(u.precursor_id, u.pair)
          if p in F and q in F and (np.abs(F[p][:, None] - F[q][None, :]) <= 34e-6 * F[p][:, None]).any(1).sum() >= 2)
okpair = ~np.isin(meta.pair.values, list(bad))
print(f'run {os.environ.get("OPENDIA_RUN")}: dropped pairs {len(bad)}; candidates kept {okpair.mean():.3f}', flush=True)
exec(open(os.path.dirname(os.path.abspath(__file__)) + '/score2.py').read().split("print(f'run", 1)[0].split('import oswlike_score as S', 1)[1])
full = {k: getattr(S, k).copy() for k in ('X0', 'prec', 'dec', 'fold', 'pair', 'entr', 'right')}
for Fm in ('O', 'P', 'S'):
    keep = okpair & ((fam == 'T') | (fam == Fm))
    for k, v in full.items(): setattr(S, k, v[keep])
    S.m = meta[keep].reset_index(drop=True)
    for kind, seed in (('lda', 0), ('mlp', 0), ('mlp', 1), ('mlp', 2)):
        sc = S.crossfit(S.ALL, kind=kind, seed=seed * 10)
        S.report(f'{Fm} clean {kind.upper()} seed{seed}', sc)
        print(f'{"":44s} S_e/S_d at q<=0.01/0.02/0.05 cuts: {tail(sc)}', flush=True)
