"""Predeclared mitigation test for the MLP's tail excess (design diagnostic; raws 01 and 02).

Variants (same features and SS loop, out-of-fold):
  A  LDA (reference, from the cache)
  B  MLP seed 0 (from the cache)
  C  MLP 5-seed ensemble: mean of fold-normalised scores, seeds 0..4
  D  rank average of A and C
Pass rule, set before running: q<=0.01 FDP_c <= 0.015 on BOTH runs, and N at FDP_c 1% within 3% of B or better.
Also reported: S_e/S_d at the q<=0.01 / 0.02 / 0.05 cuts (see calib_tail.py).
"""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, '.')
RUN = os.environ.get('OPENDIA_RUN', '01')
import oswlike_score as S
z = np.load(f'{S.S}/calib_scores.npz'); lda, mlp0 = z['lda'], z['mlp']
ens_path = f'{S.S}/calib_ens.npy'
if os.path.exists(ens_path):
    ens = np.load(ens_path)
else:
    ens = np.mean([mlp0] + [S.crossfit(S.ALL, kind='mlp', seed=s) for s in (1, 2, 3, 4)], axis=0); np.save(ens_path, ens)
rk = lambda x: pd.Series(x).rank(pct=True).values
mix = rk(lda) + rk(ens)

def ratio(sc):
    w = S.winners(np.arange(len(sc)), sc); kept, q = S.tdc_q(w, sc)
    dw = w[S.dec[w]]; ew = w[(~S.dec[w]) & S.entr[w]]; out = []
    for qq in (0.01, 0.02, 0.05):
        t = sc[kept[(~S.dec[kept]) & (q <= qq)]].min()
        out.append(f'{((sc[ew] >= t).sum() / len(ew)) / max(1e-12, (sc[dw] >= t).sum() / len(dw)):.2f}')
    return '/'.join(out)

print(f'raw {RUN}')
for name, sc in (('A LDA', lda), ('B MLP seed 0', mlp0), ('C MLP 5-seed ensemble', ens), ('D rank avg LDA+ens', mix)):
    S.report(name, sc)
    print(f'{"":44s} S_e/S_d at q<=0.01/0.02/0.05 cuts: {ratio(sc)}', flush=True)
