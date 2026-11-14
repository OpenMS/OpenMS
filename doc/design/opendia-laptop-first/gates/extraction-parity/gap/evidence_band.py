"""As evidence.py, but entrapment winners restricted to R's score band (R p25 .. cut), where entrapments that match real
co-eluting signal live; 200-resample bootstrap 90 % intervals of AUC. Evaluation only."""
import numpy as np
src = open('/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/gap-20261114/evidence.py').read().split('def auc(')[0]
exec(src.replace("E = w[(w.is_decoy == 0) & w.precursor_id.map(ent).eq(1) & (w.s < cut)]",
 "R0 = wr.loc[sorted(near_is_winner - accepted)]\nE = w[(w.is_decoy == 0) & w.precursor_id.map(ent).eq(1) & (w.s < cut) & (w.s >= R0.s.median())]"))
from sklearn.metrics import roc_auc_score
rng = np.random.default_rng(1)
def auc_ci(a, b, f):
    va, vb = F.loc[a, f].values, F.loc[b, f].values
    va, vb = va[np.isfinite(va)], vb[np.isfinite(vb)]
    est = roc_auc_score(np.r_[np.ones(len(va)), np.zeros(len(vb))], np.r_[va, vb])
    bs = []
    for _ in range(200):
        ia = rng.integers(0, len(va), len(va)); ib = rng.integers(0, len(vb), len(vb))
        bs.append(roc_auc_score(np.r_[np.ones(len(ia)), np.zeros(len(ib))], np.r_[va[ia], vb[ib]]))
    lo, hi = np.quantile(bs, [.05, .95]); return f'{est:.3f} ({lo:.3f}..{hi:.3f})'
print(f'R {len(R)}; entrapment winners in band [{R.s.quantile(0.25):.1f}, {cut:.1f}): {len(E)}, score p50 {E.s.median():.1f}')
print(f'{"feature":14s} {"AUC R>dec":>22s} {"AUC E>dec":>22s}')
for f in F.columns:
    print(f'{f:14s} {auc_ci(grp["R"], grp["R_dec"], f):>22s} {auc_ci(grp["E"], grp["E_dec"], f):>22s}')
