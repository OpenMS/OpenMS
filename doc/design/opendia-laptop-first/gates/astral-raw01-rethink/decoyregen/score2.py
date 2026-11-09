"""Decoy-family test (k3). On one extraction of targets plus three decoy families, run the same teacher-free scorer
separately with each family as the decoy set, and compare decoy-estimated q with entrapment FDP_c, plus the tail ratio
S_e/S_d. Families: O = original pseudo-reverse with switchKR (control, the same decoys as the library), P = pseudo-reverse
without switchKR, S = shuffle without switchKR.
Predeclared reading: a family fixes calibration if, on both runs, its q<=0.01 FDP_c <= 0.015 AND its S_e/S_d at the 5%
cut lies within the 90% band of 1, while O on the same extraction reproduces the excess. If O does not reproduce the
excess, the extraction differs from the reference and the test is inconclusive.
Usage: OPENDIA_RUN=R02 (reads $JCODE_SCRATCH_DIR/xicR02)."""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/../gap2')
import oswlike_score as S
full = {k: getattr(S, k).copy() for k in ('X0', 'prec', 'dec', 'fold', 'pair', 'entr', 'right')}
fam = S.m.family.values; M0 = S.m

def tail(sc):
    w = S.winners(np.arange(len(sc)), sc); kept, q = S.tdc_q(w, sc)
    dw = w[S.dec[w]]; ew = w[(~S.dec[w]) & S.entr[w]]; out = []
    for qq in (0.01, 0.02, 0.05):
        sel = kept[(~S.dec[kept]) & (q <= qq)]
        if len(sel) == 0: out.append('n/a'); continue
        t = sc[sel].min(); ke = int((sc[ew] >= t).sum()); kd = int((sc[dw] >= t).sum())
        r = (ke / len(ew)) / max(1e-12, kd / len(dw))
        lo = max(0, ke - 1.645 * np.sqrt(ke)) / len(ew) / max(1e-12, kd / len(dw)); hi = (ke + 1.645 * np.sqrt(ke + 1)) / len(ew) / max(1e-12, kd / len(dw))
        out.append(f'{r:.2f} [{lo:.2f},{hi:.2f}] ({ke}/{kd})')
    return ' | '.join(out)

print(f'run {os.environ.get("OPENDIA_RUN")}: candidates by family', pd.Series(fam).value_counts().to_dict(), flush=True)
for F, label in (('O', 'O original pseudo-reverse+switchKR (control)'), ('P', 'P pseudo-reverse, no switchKR'), ('S', 'S shuffle, no switchKR')):
    keep = (fam == 'T') | (fam == F)
    for k, v in full.items():
        setattr(S, k, v[keep])
    S.m = M0[keep].reset_index(drop=True)
    for kind in ('lda', 'mlp'):
        sc = S.crossfit(S.ALL, kind=kind)
        S.report(f'{label} {kind.upper()}', sc)
        print(f'{"":44s} S_e/S_d at q<=0.01/0.02/0.05 cuts: {tail(sc)}', flush=True)
