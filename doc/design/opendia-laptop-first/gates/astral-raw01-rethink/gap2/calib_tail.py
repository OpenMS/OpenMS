"""Tail test: do decoy winners mimic false (entrapment) target winners where it matters? Design diagnostic.

Before competition, at a score threshold t: S_e(t) = share of entrapment target winners above t, and S_d(t) = share of
decoy winners above t. Decoys are pseudo-reversed from both human and entrapment targets. Valid TDC needs
S_e/S_d ~ 1 (a little below, since some entrapment winners may sit on real but foreign signal). Above 1, false
targets beat decoys, so q is anti-conservative.
The ratio is reported at the score thresholds where TDC q = 0.5/1/2/5%, and split by the subscore group that most
often carries the excess: entrapment winners that are 'library-like' (lib_dotprod in its top 10%) versus the rest.
"""
import os, sys, numpy as np
sys.path.insert(0, '.')
RUN = os.environ.get('OPENDIA_RUN', '01')
import oswlike_score as S

def tail(name, sc):
    w = S.winners(np.arange(len(sc)), sc); kept, q = S.tdc_q(w, sc)
    dw = w[S.dec[w]]; ew = w[(~S.dec[w]) & S.entr[w]]
    print(f'raw {RUN} {name}: decoy winners {len(dw)}, entrapment winners {len(ew)}')
    for qq in (0.005, 0.01, 0.02, 0.05):
        tk = kept[(~S.dec[kept]) & (q <= qq)]
        t = sc[tk].min()
        ke = int((sc[ew] >= t).sum()); kd = int((sc[dw] >= t).sum())
        ratio = (ke / len(ew)) / max(1e-12, kd / len(dw))
        # Poisson 90% band on the ratio, treating kd as known
        lo = max(0, ke - 1.645 * np.sqrt(ke)) / len(ew) / max(1e-12, kd / len(dw))
        hi = (ke + 1.645 * np.sqrt(ke + 1)) / len(ew) / max(1e-12, kd / len(dw))
        print(f'  t(q<={qq:<5})={t:7.2f}  entrapment above {ke:4d}  decoys above {kd:4d}  S_e/S_d {ratio:.2f} [{lo:.2f}, {hi:.2f}]')

tail('LDA', S.crossfit(S.ALL))
tail('MLP seed 0', S.crossfit(S.ALL, kind='mlp', seed=0))
