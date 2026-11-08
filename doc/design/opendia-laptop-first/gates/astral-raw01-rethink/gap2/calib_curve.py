"""Is the raw 02 anti-conservativeness a calibration bias or sampling noise? Design diagnostic.

For the LDA score, and the MLP (seed 0) when OPENDIA_MLP=1: along the TDC-kept target list, compare decoy-estimated
FDR ((D+1)/T) with entrapment FDP_c, reporting the entrapment count and an exact Poisson 90% interval for FDP_c. A
calibration bias shows as FDP_c above the decoy estimate across many thresholds. Noise shows as intervals that
contain the decoy estimate.
"""
import os, sys, numpy as np
from scipy.stats import chi2
sys.path.insert(0, '.')
RUN = os.environ.get('OPENDIA_RUN', '01')
import oswlike_score as S

def curve(name, sc):
    w = S.winners(np.arange(len(sc)), sc); kept, q = S.tdc_q(w, sc)
    d = S.dec[kept]; T = np.cumsum(~d); D = np.cumsum(d)
    tgt = kept[~d]; e = np.cumsum(S.entr[tgt])
    f = (1 + 1 / S.R_ENT)
    print(f'raw {RUN} {name}')
    for qq in (0.005, 0.01, 0.02, 0.03, 0.05):
        n = int(((~d) & (q <= qq)).sum())
        if n == 0: continue
        k = int(e[n - 1]); est = (D[np.searchsorted(T, n)] + 1) / n
        lo = chi2.ppf(0.05, 2 * k) / 2 if k else 0.0; hi = chi2.ppf(0.95, 2 * k + 2) / 2
        print(f'  q<={qq:<5} N {n:5d}  decoy FDR {est:.4f}  entrapment {k:3d}  FDP_c {k * f / n:.4f}  90% CI [{lo * f / n:.4f}, {hi * f / n:.4f}]')

curve('LDA', S.crossfit(S.ALL))
if os.environ.get('OPENDIA_MLP') == '1':
    curve('MLP seed 0', S.crossfit(S.ALL, kind='mlp', seed=0))
