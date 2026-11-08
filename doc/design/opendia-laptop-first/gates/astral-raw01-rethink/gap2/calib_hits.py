"""Subscore medians of entrapment against decoy winners above the q<=0.05 cut (cached LDA/MLP scores), with a pooled
Mann-Whitney p per subscore. This is a small-n design diagnostic. It looks for the feature where false targets look better than decoys."""
import os, sys, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
sys.path.insert(0, '.')
RUN = os.environ.get('OPENDIA_RUN', '01')
import oswlike_score as S
z = np.load(f'{S.S}/calib_scores.npz')
for name in ('lda', 'mlp'):
    sc = z[name]; w = S.winners(np.arange(len(sc)), sc); kept, q = S.tdc_q(w, sc)
    t = sc[kept[(~S.dec[kept]) & (q <= 0.05)]].min()
    e = w[(~S.dec[w]) & S.entr[w] & (sc[w] >= t)]; d = w[S.dec[w] & (sc[w] >= t)]
    rows = []
    for j, n in enumerate(S.NAMES):
        a, b = S.X0[e, j], S.X0[d, j]; a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
        rows.append((n, np.median(a), np.median(b), mannwhitneyu(a, b).pvalue))
    r = pd.DataFrame(rows, columns=['subscore', 'entr_med', 'decoy_med', 'p']).sort_values('p')
    print(f'raw {RUN} {name.upper()} above q<=0.05 cut: entrapment {len(e)}, decoys {len(d)}')
    print(r.head(6).round(4).to_string(index=False))
