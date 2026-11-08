"""Why is decoy-estimated q anti-conservative for the OSW-like scorer (raws 01 and 02)? Design diagnostic.

Under target-decoy competition, q is valid only if a decoy winner is distributed like a FALSE target winner.
Entrapment target winners are known false targets (proxy). So compare, per subscore, entrapment target winners against
decoy winners (both at the rank-1 peak group chosen by the LDA score):
  - overall AUC (0.5 = decoys mimic false targets; > 0.5 = false targets look better than decoys -> anti-conservative)
  - the same inside the score band [decoy-winner 90th pct, 1% cut], where q<=0.01 is decided.
Also reported: the human-target vs entrapment AUC (should be large, i.e. features separate true from false).
Entrapment is used here only to DIAGNOSE; any feature choice based on it must be confirmed on an untouched run.
"""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, '.')
RUN = os.environ.get('OPENDIA_RUN', '01')
import oswlike_score as S

def auc(a, b):
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    if len(a) == 0 or len(b) == 0: return np.nan
    r = pd.Series(np.r_[a, b]).rank().values
    return (r[:len(a)].sum() - len(a) * (len(a) + 1) / 2) / (len(a) * len(b))

sc = S.crossfit(S.ALL)
w = S.winners(np.arange(len(sc)), sc)
kept, q = S.tdc_q(w, sc)
cut = sc[kept[(~S.dec[kept]) & (q <= 0.01)]].min()
dw = w[S.dec[w]]; ew = w[(~S.dec[w]) & S.entr[w]]
lo = np.quantile(sc[dw], 0.90)
band = lambda idx: idx[(sc[idx] >= lo) & (sc[idx] < cut + 2)]
print(f'raw {RUN}: decoy winners {len(dw)}, entrapment target winners {len(ew)}, 1% cut {cut:.2f}, band [{lo:.2f}, {cut + 2:.2f}]')
print(f'  score AUC entrapment vs decoy winners: all {auc(sc[ew], sc[dw]):.3f} | band {auc(sc[band(ew)], sc[band(dw)]):.3f} (n {len(band(ew))} vs {len(band(dw))})')
rows = []
for j, n in enumerate(S.NAMES):
    rows.append((n, auc(S.X0[ew, j], S.X0[dw, j]), auc(S.X0[band(ew), j], S.X0[band(dw), j])))
r = pd.DataFrame(rows, columns=['subscore', 'auc_all', 'auc_band'])
r['dev_all'] = (r.auc_all - 0.5).abs()
print(r.sort_values('dev_all', ascending=False).round(3).to_string(index=False))
r.to_csv(f'calib_auc_raw{RUN}.tsv', sep='\t', index=False)
