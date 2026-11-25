#!/usr/bin/env python3
"""Offline per-precursor score combinations of the network and GBM calibrated scores (both -log decoy survival per fold).
Each precursor: n = network kept-window score, g = GBM kept-window score. Combined c = f(n, g), then the tool's pooled
paired TDC (per-fold decoy median/MAD standardisation, ties to decoy) and entrapment metrics. Diagnosis only."""
import sqlite3, sys, numpy as np, pandas as pd
lib = sys.argv[1]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID precursor_id, p.TRAML_ID t, p.DECOY d,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
pq['key'] = np.where(pq.d == 1, pq.t.str.replace('^DECOY_', '', regex=True), pq.t)
P = pq.set_index('precursor_id')
COMBOS = {
    'net': lambda n, g: n,
    'gbm': lambda n, g: g,
    'mean': lambda n, g: 0.5 * (n + g),
    'min': lambda n, g: np.minimum(n, g),
    'w.25': lambda n, g: 0.25 * n + 0.75 * g,
}
def calibrate_tail(x, col):
    """Replica of CrossFitScorer::calibrateTail on one kept score per precursor (per fold, decoy winners)."""
    out = np.empty(len(x))
    for fo in np.unique(x.fold):
        m = x.fold.values == fo
        v = x[col].values[m]
        d = np.sort(v[x.dec.values[m] == 1]); n = len(d)
        top = max(n // 100, 10); t = d[n - top]
        beta = max(np.mean(d[n - top:] - t), 1e-9)
        nls = lambda i: -np.log((n - i + 1) / (n + 1))
        xs, ys = [], []
        for i in range(0, n - top + 1):
            if xs and d[i] == xs[-1]: continue
            xs.append(d[i]); ys.append(nls(i))
        xs, ys = np.array(xs), np.array(ys); yt = nls(n - top)
        slope = (ys[1] - ys[0]) / max(xs[1] - xs[0], 1e-12) if len(xs) > 1 else 0.0
        j = np.searchsorted(xs, v, side='right')
        mid = np.where(j >= len(xs), ys[-1], ys[np.clip(j - 1, 0, len(xs) - 1)] + (ys[np.clip(j, 0, len(xs) - 1)] - ys[np.clip(j - 1, 0, len(xs) - 1)])
                       * (v - xs[np.clip(j - 1, 0, len(xs) - 1)]) / np.maximum(xs[np.clip(j, 0, len(xs) - 1)] - xs[np.clip(j - 1, 0, len(xs) - 1)], 1e-300))
        c = np.where(v > t, yt + (v - t) / beta, np.where(v < xs[0], ys[0] + (v - xs[0]) * slope, mid))
        out[m] = c
    return out
def tdc(x, col):
    x = x.copy()
    z = np.empty(len(x))
    for fo in np.unique(x.fold):
        m = x.fold.values == fo
        d = x[col].values[m & (x.dec.values == 1)]
        med = np.median(d); mad = 1.4826 * np.median(np.abs(d - med))
        z[m] = (x[col].values[m] - med) / mad
    x['z'] = z
    pw = x.sort_values(['z', 'dec'], ascending=[False, False]).drop_duplicates('key')
    d = np.cumsum(pw.dec.values); t = np.cumsum(1 - pw.dec.values)
    q = np.minimum(1, (d + 1) / np.maximum(t, 1))
    return pw.assign(q=np.minimum.accumulate(q[::-1])[::-1])
def metrics(pw):
    t = pw[pw.dec == 0].sort_values(['q', 'z'], ascending=[True, False])
    acc = t[t.q <= 0.01]
    fdp = np.cumsum(t.e.values) * 2 / np.arange(1, len(t) + 1)
    ok = np.where(fdp <= 0.01)[0]
    return len(acc), int(acc.e.sum()), int(ok.max() + 1) if len(ok) else 0
rows = []
for path in sys.argv[2:]:
    f = pd.read_csv(path, sep='\t', usecols=['precursor_id', 'is_decoy', 'fold', 'score', 'network_score'])
    n = f[f.network_score > -1e29].set_index('precursor_id')[['network_score', 'fold', 'is_decoy']]
    g = f[f.score > -1e29].set_index('precursor_id')[['score']]
    b = n.join(g, how='inner').join(P[['key', 'e']]).rename(columns={'is_decoy': 'dec'})
    for name, fn in COMBOS.items():
        b['c'] = fn(b.network_score.values, b.score.values)
        if name not in ('net', 'gbm'): b['c'] = calibrate_tail(b, 'c')
        N, E, N1 = metrics(tdc(b, 'c'))
        rows.append(dict(run=path.split('/')[0], combo=name, N=N, E=E, combined=2 * E / max(N, 1), N1=N1))
R = pd.DataFrame(rows)
R['raw'] = R.run.str[:5]
print(R.pivot(index='run', columns='combo', values='N').to_string())
print('per-raw median N@q<=0.01 / N@FDP_c1%:')
print(R.groupby(['raw', 'combo']).agg(N=('N', 'median'), N1=('N1', 'median'), c=('combined', 'median')).unstack('combo').to_string(float_format=lambda v: f'{v:.4f}' if v < 1 else f'{v:.0f}'))
print(R.pivot(index='run', columns='combo', values='N1').to_string())
print(R.pivot(index='run', columns='combo', values='combined').to_string(float_format=lambda v: f'{v:.4f}'))
print('pooled:')
for name, gq in R.groupby('combo', sort=False):
    print(f'  {name:6s} N {gq.N.sum()} E {gq.E.sum()} combined {2 * gq.E.sum() / gq.N.sum():.4f} worst {gq.combined.max():.4f} runs>1% {int((gq.combined > 0.01).sum())} '
          f'| N@FDP_c1% median {gq.N1.median():.0f} sum {gq.N1.sum()}')
