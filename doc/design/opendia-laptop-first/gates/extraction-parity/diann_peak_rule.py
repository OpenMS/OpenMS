"""Which apex and boundary rule reproduces DIA-NN's RT / RT.Start / RT.Stop on DIA-NN's own exported traces?
raw01, q <= 1%, evaluation only. Rules are scored by exact scan agreement.

Apex candidates (within DIA-NN's [start, stop] and within +-6 scans of its apex):
  sum_raw     raw max of the sum of DIA-NN's exported fragments
  sum_s121    [1 2 1]/4 smoothed sum
  best_raw    raw max of the fragment with the largest area in [start, stop]
  best_s121   smoothed best fragment
  wsum_corr   sum of fragments weighted by max(0, r) with the best fragment over +-6 scans, raw max
  libtop3     sum of the 3 fragments with the highest library intensity, raw max
Boundary rules on a trace T (walk outward from DIA-NN's apex, at most R = 6 scans):
  stop at the first scan whose value <= f * T(apex) (that scan is the boundary) or where T rises (valley), f grid.
"""
import sqlite3, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; X = B + '/extraction-parity-20261112'
out = open(X + '/rtmono/DIANN_PEAK_RULE_raw01.txt', 'w')
def say(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')

rep = pd.read_parquet(X + '/diann/report_raw01.parquet', columns=['Precursor.Id', 'Q.Value', 'RT', 'RT.Start', 'RT.Stop'])
rep = rep[rep['Q.Value'] <= 0.01].set_index('Precursor.Id')
idm = pd.read_csv(X + '/idmap.tsv', sep='\t').set_index('Precursor.Id')
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
lib = pd.read_sql("""select tp.PRECURSOR_ID precursor_id, t.ANNOTATION annotation, t.LIBRARY_INTENSITY li
  from TRANSITION t join TRANSITION_PRECURSOR_MAPPING tp on tp.TRANSITION_ID=t.ID""", c)
lib = lib[lib.precursor_id.isin(idm.loc[rep.index].precursor_id)]
libi = {(p, a): v for p, a, v in zip(lib.precursor_id, lib.annotation, lib.li)}
dn = pd.read_parquet(X + '/diann/report_raw01_xic/raw01.xic.parquet')
dn = dn[dn.pr.isin(set(rep.index)) & ~dn.feature.isin(['index', 'ms1'])].copy()
dn['rt_s'] = dn.rt.astype(float) * 60
R = 6
def s121(v):
    return np.convolve(np.pad(v, 1, mode='edge'), [0.25, 0.5, 0.25], mode='valid')

stats = []
Fs = [0.0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3]
for pr, g in dn.groupby('pr'):
    P = g.pivot_table(index='rt_s', columns='feature', values='value', aggfunc='sum').fillna(0)
    t = P.index.values; F = P.values; feats = list(P.columns)
    r = rep.loc[pr]
    ia = int(np.argmin(np.abs(t - r.RT * 60))); i0 = int(np.argmin(np.abs(t - r['RT.Start'] * 60))); i1 = int(np.argmin(np.abs(t - r['RT.Stop'] * 60)))
    if abs(t[ia] - r.RT * 60) > 0.01 or abs(t[i0] - r['RT.Start'] * 60) > 0.01 or abs(t[i1] - r['RT.Stop'] * 60) > 0.01:
        continue
    lo, hi = max(0, ia - R), min(len(t) - 1, ia + R)
    S = F.sum(1)
    best = F[:, int(np.argmax(F[i0:i1 + 1].sum(0)))]
    seg = slice(lo, hi + 1)
    bc = best[seg] - best[seg].mean()
    w = []
    for k in range(F.shape[1]):
        x = F[seg, k] - F[seg, k].mean(); d = np.sqrt((x * x).sum() * (bc * bc).sum())
        w.append(max(0.0, (x * bc).sum() / d) if d > 0 else 0.0)
    W = F @ np.array(w)
    p = int(idm.loc[pr].precursor_id)
    li = np.array([libi.get((p, a), 0.0) for a in feats])
    top3 = F[:, np.argsort(-li)[:3]].sum(1)
    traces = dict(sum_raw=S, sum_s121=s121(S), best_raw=best, best_s121=s121(best), wsum_corr=W, libtop3=top3)
    row = dict(pr=pr, n_frag=F.shape[1], left=ia - i0, right=i1 - ia)
    for k, T in traces.items():
        row['apex_' + k] = int(np.argmax(T[seg])) + lo - ia
    # boundary rules evaluated from DIA-NN's apex
    for k in ['sum_raw', 'sum_s121', 'best_raw', 'wsum_corr']:
        T = traces[k]; top = T[ia]
        for f in Fs:
            def walk(step):
                j = ia
                for _ in range(R):
                    n = j + step
                    if n < 0 or n >= len(T):
                        break
                    if T[n] > T[j]:      # valley: rising again, boundary is the minimum j
                        break
                    j = n
                    if T[j] <= f * top:
                        break
                return abs(j - ia)
            row[f'L_{k}_{f}'] = walk(-1); row[f'R_{k}_{f}'] = walk(1)
    stats.append(row)
S = pd.DataFrame(stats)
say(f'q<=1% precursors with DIA-NN apex/start/stop on export scans: {len(S)}; DIA-NN width max {int((S.left + S.right).max())} scans, side max L {S.left.max()} R {S.right.max()}')
say('apex rule: share exactly at DIA-NN apex | within 1 scan')
for k in ['sum_raw', 'sum_s121', 'best_raw', 'best_s121', 'wsum_corr', 'libtop3']:
    d = S['apex_' + k]
    say(f'  {k:10s} {(d == 0).mean():.4f} | {(d.abs() <= 1).mean():.4f}')
say('boundary rule (from DIA-NN apex, cap 6 scans): share of sides exactly equal (left, right) and mean signed error (scans)')
best_rule = None
for k in ['sum_raw', 'sum_s121', 'best_raw', 'wsum_corr']:
    for f in Fs:
        L = S[f'L_{k}_{f}']; Rr = S[f'R_{k}_{f}']
        el, er = (L == S.left).mean(), (Rr == S.right).mean()
        say(f'  {k:10s} f={f:<5} L {el:.3f} R {er:.3f} | bias L {(L - S.left).mean():+.2f} R {(Rr - S.right).mean():+.2f} | within 1: L {((L - S.left).abs() <= 1).mean():.3f} R {((Rr - S.right).abs() <= 1).mean():.3f}')
        if best_rule is None or el + er > best_rule[0]:
            best_rule = (el + er, k, f)
say('best rule by exact agreement:', best_rule)
S.to_csv(X + '/rtmono/diann_peak_rule_rows_raw01.tsv', sep='\t', index=False)
