#!/usr/bin/env python3
# Port of fitMonotoneRTMap (OpenDIARTMap.cpp) with three outlier-rejection variants, evaluated on the rtmono
# stage-3 anchors against DIA-NN q<=1% apexes (simulated window: map +- 4 sigma, as the search uses).
#   insample : |y - full isotonic fit| (what rtmono ran)
#   cv       : 10-fold interleaved held-out residual (fold maps extrapolate with Theil-Sen over sqrt(n) ends)
#   cvlocal  : held-out residual, but a held-out point beyond the fold map's end knots is predicted from the
#              median pairwise slope of its 3 nearest training anchors on the inner side (a local slope)
import sys, sqlite3, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
X = B + '/extraction-parity-20261112'

def med(v): return float(np.median(v))
def theil(p):
    s = [(p[j][1] - p[i][1]) / (p[j][0] - p[i][0]) for i in range(len(p)) for j in range(i + 1, len(p)) if p[j][0] > p[i][0]]
    return med(s) if s else 0.0
def knots(p):
    blocks = []
    for x, y in p:
        blocks.append([[x], [y], y])
        while len(blocks) > 1 and blocks[-2][2] > blocks[-1][2]:
            l = blocks.pop(); b = blocks[-1]; b[0] += l[0]; b[1] += l[1]; b[2] = med(b[1])
    kx, ky = [], []
    for bx, by, lv in blocks:
        x = med(bx)
        if kx and x <= kx[-1]:
            ky[-1] = max(ky[-1], lv); continue
        kx.append(x); ky.append(max(lv, ky[-1]) if ky else lv)
    m = min(max(3, int(round(np.sqrt(len(p))))), len(p))
    return np.array(kx), np.array(ky), max(theil(p[:m]), 0.0), max(theil(p[-m:]), 0.0)
def ev(M, x):
    kx, ky, sl, sh = M; x = np.asarray(x, float)
    y = np.interp(x, kx, ky)
    y = np.where(x < kx[0], ky[0] + sl * (x - kx[0]), y)
    return np.where(x > kx[-1], ky[-1] + sh * (x - kx[-1]), y)
def local_end(train, x):
    tx = np.array([t[0] for t in train]); ty = np.array([t[1] for t in train])
    if x < tx[0]:
        q = list(zip(tx[:3], ty[:3]))
    else:
        q = list(zip(tx[-3:], ty[-3:]))
    s = max(theil(q), 0.0)
    ax, ay = (q[0] if x < tx[0] else q[-1])
    return ay + s * (x - ax)
def fit(anchors, variant):
    a = sorted(anchors)
    for rnd in range(4):
        n = len(a); res = np.zeros(n)
        for f in range(10):
            tr = [a[i] for i in range(n) if i % 10 != f]; M = knots(tr)
            for i in range(f, n, 10):
                x, y = a[i]
                if variant == 'cvlocal' and (x < M[0][0] or x > M[0][-1]):
                    res[i] = abs(y - local_end(tr, x))
                else:
                    res[i] = abs(y - ev(M, x))
        M = knots(a); sigma = 1.4826 * med(res); cut = max(3 * sigma, 3)
        r = np.abs(np.array([y for _, y in a]) - ev(M, [x for x, _ in a])) if variant == 'insample' else res
        kept = [a[i] for i in range(n) if r[i] <= cut]
        dropped = [a[i] for i in range(n) if r[i] > cut]
        if len(kept) == n: break
        a = kept
    return M, sigma, a

c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
lrt = pd.read_sql('select ID precursor_id, LIBRARY_RT lrt from PRECURSOR', c).set_index('precursor_id').lrt
d = pd.read_csv(B + '/diann-ourlib-20261112/diann_ourlib_raw01.tsv.gz', sep='\t', usecols=['precursor_id', 'RT', 'Q.Value', 'e'])
d = d[d['Q.Value'] <= 0.01].copy(); d['rt_s'] = d.RT * 60; d['lrt'] = d.precursor_id.map(lrt)
an = pd.read_csv(X + '/rtmono/cal/anchors.tsv', sep='\t')
pts = list(zip(an.library_rt, an.run_rt))
bins = [-1, 0.02, 0.04, 0.06, 0.1, 0.3, 0.6, 0.8, 0.87, 1.0]
out = open(X + '/rtmono/REJECT_VARIANTS_raw01.txt', 'w')
def say(s): print(s); out.write(s + '\n')
say(f'stage-3 anchors {len(pts)}; DIA-NN q<=1% {len(d)}')
for v in ['insample', 'cv', 'cvlocal']:
    M, sigma, kept = fit(pts, v)
    keptset = set(kept); drop = [p for p in sorted(pts) if p not in keptset]
    pred = ev(M, d.lrt.values); err = np.abs(d.rt_s.values - pred); tol = 4 * sigma
    out_w = err > tol
    by = pd.Series(err).groupby(pd.cut(d.lrt.values, bins), observed=False).median().round(1).to_dict()
    miss = pd.Series(out_w).groupby(pd.cut(d.lrt.values, bins), observed=False).sum().to_dict()
    say(f'== {v}: used {len(kept)} sigma {sigma:.2f} s tol {tol:.1f} s knots {len(M[0])} slopes {M[2]:.0f}/{M[3]:.0f}; '
        f'DIA-NN apex outside window {int(out_w.sum())} (entrapments {int((out_w & (d.e.values == 1)).sum())})')
    say('   dropped (lib, run): ' + ', '.join(f'({x:.4f},{y:.0f})' for x, y in drop))
    say('   misses by lib bin: ' + str({str(k): int(val) for k, val in miss.items()}))
    say('   median |apex - map| by lib bin: ' + str({str(k): val for k, val in by.items()}))
