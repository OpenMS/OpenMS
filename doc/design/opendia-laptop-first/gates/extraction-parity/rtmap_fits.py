"""RT map fits on the rtedge anchors (raw01). Label-free criterion: 10-fold cross-validated absolute residual of held-out
anchors, overall and for the 30 lowest / highest library-RT anchors. DIA-NN apex residuals are evaluation only."""
import numpy as np, pandas as pd, sqlite3
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; X = B + '/extraction-parity-20261112'
a = pd.read_csv(X + '/rtedge/cal/anchors.tsv', sep='\t').sort_values(['library_rt', 'run_rt'])
xa, ya = a.library_rt.values, a.run_rt.values

def interp_lin(kx, ky, x, span=3):
    x = np.asarray(x, float); y = np.interp(x, kx, ky); n = len(kx); s = min(span, n - 1)
    sl = max((ky[s] - ky[0]) / max(kx[s] - kx[0], 1e-12), 0); sh = max((ky[-1] - ky[-1 - s]) / max(kx[-1] - kx[-1 - s], 1e-12), 0)
    return np.where(x < kx[0], ky[0] + sl * (x - kx[0]), np.where(x > kx[-1], ky[-1] + sh * (x - kx[-1]), y))

def binned(x, y, bins):
    n = len(x); kx, ky = [], []
    for b in range(bins):
        lo, hi = b * n // bins, (b + 1) * n // bins
        if hi > lo: kx.append(np.median(x[lo:hi])); ky.append(np.median(y[lo:hi]))
    return np.array(kx), np.array(ky)

def current(x, y):
    bins = int(np.clip(len(x) // 10, 3, 20)); keep = np.ones(len(x), bool)
    for _ in range(4):
        kx, ky = binned(x[keep], y[keep], bins); r = np.abs(y - interp_lin(kx, ky, x))
        sig = 1.4826 * np.median(r[keep]); k2 = r <= max(3 * sig, 3)
        if (k2 == keep).all(): break
        keep = k2
    return kx, ky

def pava_l1(y):
    # isotonic (nondecreasing) regression with L1 loss: blocks hold their median
    blocks = []
    for v in y:
        blocks.append([v])
        while len(blocks) > 1 and np.median(blocks[-2]) > np.median(blocks[-1]):
            last = blocks.pop(); blocks[-1].extend(last)
    return np.concatenate([[np.median(b)] * len(b) for b in blocks])

def isotonic(x, y, smooth=0):
    f = pava_l1(y)
    if smooth:
        f = pd.Series(f).rolling(smooth, center=True, min_periods=1).median().values
        f = np.maximum.accumulate(f)
    # knots at unique x (median of fitted values)
    d = pd.DataFrame({'x': x, 'f': f}).groupby('x').f.median()
    return d.index.values, d.values

def binned_ends(x, y):
    # current fit, but the end bins split by halving while they hold > ~sqrt(n) anchors, so tails get their own knots
    bins = int(np.clip(len(x) // 10, 3, 20)); n = len(x); edges = [0]
    for b in range(1, bins): edges.append(b * n // bins)
    edges.append(n)
    m = int(np.sqrt(n))
    while edges[1] - edges[0] > m: edges.insert(1, (edges[0] + edges[1]) // 2)
    while edges[-1] - edges[-2] > m: edges.insert(-1, (edges[-2] + edges[-1]) // 2)
    keep = np.ones(n, bool)
    for _ in range(4):
        kx, ky = [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            s = np.arange(lo, hi); s = s[keep[s]]
            if len(s): kx.append(np.median(x[s])); ky.append(np.median(y[s]))
        kx, ky = np.array(kx), np.maximum.accumulate(np.array(ky)); r = np.abs(y - interp_lin(kx, ky, x))
        sig = 1.4826 * np.median(r[keep]); k2 = r <= max(3 * sig, 3)
        if (k2 == keep).all(): break
        keep = k2
    return kx, ky

fits = {'current': current, 'isotonic_L1': lambda x, y: isotonic(x, y), 'isotonic_L1_smooth27': lambda x, y: isotonic(x, y, 27),
        'binned_split_ends': binned_ends}
rng = np.random.default_rng(0); fold = rng.permutation(len(xa)) % 10
low = np.argsort(xa)[:30]; high = np.argsort(xa)[-30:]
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
lrt = pd.read_sql('select ID precursor_id, LIBRARY_RT lrt from PRECURSOR', c).set_index('precursor_id').lrt
dn = pd.read_csv(B + '/diann-ourlib-20261112/diann_ourlib_raw01.tsv.gz', sep='\t', usecols=['precursor_id', 'RT', 'Q.Value']).drop_duplicates('precursor_id')
dn = dn[dn['Q.Value'] <= 0.01]; dx = lrt.reindex(dn.precursor_id).values; dy = dn.RT.values * 60
bins = [-1, 0.01, 0.02, 0.03, 0.05, 0.1, 0.5, 0.8, 0.87, 1]
out = open(X + '/rtedge/RTMAP_FITS_raw01.txt', 'w')
def say(s): print(s); out.write(s + '\n')
say(f'anchors {len(xa)}; library RT {xa.min():.4f}..{xa.max():.4f}')
for name, f in fits.items():
    pred = np.empty(len(xa))
    for k in range(10):
        tr = fold != k; kx, ky = f(xa[tr], ya[tr]); pred[~tr] = interp_lin(kx, ky, xa[~tr])
    r = np.abs(ya - pred)
    kx, ky = f(xa, ya); dr = np.abs(dy - interp_lin(kx, ky, dx))
    by = pd.Series(dr).groupby(pd.cut(dx, bins)).median().round(1).to_dict()
    say(f'{name}: CV |res| median {np.median(r):.1f} p95 {np.quantile(r, .95):.1f} | lowest-30 median {np.median(r[low]):.1f} | highest-30 {np.median(r[high]):.1f} s'
        f' || DIA-NN q<=1% |apex - map| median {np.median(dr):.1f}, p95 {np.quantile(dr, .95):.1f}; by library RT {by}')
out.close()

# Isotonic L1 knots, extrapolated with the slope of a Theil-Sen line through the sqrt(n) end anchors on each side.
def theil_sen(x, y):
    i, j = np.triu_indices(len(x), 1); d = x[j] - x[i]; ok = d > 1e-9
    return np.median((y[j] - y[i])[ok] / d[ok]) if ok.any() else 0.0
def iso_ts(x, y):
    o = np.argsort(x, kind='stable'); x, y = x[o], y[o]
    kx, ky = isotonic(x, y); m = max(int(np.sqrt(len(x))), 3)
    return kx, ky, max(theil_sen(x[:m], y[:m]), 0), max(theil_sen(x[-m:], y[-m:]), 0)
def ev(kx, ky, sl, sh, x):
    x = np.asarray(x, float); y = np.interp(x, kx, ky)
    return np.where(x < kx[0], ky[0] + sl * (x - kx[0]), np.where(x > kx[-1], ky[-1] + sh * (x - kx[-1]), y))
out = open(X + '/rtedge/RTMAP_FITS_raw01.txt', 'a')
pred = np.empty(len(xa))
for k in range(10):
    tr = fold != k; p = iso_ts(xa[tr], ya[tr]); pred[~tr] = ev(*p, xa[~tr])
r = np.abs(ya - pred); p = iso_ts(xa, ya); dr = np.abs(dy - ev(*p, dx))
by = pd.Series(dr).groupby(pd.cut(dx, bins)).median().round(1).to_dict()
say(f'isotonic_L1+TheilSen_ends: CV |res| median {np.median(r):.1f} p95 {np.quantile(r, .95):.1f} | lowest-30 median {np.median(r[low]):.1f} | highest-30 {np.median(r[high]):.1f} s'
    f' || DIA-NN q<=1% |apex - map| median {np.median(dr):.1f}, p95 {np.quantile(dr, .95):.1f}; by library RT {by}; end slopes {p[2]:.0f} {p[3]:.0f}')
# robustness: one false anchor placed at the lowest library RT, 300 s late
xb = np.append(xa, xa.min() - 0.001); yb = np.append(ya, 500.0); p = iso_ts(xb, yb)
say(f'  with one false early anchor (lib {xa.min()-0.001:.4f} at 500 s): map at lib 0.005 = {ev(*p, [0.005])[0]:.0f} s (clean {ev(*iso_ts(xa, ya), [0.005])[0]:.0f} s)')
out.close()

# Block-centre knots: each PAVA block becomes one knot at (median library RT, block median run RT); fewer, stable knots.
def iso_blocks(x, y):
    o = np.argsort(x, kind='stable'); x, y = x[o], y[o]
    f = pava_l1(y); kx, ky = [], []
    start = 0
    for i in range(1, len(f) + 1):
        if i == len(f) or f[i] != f[start]:
            kx.append(np.median(x[start:i])); ky.append(f[start]); start = i
    kx, ky = np.array(kx), np.array(ky); m = max(int(np.sqrt(len(x))), 3)
    return kx, ky, max(theil_sen(x[:m], y[:m]), 0), max(theil_sen(x[-m:], y[-m:]), 0)
out = open(X + '/rtedge/RTMAP_FITS_raw01.txt', 'a')
pred = np.empty(len(xa))
for k in range(10):
    tr = fold != k; p = iso_blocks(xa[tr], ya[tr]); pred[~tr] = ev(*p, xa[~tr])
r = np.abs(ya - pred); p = iso_blocks(xa, ya); dr = np.abs(dy - ev(*p, dx)); res = ya - ev(*p, xa)
by = pd.Series(dr).groupby(pd.cut(dx, bins)).median().round(1).to_dict()
say(f'isotonic_L1_blocks+TheilSen_ends: knots {len(p[0])} | CV |res| median {np.median(r):.1f} p95 {np.quantile(r, .95):.1f} | lowest-30 {np.median(r[low]):.1f} | highest-30 {np.median(r[high]):.1f} s'
    f' || DIA-NN q<=1% |apex - map| median {np.median(dr):.1f}, p95 {np.quantile(dr, .95):.1f}; by library RT {by}; in-sample sigma {1.4826*np.median(np.abs(res)):.2f} s')
out.close()

# Simulated search window: monotone map (blocks + Theil-Sen ends), sigma from leave-one-out residuals, half width 4 sigma.
out = open(X + '/rtedge/RTMAP_FITS_raw01.txt', 'a')
loo = np.array([ya[i] - ev(*iso_blocks(np.delete(xa, i), np.delete(ya, i)), [xa[i]])[0] for i in range(len(xa))])
sig = 1.4826 * np.median(np.abs(loo)); keep = np.abs(loo) <= max(3 * sig, 3)
p = iso_blocks(xa[keep], ya[keep])
loo2 = np.array([ya[keep][i] - ev(*iso_blocks(np.delete(xa[keep], i), np.delete(ya[keep], i)), [xa[keep][i]])[0] for i in range(keep.sum())])
sig2 = 1.4826 * np.median(np.abs(loo2)); tol = 4 * sig2
ent = pd.read_csv(B + '/diann-ourlib-20261112/diann_ourlib_raw01.tsv.gz', sep='\t', usecols=['precursor_id', 'Q.Value', 'e']).drop_duplicates('precursor_id')
ee = ent.set_index('precursor_id').e.reindex(dn.precursor_id).values
for name, centre, t in [('current', interp_lin(*current(xa, ya), dx), 84.66), ('monotone', ev(*p, dx), tol)]:
    outside = np.abs(dy - centre) > t
    say(f'search window sim [{name}]: half width {t:.1f} s; DIA-NN q<=1% apex outside {outside.sum()} (entrapments {int(ee[outside].sum())})')
say(f'  monotone: LOO sigma {sig:.2f} s -> kept {keep.sum()} of {len(xa)}; refit LOO sigma {sig2:.2f} s; knots {len(p[0])}')
out.close()
out = open(X + '/rtedge/RTMAP_FITS_raw01.txt', 'a')
c2 = ev(*p, dx); outside = np.abs(dy - c2) > tol
say('  monotone: outside by library RT ' + str(pd.Series(outside).groupby(pd.cut(dx, bins)).sum().to_dict()))
say('  monotone: their |apex - map| s ' + str(np.round(np.sort(np.abs(dy - c2)[outside]), 0).tolist()))
out.close()
