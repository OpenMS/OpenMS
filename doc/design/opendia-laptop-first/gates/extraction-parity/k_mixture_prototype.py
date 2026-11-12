import numpy as np, pandas as pd
exec(open('k_estimate.py').read().split("S = us[us.precursor_id.isin(ids)]")[0])
us = pd.read_csv(X + '/rtmono/traces.tsv', sep='\t', usecols=['precursor_id', 'rt', 'slot', 'intensity'])
us = us[us.precursor_id.isin(ids)]
exec(open('k_estimate.py').read().split("def smooth(v):")[1].split("U = []; BG = []; per = []")[0].join(["def smooth(v):",""]))
A = an.set_index('precursor_id')
ks = np.arange(0.5, 8.01, 0.25)
snr = []; cosrows = []; outcos = []
for p, g in us.groupby('precursor_id'):
    M = g.pivot_table(index='rt', columns='slot', values='intensity', aggfunc='sum').fillna(0).clip(lower=0)
    t = M.index.values; F = M.values; v = F.sum(1)
    res = pick_one(v, int(np.argmin(np.abs(t - A.loc[p, 'run_rt']))))
    if res is None or v[res[0]] <= 0: continue
    a, l, r = res; hl, hr = half(v, a, l, r, -1), half(v, a, l, r, 1)
    om = np.ones(len(v), bool); om[max(0, l - 2):r + 3] = False
    if om.sum() < 10: continue
    sd = v[om].std() / v[a]
    if not sd > 0: sd = 1.0 / v[a]  # zero background: one count-level unit
    y = v / v[a]
    row = []
    for k in ks:
        lo = max(0, a - int(np.ceil(k * hl))); hi = min(len(v) - 1, a + int(np.ceil(k * hr)))
        row.append(y[lo:hi + 1].sum() / (sd * np.sqrt(hi - lo + 1)))
    snr.append(np.array(row) / max(row))
    ref = F[a] / (np.linalg.norm(F[a]) + 1e-12)
    for j in range(len(v)):
        n = np.linalg.norm(F[j])
        if n <= 0: continue
        c = float(F[j] @ ref / n)
        u = (a - j) / hl if j < a else (j - a) / hr
        if om[j]: outcos.append((c, p))
        else: cosrows.append((u, c, p))
snr = np.array(snr)
med = np.median(snr, 0)
print('anchors', len(snr))
print('SNR-optimal k: argmax of median normalised SNR =', ks[int(np.argmax(med))], '| per-anchor argmax quantiles', pd.Series(ks[np.argmax(snr, 1)]).quantile([.1,.25,.5,.75,.9]).to_dict())
print('  median normalised SNR by k:', {float(k): round(float(m), 3) for k, m in zip(ks, med) if k in (1, 1.5, 2, 2.5, 3, 4, 5, 6, 8)})
C = pd.DataFrame(cosrows, columns=["u", "c", "pid"]); OC = pd.DataFrame(outcos, columns=["c", "pid"]); oc = OC.c.values
print('cosine to apex fragment pattern: out-of-peak scans median', round(float(np.median(oc)), 3), 'p75', round(float(np.quantile(oc, .75)), 3), 'n', len(oc))
prof = C.groupby(pd.cut(C.u, np.arange(0, 8.25, 0.5)), observed=False).c.median()
print('  in-peak-region median cosine by u:', {round(float(i.mid), 2): round(float(v), 3) for i, v in prof.items()})

# Mixture: cosine histogram at distance u = pi(u) * signal (u < 1) + (1 - pi(u)) * background (out of peak).
bins = np.linspace(-0.0001, 1.0001, 21)
hs = np.histogram(C[C.u < 1].c, bins)[0].astype(float); hs /= hs.sum()
hb = np.histogram(oc, bins)[0].astype(float); hb /= hb.sum()
def pi_of(cs):
    h = np.histogram(cs, bins)[0].astype(float)
    if h.sum() < 30: return np.nan
    h /= h.sum(); d = hs - hb
    return float(np.clip(((h - hb) @ d) / (d @ d), 0, 1))
grid = np.arange(0.5, 8.01, 0.5)
pis = {float(u): pi_of(C[(C.u >= u - 0.25) & (C.u < u + 0.25)].c) for u in grid}
print('signal share pi(u):', {k: round(v, 3) for k, v in pis.items()})
gs = sorted(pis); k_half = next((u for u in gs if pis[u] <= 0.5), np.nan)
print('k where pi(u) falls to 0.5:', k_half)
# bootstrap over anchors

def k_est(Cc, ocb, step=0.25):
    hs_ = np.histogram(Cc[Cc.u < 1].c, bins)[0].astype(float); hs_ /= hs_.sum()
    hb_ = np.histogram(ocb, bins)[0].astype(float); hb_ /= hb_.sum(); d = hs_ - hb_
    prev_u, prev_p = None, None
    for u in np.arange(1.0, 8.01, step):
        cs = Cc[(Cc.u >= u - step) & (Cc.u < u + step)].c
        h = np.histogram(cs, bins)[0].astype(float)
        if h.sum() < 30: continue
        p = float(np.clip(((h / h.sum() - hb_) @ d) / (d @ d), 0, 1))
        if p <= 0.5:
            return u if prev_p is None else prev_u + (prev_p - 0.5) / (prev_p - p) * (u - prev_u)
        prev_u, prev_p = u, p
    return np.nan
C['p'] = 0
print('k (interpolated, step 0.25):', round(k_est(C, oc), 3))
rng = np.random.default_rng(3); pids = C.pid.unique(); Cg = dict(tuple(C.groupby('pid'))); Og = dict(tuple(OC.groupby('pid')))
for n in [len(pids), 200, 100, 50]:
    ks_ = []
    for b in range(100):
        pick = rng.choice(pids, n, replace=n == len(pids))
        ks_.append(k_est(pd.concat([Cg[x] for x in pick]), pd.concat([Og[x] for x in pick if x in Og]).c.values))
    print(f'  bootstrap n={n}: k median {np.nanmedian(ks_):.3f} 5-95% {np.nanquantile(ks_, .05):.3f}-{np.nanquantile(ks_, .95):.3f}')
