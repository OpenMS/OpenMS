"""Estimate the per-side boundary multiplier k from the run's own anchors (no DIA-NN, no labels).
raw01, rtmono run. Anchors: used calibration anchors (paired target-decoy q <= 1%), their 6-fragment consensus on
the search traces (traces.tsv; the C++ version will use the anchor traces stored by the calibration).

Per anchor, picked as the picker does (smoothed consensus, valley merge, 5% floor), then per side:
  h = half width at half maximum on that side (raw consensus, interpolated, scans)
  u = distance from the apex in units of h,  y = consensus / apex
Pooled over anchors and both sides (sides kept separate as well):
  signal profile  s(u) = mean y in bins of u (0.25 wide)
  background      b    = mean y over the trace outside the picked peak (all scans not in [left, right])
  k_cross = first u where s(u) <= b   (signal density = background density, as the fragment-tolerance fit)
Also k at which the profile's cumulative area reaches 95 / 99% of its total above background.
Then each estimate is scored against DIA-NN's boundaries (evaluation only), with bootstrap over anchors.
"""
import numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; X = B + '/extraction-parity-20261112'
out = open(X + '/rtmono/K_ESTIMATE_raw01.txt', 'w')
def say(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')

an = pd.read_csv(X + '/rtmono/cal/anchors.tsv', sep='\t'); an = an[an.used == 1]
us = pd.read_csv(X + '/rtmono/traces.tsv', sep='\t', usecols=['precursor_id', 'rt', 'intensity'])
have = set(us.precursor_id); ids = [p for p in an.precursor_id if p in have]
say(f'used anchors {len(an)}, with search traces {len(ids)}')
S = us[us.precursor_id.isin(ids)].groupby(['precursor_id', 'rt']).intensity.sum()

def smooth(v):
    y = v.copy(); y[1:-1] = 0.25 * v[:-2] + 0.5 * v[1:-1] + 0.25 * v[2:]; return y
def pick_one(v, near_idx, vf=0.5, bf=0.05):
    y = smooth(v); n = len(y); regions = []
    for s in range(n):
        if not (y[s] > 0) or not (s == 0 or y[s] > y[s - 1]) or not (s + 1 == n or y[s] >= y[s + 1]):
            continue
        cur = [s, 0]
        while regions:
            top = regions[-1]; val = top[0] + int(np.argmin(y[top[0]:cur[0] + 1]))
            if y[val] > vf * min(y[top[0]], y[cur[0]]):
                cur = [top[0] if y[top[0]] >= y[cur[0]] else cur[0], top[1]]; regions.pop(); continue
            cur[1] = val; break
        regions.append(cur)
    best = None
    for r, (top, vl) in enumerate(regions):
        vr = regions[r + 1][1] if r + 1 < len(regions) else n - 1
        l = top
        while l > vl and y[l] > bf * y[top]: l -= 1
        rr = top
        while rr < vr and y[rr] > bf * y[top]: rr += 1
        if l <= near_idx <= rr:
            best = (l, rr)
    if best is None:
        return None
    l, rr = best; a = l + int(np.argmax(v[l:rr + 1]))
    return a, l, rr
def half(v, a, l, r, step):
    h = 0.5 * v[a]; i = a
    lim = l if step < 0 else r
    while i != lim and v[i + step] > h: i += step
    if i == lim:
        return max(abs(i - a), 0.5)
    frac = (v[i] - h) / (v[i] - v[i + step]) if v[i] != v[i + step] else 0.0
    return max(abs(i - a) + frac, 0.5)

U = []; BG = []; per = []
for p in ids:
    s = S.loc[p]; t = s.index.values; v = np.clip(s.values, 0, None)
    rt = float(an.set_index('precursor_id').loc[p, 'run_rt'])
    res = pick_one(v, int(np.argmin(np.abs(t - rt))))
    if res is None or v[res[0]] <= 0:
        continue
    a, l, r = res
    hl, hr = half(v, a, l, r, -1), half(v, a, l, r, +1)
    y = v / v[a]
    out_mask = np.ones(len(v), bool); out_mask[l:r + 1] = False
    if out_mask.sum() >= 10:
        BG.append(y[out_mask].mean())
    idx = np.arange(len(v))
    # only inside the valley-bounded region on each side (stop at the neighbouring valley, not into the next peak)
    for side, hh, sl in (('L', hl, slice(None, a + 1)), ('R', hr, slice(a, None))):
        ii = idx[sl]
        u = np.abs(ii - a) / hh
        U.append(pd.DataFrame(dict(p=p, side=side, u=u, y=y[ii], inpeak=(ii >= l) & (ii <= r))))
    per.append(dict(p=p, hl=hl, hr=hr, l=a - l, r=r - a))
U = pd.concat(U); bg = float(np.mean(BG)); bg_med = float(np.median(BG))
say(f'anchors profiled {len(per)}; background mean y outside peak: mean over anchors {bg:.4f}, median {bg_med:.4f}')

def estimate(Uu, bg):
    Uu = Uu[Uu.u <= 12]
    bins = np.arange(0, 12.25, 0.25)
    prof = Uu.groupby(pd.cut(Uu.u, bins), observed=False).y.mean()
    centres = bins[:-1] + 0.125
    below = np.where(prof.values <= bg)[0]
    k_cross = centres[below[0]] if len(below) else np.nan
    sig = np.clip(prof.values - bg, 0, None); cs = np.cumsum(sig) / sig.sum()
    return k_cross, centres[np.searchsorted(cs, 0.95)], centres[np.searchsorted(cs, 0.99)], prof
kc, k95, k99, prof = estimate(U, bg)
kcL, _, _, profL = estimate(U[U.side == 'L'], bg); kcR, _, _, profR = estimate(U[U.side == 'R'], bg)
say(f'k_cross (profile = background) {kc:.3f}; left {kcL:.3f}, right {kcR:.3f}; area 95% at {k95:.3f}, 99% at {k99:.3f}')
say('pooled profile mean y by u:', {round(float(c.mid), 3): round(float(v), 4) for c, v in prof.items() if c.mid <= 6})

# bootstrap over anchors
rng = np.random.default_rng(1)
pl = np.array(sorted(U.p.unique())); g = dict(tuple(U.groupby('p')))
boots = []
BGs = pd.Series(BG)
for b in range(200):
    pick = rng.choice(pl, len(pl), replace=True)
    Ub = pd.concat([g[x] for x in pick]); bb = float(BGs.sample(len(BGs), replace=True, random_state=b).mean())
    boots.append(estimate(Ub, bb)[0])
say(f'bootstrap k_cross: median {np.nanmedian(boots):.3f}, 5-95% {np.nanquantile(boots, .05):.3f}-{np.nanquantile(boots, .95):.3f}')
# sensitivity to the number of anchors
for n in [50, 100, 200]:
    ks = [estimate(pd.concat([g[x] for x in rng.choice(pl, n, replace=False)]), bg)[0] for _ in range(100)]
    say(f'  with {n} anchors: k_cross median {np.nanmedian(ks):.3f}, 5-95% {np.nanquantile(ks, .05):.3f}-{np.nanquantile(ks, .95):.3f}')
pd.DataFrame(per).to_csv(X + '/rtmono/k_estimate_anchor_rows_raw01.tsv', sep='\t', index=False)

# score vs DIA-NN (evaluation only) using the cap_generalisation rows: side_k with k = estimate
R = pd.read_csv(X + '/rtmono/cap_generalisation_rows_raw01.tsv', sep='\t')
rows = pd.read_csv(X + '/rtmono/peak_parity_rows_raw01.tsv', sep='\t').set_index('precursor_id')
R = R.join(rows[['rt_apex', 'rt_start']], on='precursor_id')
apx = R.cur_l + np.round((R.rt_apex - R.rt_start) / 1.646)
for nm, k in [('k_cross', kc), ('k95', k95), ('k99', k99), ('fixed 2.5', 2.5)]:
    L = np.maximum(R.cur_l, apx - np.ceil(k * R.hl)); Rr = np.minimum(R.cur_r, apx + np.ceil(k * R.hr))
    dl, dr = L - R.dn_l, Rr - R.dn_r
    inter = (np.minimum(Rr, R.dn_r) - np.maximum(L, R.dn_l) + 1).clip(lower=0); un = np.maximum(Rr, R.dn_r) - np.minimum(L, R.dn_l) + 1
    say(f'  vs DIA-NN q<=1% (n={len(R)}), {nm} = {k:.2f}: start within1 {(dl.abs() <= 1).mean():.3f} | stop within1 {(dr.abs() <= 1).mean():.3f} '
        f'median {dr.median():+.0f} | IoU {(inter / un).median():.2f}')
