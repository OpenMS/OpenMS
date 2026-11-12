"""Does a boundary cap scaled to the peak's own width generalise across peak widths and shapes?
raw01, DIA-NN q <= 1%, DIA-NN apex inside our matched peak. Evaluation only.

Per peak, on our raw 6-fragment consensus (current boundaries from the picker):
  hl, hr   left / right half width at half maximum (interpolated), in scans
  rules    sym_k:  each side capped at k x FWHM (= k x (hl + hr))
           side_k: each side capped at k x its own half width (left at k x hl, right at k x hr), so a tailing peak
                   keeps a longer right side; for a Gaussian, k = 2.5 is 2.94 sigma per side
Reports, by FWHM tertile and by asymmetry (hr / hl) tertile: stop/start agreement with DIA-NN, the best k per
stratum, and the share of the current area each rule keeps.
Then a simulation: Gaussian and exponentially modified Gaussian (EMG) peaks at FWHM 1.5 to 20 scans, with
Poisson-like noise, to see whether the share of true area inside the cap depends on width.
"""
import numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; X = B + '/extraction-parity-20261112'
out = open(X + '/rtmono/CAP_GENERALISATION_raw01.txt', 'w')
def say(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')

rows = pd.read_csv(X + '/rtmono/peak_parity_rows_raw01.tsv', sep='\t')
q1 = rows[(rows.grp == 'q<=1%') & rows.contains]
us = pd.read_csv(X + '/rtmono/traces.tsv', sep='\t', usecols=['precursor_id', 'rt', 'intensity'])
us = us[us.precursor_id.isin(set(q1.precursor_id))]
S = us.groupby(['precursor_id', 'rt']).intensity.sum()

def halfwidths(v, a, l, r):
    h = 0.5 * v[a]
    i = a
    while i > l and v[i - 1] > h: i -= 1
    hl = (a - i) + ((v[i] - h) / (v[i] - v[i - 1]) if i > l and v[i] != v[i - 1] else 0.0)
    i = a
    while i < r and v[i + 1] > h: i += 1
    hr = (i - a) + ((v[i] - h) / (v[i] - v[i + 1]) if i < r and v[i] != v[i + 1] else 0.0)
    return max(hl, 0.5), max(hr, 0.5)

K_SYM = [0.75, 1.0, 1.25, 1.5, 2.0]; K_SIDE = [1.5, 2.0, 2.5, 3.0, 4.0]
res = []
for r in q1.itertuples():
    s = S.loc[r.precursor_id]; t = s.index.values; v = np.clip(s.values, 0, None)
    near = lambda x: int(np.argmin(np.abs(t - x)))
    a, l, rr, d0, d1 = near(r.rt_apex), near(r.rt_start), near(r.rt_end), near(r._5), near(r._6)
    hl, hr = halfwidths(v, a, l, rr)
    area = v[l:rr + 1].sum()
    row = dict(precursor_id=r.precursor_id, hl=hl, hr=hr, fwhm=hl + hr, asym=hr / hl, dn_l=d0, dn_r=d1, cur_l=l, cur_r=rr)
    for k in K_SYM:
        cl = max(l, a - int(np.ceil(k * (hl + hr)))); cr = min(rr, a + int(np.ceil(k * (hl + hr))))
        row[f'sym{k}_l'], row[f'sym{k}_r'], row[f'sym{k}_area'] = cl, cr, v[cl:cr + 1].sum() / area if area else np.nan
    for k in K_SIDE:
        cl = max(l, a - int(np.ceil(k * hl))); cr = min(rr, a + int(np.ceil(k * hr)))
        row[f'side{k}_l'], row[f'side{k}_r'], row[f'side{k}_area'] = cl, cr, v[cl:cr + 1].sum() / area if area else np.nan
    res.append(row)
R = pd.DataFrame(res)
R['w_bin'] = pd.qcut(R.fwhm, 3, labels=['narrow', 'mid', 'wide'])
R['a_bin'] = pd.qcut(R.asym, 3, labels=['fronting/sym', 'sym', 'tailing'])
say(f'n={len(R)}; our FWHM in scans by tertile:', R.groupby('w_bin', observed=True).fwhm.median().round(2).to_dict(),
    '| hr/hl by tertile:', R.groupby('a_bin', observed=True).asym.median().round(2).to_dict())
say('DIA-NN width / our FWHM by width tertile (does DIA-NN scale with width?):',
    R.groupby('w_bin', observed=True).apply(lambda g: ((g.dn_r - g.dn_l) / g.fwhm).median(), include_groups=False).round(2).to_dict())

def score(g, rule):
    if rule == 'current':
        L, Rr, ar = g.cur_l, g.cur_r, pd.Series(1.0, index=g.index)
    else:
        L, Rr, ar = g[rule + '_l'], g[rule + '_r'], g[rule + '_area']
    dl, dr = L - g.dn_l, Rr - g.dn_r
    inter = (np.minimum(Rr, g.dn_r) - np.maximum(L, g.dn_l) + 1).clip(lower=0); un = np.maximum(Rr, g.dn_r) - np.minimum(L, g.dn_l) + 1
    return dict(stop1=round(float((dr.abs() <= 1).mean()), 3), start1=round(float((dl.abs() <= 1).mean()), 3),
                stop_med=float(dr.median()), iou=round(float((inter / un).median()), 2), area=round(float(ar.median()), 3), area_p10=round(float(ar.quantile(.1)), 3))
rules = ['current'] + [f'sym{k}' for k in K_SYM] + [f'side{k}' for k in K_SIDE]
for by in ['w_bin', 'a_bin']:
    say(f'\n-- by {by}: stop within 1 scan / start within 1 / median stop offset / IoU / median share of current area kept (p10)')
    for nm, g in R.groupby(by, observed=True):
        say(f'  [{nm}]')
        best = {}
        for ru in rules:
            sc = score(g, ru)
            say(f'    {ru:9s} stop1 {sc["stop1"]:.3f} start1 {sc["start1"]:.3f} stop_med {sc["stop_med"]:+.0f} IoU {sc["iou"]:.2f} area {sc["area"]:.3f} ({sc["area_p10"]:.3f})')
            fam = ru.rstrip('0123456789.')
            if fam != 'current' and (fam not in best or sc['iou'] + sc['stop1'] > best[fam][1]):
                best[fam] = (ru, sc['iou'] + sc['stop1'])
        say('    best k by IoU + stop1:', {f: b[0] for f, b in best.items()})

# ---- simulation: true area inside the cap, across widths and tailing ----
say('\n-- simulation: share of the TRUE peak area inside the cap (noise-free shape, cap from a noisy FWHM estimate)')
rng = np.random.default_rng(7)
def emg(x, sigma, tau):
    if tau == 0:
        return np.exp(-0.5 * (x / sigma) ** 2)
    from scipy.special import erfc
    lam = 1 / tau
    return np.exp(lam / 2 * (lam * sigma ** 2 - 2 * x)) * erfc((lam * sigma ** 2 - x) / (np.sqrt(2) * sigma))
for tau_ratio in [0.0, 0.5, 1.0, 2.0]:
    line = []
    for fwhm in [1.5, 2.0, 3.0, 5.0, 10.0, 20.0]:
        kept = {ru: [] for ru in ['sym1.25', 'side2.5']}
        for rep in range(200):
            sigma = fwhm / 2.355; tau = tau_ratio * sigma
            off = rng.uniform(-0.5, 0.5)
            x = np.arange(-int(12 * fwhm) - 10, int(12 * fwhm) + 30) + off
            fine = np.linspace(x[0], x[-1], 20000); yf = emg(fine, sigma, tau); yf /= yf.max()
            y = emg(x, sigma, tau); y = y / y.max() * 1000
            y = rng.poisson(np.clip(y, 0, None) + 5).astype(float) - 5  # counting noise + flat background
            y = np.clip(y, 0, None)
            a = int(np.argmax(y)); hl, hr = halfwidths(y, a, 0, len(y) - 1)
            for ru in kept:
                if ru == 'sym1.25':
                    lo, hi = x[a] - np.ceil(1.25 * (hl + hr)), x[a] + np.ceil(1.25 * (hl + hr))
                else:
                    lo, hi = x[a] - np.ceil(2.5 * hl), x[a] + np.ceil(2.5 * hr)
                m = (fine >= lo) & (fine <= hi)
                kept[ru].append(yf[m].sum() / yf.sum())
        line.append(f'FWHM {fwhm:>4}: sym1.25 {np.median(kept["sym1.25"]):.3f} (p10 {np.quantile(kept["sym1.25"], .1):.3f}), '
                    f'side2.5 {np.median(kept["side2.5"]):.3f} (p10 {np.quantile(kept["side2.5"], .1):.3f})')
    say(f'  tau/sigma = {tau_ratio}:')
    for s_ in line:
        say('    ' + s_)
R.to_csv(X + '/rtmono/cap_generalisation_rows_raw01.tsv', sep='\t', index=False)
