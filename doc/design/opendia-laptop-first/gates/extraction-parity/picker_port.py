"""Python port of pickXicPeaks (OpenDIAPeakPicker.cpp) on our dumped traces, validated against our reported peaks,
then picker variants scored against DIA-NN's peaks (raw01, DIA-NN q <= 1%, apex inside our RT window). Evaluation only.

Variants (each changes one thing):
  current       the C++ picker
  profile2      a fragment enters co-elution only if it has signal on >= 2 scans of the correlation window
                (Pearson r over a single nonzero point is 1 for any coincidence: an impulse is not a profile)
  cut_valley    boundaries stop at the first smoothed minimum (no shoulder merge for boundaries; ranking unchanged)
Usage: picker_port.py
"""
import numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; X = B + '/extraction-parity-20261112'
out = open(X + '/rtmono/PICKER_VARIANTS_raw01.txt', 'w')
def say(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')

def pick(M, rts, valley_fraction=0.5, boundary_fraction=0.05, h=2, max_peaks=5, profile_min=1, cut_valley=False):
    scans, F = M.shape
    raw = M.sum(1).astype(float)
    y = raw.copy()
    t = y.copy()
    for s in range(1, scans - 1):
        t[s] = 0.25 * y[s - 1] + 0.5 * y[s] + 0.25 * y[s + 1]
    y = t
    regions = []
    for s in range(scans):
        rises = s == 0 or y[s] > y[s - 1]
        falls = s + 1 == scans or y[s] >= y[s + 1]
        if not (y[s] > 0) or not rises or not falls:
            continue
        cur = [s, 0]
        while regions:
            top = regions[-1]
            seg = y[top[0]:cur[0] + 1]; valley = top[0] + int(np.argmin(seg))
            if y[valley] > valley_fraction * min(y[top[0]], y[cur[0]]):
                cur[0] = top[0] if y[top[0]] >= y[cur[0]] else cur[0]
                cur[1] = top[1]; regions.pop(); continue
            cur[1] = valley; break
        regions.append(cur)
    peaks = []
    for r, (top, vl) in enumerate(regions):
        vr = regions[r + 1][1] if r + 1 < len(regions) else scans - 1
        fl = boundary_fraction * y[top]
        left = top
        while left > vl and y[left] > fl and not (cut_valley and y[left - 1] > y[left]):
            left -= 1
        right = top
        while right < vr and y[right] > fl and not (cut_valley and y[right + 1] > y[right]):
            right += 1
        apex = left + int(np.argmax(raw[left:right + 1]))
        if not raw[apex] > 0:
            continue
        lo = min(left, apex - h if apex > h else 0); hi = max(right, min(scans - 1, apex + h))
        co = 0.0
        if hi - lo >= 2:
            W = M[lo:hi + 1].astype(float); tot = W.sum(1)
            for k in range(F):
                x = W[:, k]
                if (x > 0).sum() < profile_min:
                    continue
                o = tot - x
                vx = x.var(); vo = o.var()
                if vx > 0 and vo > 0:
                    co += max(0.0, float(np.mean((x - x.mean()) * (o - o.mean())) / np.sqrt(vx * vo)))
        peaks.append(dict(apex=apex, left=left, right=right, h=raw[apex], co=co))
    peaks.sort(key=lambda p: (-p['co'], -p['h'], p['apex']))
    return peaks[:max_peaks], len(regions)

rows = pd.read_csv(X + '/rtmono/peak_parity_rows_raw01.tsv', sep='\t')
q1 = rows[rows.grp == 'q<=1%']
us = pd.read_csv(X + '/rtmono/traces.tsv', sep='\t', usecols=['precursor_id', 'scan', 'rt', 'slot', 'intensity'])
us = us[us.precursor_id.isin(set(q1.precursor_id))]
feat = pd.read_csv(X + '/rtmono/search.features.tsv', sep='\t', usecols=['precursor_id', 'is_decoy', 'rt_apex', 'rt_start', 'rt_end', 'coelution', 'peak_rank'])
feat = feat[(feat.is_decoy == 0) & feat.precursor_id.isin(set(q1.precursor_id))]
fg = dict(tuple(feat.groupby('precursor_id'))); ug = dict(tuple(us.groupby('precursor_id')))

variants = {'current': {}, 'profile2': dict(profile_min=2), 'cut_valley': dict(cut_valley=True), 'profile2+cut_valley': dict(profile_min=2, cut_valley=True)}
agree = 0; n = 0
res = {v: [] for v in variants}
for r in q1.itertuples():
    g = ug[r.precursor_id]
    P = g.pivot_table(index='rt', columns='slot', values='intensity', aggfunc='sum').fillna(0).clip(lower=0)
    rts = P.index.values; M = P.values
    inwin = rts[0] <= r.RT <= rts[-1]
    for v, kw in variants.items():
        pk, _ = pick(M, rts, **kw)
        if v == 'current':
            rep = fg[r.precursor_id].sort_values('peak_rank')
            mine = [(rts[p['apex']], p['co']) for p in pk]
            theirs = list(zip(rep.rt_apex, rep.coelution))
            n += 1; agree += len(mine) == len(theirs) and all(abs(a - b) < 1e-3 and abs(c - d) < 1e-4 for (a, c), (b, d) in zip(mine, theirs))
        ia = int(np.argmin(np.abs(rts - r.RT))); d0 = int(np.argmin(np.abs(rts - r._5))); d1 = int(np.argmin(np.abs(rts - r._6)))
        cont = [i for i, p in enumerate(pk) if p['left'] <= ia <= p['right']]
        top = pk[0] if pk else None
        res[v].append(dict(precursor_id=r.precursor_id, e=r.e, inwin=inwin, contained=bool(cont), top_ok=bool(cont) and cont[0] == 0,
                           top_apex_d=(top['apex'] - ia) if top else np.nan,
                           dl=(pk[cont[0]]['left'] - d0) if cont else np.nan, dr=(pk[cont[0]]['right'] - d1) if cont else np.nan))
say(f'port check: our 5 reported peaks (apex RT and co-elution, in rank order) reproduced for {agree} of {n} precursors')
for v in variants:
    R = pd.DataFrame(res[v]); W = R[R.inwin]
    say(f'== {v}: DIA-NN apex in window {len(W)}: contained in one of 5 peaks {W.contained.sum()} | our rank-0 is DIA-NN\'s peak {W.top_ok.sum()} '
        f'(entrapments {int(W[W.top_ok].e.sum())} of {int(W.e.sum())}) | rank-0 apex exact {(W.top_apex_d == 0).sum()} within 1 {(W.top_apex_d.abs() <= 1).sum()}')
    C = W[W.contained]
    say(f'   boundaries of the containing peak vs DIA-NN, scans: start exact {(C.dl == 0).mean():.3f} median {C.dl.median():+.0f} | '
        f'stop exact {(C.dr == 0).mean():.3f} median {C.dr.median():+.0f} p90 {C.dr.quantile(.9):+.0f}')
