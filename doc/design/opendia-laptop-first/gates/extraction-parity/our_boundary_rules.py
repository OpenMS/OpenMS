"""Why our right boundary runs long, on OUR traces (raw01, q <= 1%, DIA-NN apex inside our matched peak).
Walk outward from OUR apex on our [1 2 1]-smoothed 6-fragment consensus with several rules and compare with
DIA-NN's RT.Start / RT.Stop in scans:
  current      our reported rt_start / rt_end (valley merge at 0.5 + 5% floor, no span cap)
  valley       stop at the first local minimum (no merging of shoulders) or 5% floor
  valley_cap   valley, capped at DIA-NN's scan window radius (6, from its log; DIA-NN derives it from the peak width)
  f10_cap      valley or 10% floor, cap 6
Also lists the q<=1% precursors whose DIA-NN apex is not inside any of our peaks, by cause.
"""
import numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; X = B + '/extraction-parity-20261112'
out = open(X + '/rtmono/OUR_BOUNDARY_RULES_raw01.txt', 'w')
def say(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')
rows = pd.read_csv(X + '/rtmono/peak_parity_rows_raw01.tsv', sep='\t')
q1 = rows[rows.grp == 'q<=1%']
us = pd.read_csv(X + '/rtmono/traces.tsv', sep='\t', usecols=['precursor_id', 'rt', 'intensity'])
us = us[us.precursor_id.isin(set(q1.precursor_id))]
S = us.groupby(['precursor_id', 'rt']).intensity.sum()
def s121(v): return np.convolve(np.pad(v, 1, mode='edge'), [0.25, 0.5, 0.25], mode='valid')
def near(t, x): i = int(np.argmin(np.abs(t - x))); return i
res = []
for r in q1[q1.contains].itertuples():
    s = S.loc[r.precursor_id]; t = s.index.values; v = s121(s.values)
    ia = near(t, r.rt_apex); d0 = near(t, r._5); d1 = near(t, r._6); o0 = near(t, r.rt_start); o1 = near(t, r.rt_end)
    def walk(step, f, cap):
        j = ia; n_ = 0
        while True:
            n = j + step
            if n < 0 or n >= len(v) or (cap and n_ >= cap) or v[n] > v[j]:
                return j
            j = n; n_ += 1
            if v[j] <= f * v[ia]:
                return j
    row = dict(precursor_id=r.precursor_id, dn_l=d0, dn_r=d1, cur_l=o0, cur_r=o1)
    for nm, f, cap in [('valley', 0.05, 0), ('valley_cap', 0.05, 6), ('f10_cap', 0.10, 6)]:
        row[nm + '_l'] = walk(-1, 0.0 + f, cap); row[nm + '_r'] = walk(1, f, cap)
    # internal valleys of our current peak between DIA-NN stop and our end (merged shoulders)
    seg = v[d1:o1 + 1]
    row['merged_shoulder'] = o1 > d1 + 1 and any(seg[k] < seg[k - 1] and seg[k] < seg[k + 1] for k in range(1, len(seg) - 1))
    res.append(row)
R = pd.DataFrame(res)
say(f'n={len(R)} (q<=1%, DIA-NN apex inside our matched peak)')
for nm in ['cur', 'valley', 'valley_cap', 'f10_cap']:
    dl = R[nm + '_l'] - R.dn_l; dr = R[nm + '_r'] - R.dn_r
    inter = (np.minimum(R[nm + '_r'], R.dn_r) - np.maximum(R[nm + '_l'], R.dn_l) + 1).clip(lower=0)
    union = np.maximum(R[nm + '_r'], R.dn_r) - np.minimum(R[nm + '_l'], R.dn_l) + 1
    say(f'  {nm:11s} start exact {(dl == 0).mean():.3f} within1 {(dl.abs() <= 1).mean():.3f} median {dl.median():+.0f} | '
        f'stop exact {(dr == 0).mean():.3f} within1 {(dr.abs() <= 1).mean():.3f} median {dr.median():+.0f} p90 {dr.quantile(.9):+.0f} | '
        f'width median {(R[nm + "_r"] - R[nm + "_l"]).median():.0f} (DIA-NN {(R.dn_r - R.dn_l).median():.0f}) | scan IoU median {(inter / union).median():.2f}')
say(f'  our end > DIA-NN stop + 1 with an internal valley in between (merged shoulder): {R.merged_shoulder.mean():.3f}')

# q<=1% not contained: cause
nc = q1[~q1.contains].copy()
g = us.groupby('precursor_id').rt.agg(['min', 'max'])
nc = nc.join(g, on='precursor_id')
nc['dn_apex_in_window'] = (nc.RT >= nc['min']) & (nc.RT <= nc['max'])
say(f'\nq<=1% with DIA-NN apex NOT inside any of our (<=5) peaks: {len(nc)}; DIA-NN apex outside our RT window: {(~nc.dn_apex_in_window).sum()} '
    f'(entrapments {int(nc[~nc.dn_apex_in_window].e.sum())}); inside window: {nc.dn_apex_in_window.sum()} (entrapments {int(nc[nc.dn_apex_in_window].e.sum())})')
iw = nc[nc.dn_apex_in_window]
say('  inside window: our peaks_in_window', iw.peaks_in_window.describe()[['min', '50%', 'max']].to_dict(),
    '| nearest-peak apex offset scans', iw.d_apex_sc.abs().round().describe()[['min', '50%', 'max']].to_dict())
R.to_csv(X + '/rtmono/our_boundary_rules_rows_raw01.tsv', sep='\t', index=False)
