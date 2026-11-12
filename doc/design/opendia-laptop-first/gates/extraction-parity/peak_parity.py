"""Peak parity: OpenDIA xic peaks vs DIA-NN 2.2.0 peaks, raw01, same precursors.

For every DIA-NN precursor (report of the --xic run, q <= 0.5) take our target candidate peaks (search.features.tsv,
up to xic_max_peaks per precursor) and compare with DIA-NN's RT (apex), RT.Start, RT.Stop and FWHM.
Positions are converted to scan indices on our own native-scan grid for that precursor (traces.tsv), so offsets
are reported in scans (cycle ~1.65 s). Labels (entrapment) are used only to split the evaluation.
Usage: peak_parity.py <run dir> <outfile>
"""
import sys, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; X = B + '/extraction-parity-20261112'
RUN = sys.argv[1] if len(sys.argv) > 1 else 'rtmono'
OUT = sys.argv[2] if len(sys.argv) > 2 else RUN + '/PEAK_PARITY_raw01.txt'
out = open(X + '/' + OUT, 'w')
def say(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')

idm = pd.read_csv(X + '/idmap.tsv', sep='\t')
rep = pd.read_parquet(X + '/diann/report_raw01.parquet', columns=['Precursor.Id', 'Q.Value', 'RT', 'RT.Start', 'RT.Stop', 'FWHM'])
rep = rep.merge(idm[['Precursor.Id', 'precursor_id']], on='Precursor.Id')
ent = pd.read_csv(B + '/diann-ourlib-20261112/diann_ourlib_raw01.tsv.gz', sep='\t', usecols=['precursor_id', 'e']).drop_duplicates('precursor_id')
rep = rep.merge(ent, on='precursor_id', how='left')
for c in ['RT', 'RT.Start', 'RT.Stop', 'FWHM']:
    rep[c] = rep[c] * 60.0
rep['grp'] = pd.cut(rep['Q.Value'], [-1, 0.01, 0.05, 1], labels=['q<=1%', '1-5%', '>5%'])
ids = set(rep.precursor_id)

cols = ['precursor_id', 'is_decoy', 'rt_start', 'rt_apex', 'rt_end', 'fwhm_seconds', 'coelution', 'peak_rank', 'peaks_in_window', 'apex_intensity']
parts = []
for ch in pd.read_csv(X + '/' + RUN + '/search.features.tsv', sep='\t', usecols=cols, chunksize=2_000_000):
    parts.append(ch[(ch.is_decoy == 0) & ch.precursor_id.isin(ids)])
f = pd.concat(parts)

# Native-scan grid per precursor (scan RTs of the dumped traces).
t = pd.read_csv(X + '/' + RUN + '/traces.tsv', sep='\t', usecols=['precursor_id', 'rt', 'slot'])
t = t[t.slot == 0]
grid = {p: np.sort(g.rt.values) for p, g in t.groupby('precursor_id')}
def idx(p, rt):
    g = grid.get(p)
    if g is None or len(g) < 2:
        return np.nan
    # fractional scan index (linear between scans), extrapolated by the median cycle
    cyc = np.median(np.diff(g))
    if rt < g[0]:
        return (rt - g[0]) / cyc
    if rt > g[-1]:
        return len(g) - 1 + (rt - g[-1]) / cyc
    return np.interp(rt, g, np.arange(len(g)))

# DIA-NN boundaries on the grid: are they scan RTs?
snap = []
for r in rep.itertuples():
    g = grid.get(r.precursor_id)
    if g is None:
        continue
    for v in (r.RT, r._4, r._5):
        snap.append(np.min(np.abs(g - v)))
snap = np.array(snap)
say(f'run {RUN}: DIA-NN precursors {len(rep)}, with our trace grid {sum(p in grid for p in rep.precursor_id)}, with our target peaks {f.precursor_id.nunique()}')
say('DIA-NN RT / RT.Start / RT.Stop distance to nearest of our scans, s quantiles', pd.Series(snap).quantile([.5, .9, .99, 1]).round(4).to_dict())

f = f.merge(rep[['precursor_id', 'RT', 'RT.Start', 'RT.Stop', 'FWHM', 'grp', 'e']], on='precursor_id')
f['dapex'] = f.rt_apex - f.RT
f['inter'] = (np.minimum(f.rt_end, f['RT.Stop']) - np.maximum(f.rt_start, f['RT.Start'])).clip(lower=0)
f['union'] = np.maximum(f.rt_end, f['RT.Stop']) - np.minimum(f.rt_start, f['RT.Start'])
f['iou'] = f.inter / f.union
f['contains'] = (f.rt_start <= f.RT) & (f.rt_end >= f.RT)
# Our peak "is DIA-NN's peak" if its boundaries contain DIA-NN's apex (else nearest apex).
f['key'] = np.where(f.contains, 0, 1) * 1e6 + f.dapex.abs()
best = f.sort_values('key').groupby('precursor_id').head(1).set_index('precursor_id')
top = f[f.peak_rank == 0].set_index('precursor_id')

def scans(p, a, b):
    return idx(p, a) - idx(p, b)
best['d_apex_sc'] = [scans(p, a, b) for p, a, b in zip(best.index, best.rt_apex, best.RT)]
best['d_start_sc'] = [scans(p, a, b) for p, a, b in zip(best.index, best.rt_start, best['RT.Start'])]
best['d_stop_sc'] = [scans(p, a, b) for p, a, b in zip(best.index, best.rt_end, best['RT.Stop'])]
best['w_ours'] = best.rt_end - best.rt_start; best['w_dn'] = best['RT.Stop'] - best['RT.Start']
best['dn_w_sc'] = [scans(p, a, b) for p, a, b in zip(best.index, best['RT.Stop'], best['RT.Start'])]
best['our_w_sc'] = [scans(p, a, b) for p, a, b in zip(best.index, best.rt_end, best.rt_start)]
best['is_top'] = [p in top.index and abs(top.loc[p, 'rt_apex'] - best.loc[p, 'rt_apex']) < 1e-6 for p in best.index]

def q(s, qs=(.1, .25, .5, .75, .9)):
    return {k: round(float(v), 2) for k, v in s.quantile(list(qs)).items()}
def report(name, b, allrep):
    n = len(allrep)
    say(f'\n== {name}: DIA-NN precursors {n}, with any of our target peaks {len(b)} ==')
    say(f'  DIA-NN apex inside one of our peaks: {b.contains.sum()} ({b.contains.sum() / n:.4f}); that peak is our rank-0: {(b.contains & b.is_top).sum()} ({(b.contains & b.is_top).sum() / max(1, b.contains.sum()):.4f} of contained)')
    a = b.d_apex_sc.round().abs()
    say('  apex offset of matched peak, scans: ' + ' | '.join(f'{k}: {int(v)}' for k, v in
        {'0': (a == 0).sum(), '1': (a == 1).sum(), '2': (a == 2).sum(), '3-6': ((a >= 3) & (a <= 6)).sum(), '>6': (a > 6).sum()}.items()))
    say('  signed apex offset (ours - DIA-NN), scans, quantiles', q(b.d_apex_sc))
    c = b[b.contains]
    say(f'  among peaks containing the DIA-NN apex (n={len(c)}):')
    say('    start offset (ours - DIA-NN), scans', q(c.d_start_sc), ' stop offset', q(c.d_stop_sc))
    say('    width scans: ours', q(c.our_w_sc), ' DIA-NN', q(c.dn_w_sc))
    say('    width ratio ours/DIA-NN', q(c.w_ours / c.w_dn), ' IoU', q(c.iou))
    say(f'    our peak contains DIA-NN [start, stop]: {((c.rt_start <= c["RT.Start"] + 1e-6) & (c.rt_end >= c["RT.Stop"] - 1e-6)).mean():.4f}; '
        f'DIA-NN contains ours: {((c.rt_start >= c["RT.Start"] - 1e-6) & (c.rt_end <= c["RT.Stop"] + 1e-6)).mean():.4f}')
    say('    FWHM s: ours', q(c.fwhm_seconds), ' DIA-NN', q(c.FWHM), ' ratio', q(c.fwhm_seconds / c.FWHM))
    say('    our peaks_in_window', q(c.peaks_in_window), ' coelution', q(c.coelution))
    return c

for g in ['q<=1%', '1-5%', '>5%']:
    sub = rep[rep.grp == g]
    b = best[best.grp == g]
    report(g, b, sub)
    if g == 'q<=1%':
        for e, nm in [(0, 'q<=1% non-entrapment'), (1, 'q<=1% entrapment')]:
            report(nm, b[b.e == e], sub[sub.e == e])
best.reset_index()[['precursor_id', 'grp', 'e', 'RT', 'RT.Start', 'RT.Stop', 'FWHM', 'rt_apex', 'rt_start', 'rt_end', 'fwhm_seconds',
                    'peak_rank', 'peaks_in_window', 'coelution', 'contains', 'is_top', 'd_apex_sc', 'd_start_sc', 'd_stop_sc', 'iou']].to_csv(
    X + '/' + RUN + '/peak_parity_rows_raw01.tsv', sep='\t', index=False)
