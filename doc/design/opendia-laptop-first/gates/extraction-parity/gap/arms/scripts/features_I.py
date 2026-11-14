"""PLAN_GAP arm I: explaining away shared fragment signal (step 2 of DESIGN_INTERFERENCE_GRAPH.md, a one-pass fit, not
the full belief propagation). Label-free and symmetric for targets and decoys. featI.npy [n, 4] float32.
Graph (as interference_components.py, on the scored population with the 12-fragment library):
  windows a, b of different precursors are linked when their precursors are in the same DIA isolation window, their
  peak ranges [rt_start, rt_end] overlap, and >= 1 fragment of a matches one of b within the run's fitted fragment
  tolerance (cal.json ms2_mz.tolerance_half_width_ppm). Per (a, claimant precursor) only the claimant window with the
  nearest apex is kept (two windows of one precursor see the same signal).
Fit per directed link (a <- b), from the 12-fragment arrays (xic12 [n,12,32], lib12 [n,12]):
  - b's abundance: scale_b = median over b's fragments NOT shared with a (lib > 0, central area > 0) of area / lib,
    central bins 12:20 of b's own window; no unshared fragment with signal -> b explains nothing;
  - b's elution in a's window: frac = (b's summed profile over the bins of a's central range, shifted by the apex
    difference) / (b's summed profile over its own central bins);
  - b's expected area on a's shared fragment f (matched to b's fragment g): scale_b * lib_b[g] * frac.
  Contributions of all claimants are summed per fragment of a and capped at a's observed area.
Columns:
  i_explained     explained share of a's central area (0 without signal)
  i_libcos_resid  cosine of sqrt(residual area) and sqrt(lib12) over the 12 fragments
  i_dlibcos       i_libcos_resid minus the same cosine on the unexplained-away areas
  i_nclaim        log(1 + claimant precursors with a nonzero contribution)
Usage: features_I.py <raw01|raw02> <root> <library12_pop.pqp> <iso_windows.tsv> <search.features.tsv>."""
import json, sqlite3, sys, time, warnings, numpy as np, pandas as pd
RAW, O, LIB12, ISO, FEAT = sys.argv[1:6]
C0, C1, K, BIN, CH = 12, 20, 12, 30.0 / 32, 2_000_000
warnings.filterwarnings('ignore', category=RuntimeWarning)
t0 = time.time()
TOL = float(json.load(open(f'{O}/{RAW}/cal/cal.json'))['ms2_mz']['tolerance_half_width_ppm']) * 1e-6
iso = pd.read_csv(ISO, sep='\t')[['lo', 'hi']].sort_values('lo').values
pop = pd.read_csv(f'{O}/{RAW}/population.tsv', sep='\t', usecols=['sample_id', 'precursor_id', 'rt_apex'])
n = len(pop)
f = pd.read_csv(FEAT, sep='\t', usecols=['candidate_row', 'rt_start', 'rt_end']).set_index('candidate_row')
pop['rt_start'] = f.rt_start.reindex(pop.sample_id).values; pop['rt_end'] = f.rt_end.reindex(pop.sample_id).values
assert pop.rt_start.notna().all()
c = sqlite3.connect(f'file:{LIB12}?mode=ro', uri=True)
pmz = pd.read_sql('select ID pid, PRECURSOR_MZ mz from PRECURSOR', c).set_index('pid').mz
t = pd.read_sql("""select m.PRECURSOR_ID pid, t.PRODUCT_MZ mz from TRANSITION t join TRANSITION_PRECURSOR_MAPPING m
  on m.TRANSITION_ID=t.ID order by m.PRECURSOR_ID, t.ID""", c)
t['slot'] = t.groupby('pid').cumcount(); t = t[t.slot < K]
uniq, inv = np.unique(pop.precursor_id.values, return_inverse=True)
FMu = np.full((len(uniq), K), np.nan)
i = np.searchsorted(uniq, t.pid.values); ok = (i < len(uniq)) & (uniq[np.minimum(i, len(uniq) - 1)] == t.pid.values)
FMu[i[ok], t.slot.values[ok]] = t.mz.values[ok]
mz = pmz.reindex(pop.precursor_id).values
k = np.searchsorted(iso[:, 0], mz, side='right') - 1
pop['iso'] = np.where((k >= 0) & (mz < iso[np.clip(k, 0, None), 1]), k, -1)

# Per-window central areas, library and cumulative summed profile.
X = np.load(f'{O}/{RAW}/xic12.npy', mmap_mode='r'); L = np.asarray(np.load(f'{O}/{RAW}/lib12.npy', mmap_mode='r'), np.float64)
A = np.zeros((n, K)); CS = np.zeros((n, 33))
for b in range(0, n, 100000):
    x = np.asarray(X[b:b + 100000], np.float64); A[b:b + 100000] = x[:, :, C0:C1].sum(2)
    CS[b:b + 100000, 1:] = np.cumsum(x.sum(1), 1)
print(f'{RAW}: {n} windows, {len(iso)} isolation windows, tolerance {TOL * 1e6:.2f} ppm, arrays read ({time.time() - t0:.0f} s)', flush=True)

# Links.
pid = pop.precursor_id.values; apex = pop.rt_apex.values
ra, rb, ma, mb = [], [], [], []  # directed: a <- b with slot maps a->b and b->a
for _, g in pop[pop.iso >= 0].groupby('iso'):
    g = g.sort_values('rt_start'); idx = g.index.values; st, en = g.rt_start.values, g.rt_end.values
    last = np.searchsorted(st, en, side='right'); cnt = np.maximum(last - np.arange(len(g)) - 1, 0)
    I = np.repeat(np.arange(len(g)), cnt); J = I + 1 + (np.arange(len(I)) - np.repeat(np.cumsum(cnt) - cnt, cnt))
    I, J = idx[I], idx[J]; keep = pid[I] != pid[J]; I, J = I[keep], J[keep]
    for s0 in range(0, len(I), CH):
        a, b = I[s0:s0 + CH], J[s0:s0 + CH]
        fa, fb = FMu[inv[a]][:, :, None], FMu[inv[b]][:, None, :]
        M = np.abs(fa - fb) <= TOL * fa  # [pairs, slot a, slot b]
        hit = M.any((1, 2))
        if not hit.any(): continue
        M = M[hit]; a, b = a[hit], b[hit]
        m_ab = np.where(M.any(2), M.argmax(2), -1).astype(np.int8)  # a slot -> b slot
        m_ba = np.where(M.any(1), M.argmax(1), -1).astype(np.int8)  # b slot -> a slot
        ra += [a, b]; rb += [b, a]; ma += [m_ab, m_ba]; mb += [m_ba, m_ab]
if not ra:  # no window shares any fragment signal (post-run fix; never reached on the real runs)
    ra, rb = [np.zeros(0, int)], [np.zeros(0, int)]; ma, mb = [np.zeros((0, K), np.int8)], [np.zeros((0, K), np.int8)]
ra, rb = np.concatenate(ra), np.concatenate(rb); ma, mb = np.concatenate(ma), np.concatenate(mb)
print(f'directed links {len(ra)} ({time.time() - t0:.0f} s)', flush=True)
# One claimant window per (a, claimant precursor): nearest apex.
d = np.abs(apex[ra] - apex[rb])
o = np.lexsort((d, pid[rb], ra)); ra, rb, ma, mb, d = ra[o], rb[o], ma[o], mb[o], d[o]
first = np.r_[True, (ra[1:] != ra[:-1]) | (pid[rb][1:] != pid[rb][:-1])][:len(ra)]
ra, rb, ma, mb = ra[first], rb[first], ma[first], mb[first]
print(f'links after one window per claimant precursor {len(ra)}; windows with a claimant {len(np.unique(ra))} of {n}', flush=True)

contrib = np.zeros((n, K)); ncl = np.zeros(n)
for s0 in range(0, len(ra), CH):
    a, b, m_ab, m_ba = ra[s0:s0 + CH], rb[s0:s0 + CH], ma[s0:s0 + CH].astype(int), mb[s0:s0 + CH].astype(int)
    Ab, Lb = A[b], L[b]
    shared_b = m_ba >= 0  # b's slots matched to a fragment of a
    use = (~shared_b) & (Lb > 0) & (Ab > 0)
    scale = np.where(use.any(1), np.nanmedian(np.where(use, Ab / np.maximum(Lb, 1e-12), np.nan), 1), 0.0)
    s = np.rint((apex[a] - apex[b]) / BIN).astype(int)
    lo, hi = np.clip(C0 + s, 0, 32), np.clip(C1 + s, 0, 32)
    own = CS[b, C1] - CS[b, C0]; r_ = np.arange(len(b))
    frac = np.where(own > 0, (CS[b, hi] - CS[b, lo]) / np.maximum(own, 1e-12), 0.0)
    e = np.where(m_ab >= 0, scale[:, None] * Lb[r_[:, None], np.maximum(m_ab, 0)] * frac[:, None], 0.0)
    np.add.at(contrib, a, e)
    np.add.at(ncl, a, (e.sum(1) > 0).astype(float))
expl = np.minimum(A, contrib); res = A - expl
tot = A.sum(1)
def cos(x): return (np.sqrt(x) * np.sqrt(L)).sum(1) / np.sqrt(x.sum(1) * L.sum(1) + 1e-12)
c_res, c_obs = cos(res), cos(A)
out = np.c_[np.where(tot > 0, expl.sum(1) / np.maximum(tot, 1e-12), 0.0), c_res, c_res - c_obs, np.log1p(ncl)].astype(np.float32)
np.save(f'{O}/{RAW}/featI.npy', out)
print(RAW, 'featI', out.shape, 'finite', bool(np.isfinite(out).all()), 'share explained > 0', round(float((out[:, 0] > 0).mean()), 3),
      'p50/p90/p99', np.percentile(out, [50, 90, 99], axis=0).round(3).T.tolist(), f'({time.time() - t0:.0f} s)')
