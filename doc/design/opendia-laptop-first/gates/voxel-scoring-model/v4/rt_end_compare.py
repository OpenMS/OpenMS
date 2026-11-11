"""RT-map end-behaviour comparison. Selection is label-free: maps are fitted on calibration anchors (anchors.tsv used=1)
and scored on held-out anchors. The test anchors are the 15 % of anchors with the lowest and the highest library RT,
held out together, so every map has to extrapolate to them. DIA-NN RT of the reference precursors outside the anchor
library-RT range is a second, reference-only check (DIA-NN is not used to choose).
Maps (all monotone):
  C   current: binned-median knots (fitRTMap), constant outside the knots
  L3  C + linear ends through end knot and 3rd knot from the end (extend_cal.py)
  LF  C + linear ends with slope from a least-squares line on the outermost 20 % of fit anchors at each end
  Q   C + quadratic ends from a least-squares quadratic on the outermost 20 % at each end, slope clamped >= 0
  PCH monotone cubic (PCHIP) through the binned knots, linear ends with the PCHIP end slope
  LO  LOWESS (frac 0.3) on all fit anchors, evaluated on a grid, linear ends with the slope of the last 5 % of the grid
  E   C + end knots from the extreme 15 anchors (rtx_maps M2 in DIAG_B_RT.txt: median of the 15 lowest / highest), then L3 ends"""
import ast, json, sqlite3, numpy as np, pandas as pd
from scipy.interpolate import PchipInterpolator


def lowess(y, x, frac=0.3, grid=200):
  """Tricube local-linear smoother on a grid over x (no robustness iterations); returns [grid, 2] sorted."""
  o = np.argsort(x); x, y = x[o], y[o]; k = max(int(frac * len(x)), 5); g = np.linspace(x[0], x[-1], grid); out = []
  for t in g:
    d = np.abs(x - t); h = np.partition(d, k - 1)[k - 1] + 1e-12; w = np.clip(1 - (d / h) ** 3, 0, None) ** 3
    W = w.sum(); mx = (w * x).sum() / W; my = (w * y).sum() / W; sxx = (w * (x - mx) ** 2).sum()
    b = (w * (x - mx) * (y - my)).sum() / sxx if sxx > 0 else 0.0; out.append(my + b * (t - mx))
  return np.c_[g, out]
R = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
CAL = {'raw01': R + '/astral-design-raw01-20261107/window300k', 'raw02': R + '/astral-entrapment-raw02-20261108/cal'}
c = sqlite3.connect('file:' + R + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
lib = pd.read_sql('select ID precursor_id, LIBRARY_RT lrt from PRECURSOR', c).set_index('precursor_id').lrt
ref = pd.read_parquet(R + '/astral-diann-gap-20261108/ref_exact.parquet')


def binned(x, y, bins=20):
  """fitRTMap: equal-count bin medians, 3 sigma rejection, up to 4 rounds."""
  o = np.argsort(x); x, y = x[o], y[o]
  for _ in range(4):
    e = np.array_split(np.arange(len(x)), bins)
    kx = np.array([np.median(x[i]) for i in e]); ky = np.array([np.median(y[i]) for i in e])
    r = np.abs(y - np.interp(x, kx, ky)); s = 1.4826 * np.median(r); keep = r <= max(3 * s, 3)
    if keep.all(): break
    x, y = x[keep], y[keep]
  return kx, ky, x, y


def with_ends(kx, ky, sl, sh):
  return lambda t: np.where(t < kx[0], ky[0] + sl * (t - kx[0]), np.where(t > kx[-1], ky[-1] + sh * (t - kx[-1]), np.interp(t, kx, ky)))


def maps(x, y):
  kx, ky, xs, ys = binned(x, y); out = {'C': lambda t: np.interp(t, kx, ky)}
  out['L3'] = with_ends(kx, ky, (ky[3] - ky[0]) / (kx[3] - kx[0]), (ky[-1] - ky[-4]) / (kx[-1] - kx[-4]))
  n = max(int(0.2 * len(xs)), 10); lo, hi = slice(0, n), slice(len(xs) - n, len(xs))
  sl = np.polyfit(xs[lo], ys[lo], 1)[0]; sh = np.polyfit(xs[hi], ys[hi], 1)[0]
  out['LF'] = with_ends(kx, ky, max(sl, 0), max(sh, 0))
  ql, qh = np.polyfit(xs[lo], ys[lo], 2), np.polyfit(xs[hi], ys[hi], 2)
  def quad(t):
    base = np.interp(t, kx, ky)
    # anchor each quadratic at the end knot, clamp to non-decreasing
    fl = lambda u: ky[0] + np.polyval(ql, u) - np.polyval(ql, kx[0]); fh = lambda u: ky[-1] + np.polyval(qh, u) - np.polyval(qh, kx[-1])
    vl = np.minimum(fl(t), ky[0]); vh = np.maximum(fh(t), ky[-1])
    return np.where(t < kx[0], vl, np.where(t > kx[-1], vh, base))
  out['Q'] = quad
  p = PchipInterpolator(kx, ky, extrapolate=False); d = p.derivative()
  out['PCH'] = lambda t: np.where(t < kx[0], ky[0] + d(kx[0]) * (t - kx[0]), np.where(t > kx[-1], ky[-1] + d(kx[-1]) * (t - kx[-1]), p(np.clip(t, kx[0], kx[-1]))))
  sm = lowess(ys, xs, frac=0.3); gx, gy = sm[:, 0], np.maximum.accumulate(sm[:, 1])
  k = max(int(0.05 * len(gx)), 3)
  sll = (gy[k] - gy[0]) / max(gx[k] - gx[0], 1e-9); slh = (gy[-1] - gy[-1 - k]) / max(gx[-1] - gx[-1 - k], 1e-9)
  out['LO'] = lambda t: np.where(t < gx[0], gy[0] + sll * (t - gx[0]), np.where(t > gx[-1], gy[-1] + slh * (t - gx[-1]), np.interp(t, gx, gy)))
  ekx = np.r_[np.median(xs[:15]), kx, np.median(xs[-15:])]; eky = np.r_[np.median(ys[:15]), ky, np.median(ys[-15:])]
  ekx, idx = np.unique(ekx, return_index=True); eky = np.maximum.accumulate(eky[idx])
  out['E'] = with_ends(ekx, eky, (eky[3] - eky[0]) / (ekx[3] - ekx[0]), (eky[-1] - eky[-4]) / (ekx[-1] - ekx[-4]))
  return out, (xs.min(), xs.max())


for raw in ('raw01', 'raw02'):
  an = pd.read_csv(CAL[raw] + '/anchors.tsv', sep='\t').query('used == 1').sort_values('library_rt')
  x, y = an.library_rt.values, an.run_rt.values; n = len(x); h = int(0.15 * n)
  test = np.r_[np.arange(h), np.arange(n - h, n)]; fit = np.arange(h, n - h)
  M, (lo, hi) = maps(x[fit], y[fit])
  Mall, (alo, ahi) = maps(x, y)
  r = ref[ref.Run.str.endswith('_' + raw[-2:]) & ~ref.entrap].drop_duplicates('precursor_id'); rl = r.precursor_id.map(lib).values; ro = r.RT.values * 60
  outside = (rl < alo) | (rl > ahi)
  sig = json.load(open(CAL[raw] + '/cal.json'))['rt']['sigma_s']; tol = 4 * sig
  print(f'\n== {raw}: anchors {n}, fit on middle {len(fit)} (library RT {lo:.3f}..{hi:.3f}), test {len(test)} end anchors; tolerance {tol:.0f} s')
  print(f'{"map":4s} | held-out end anchors: median |res| s, q90, within tol | DIA-NN ref outside full-anchor range (n {outside.sum()}): median |res|, within tol (maps fitted on all anchors)')
  for k in M:
    e = np.abs(y[test] - M[k](x[test])); f = np.abs(ro[outside] - Mall[k](rl[outside]))
    print(f'{k:4s} | {np.median(e):6.1f} {np.quantile(e, 0.9):6.1f} {np.mean(e <= tol):.3f} | {np.median(f):6.1f} {np.mean(f <= tol):.3f}')


# Part 2: distance-dependent tolerance, label-free. On the held-out end anchors (L3 map fitted on the middle), regress
# the q95 of |residual| on the library-RT distance beyond the fit range: tol(d) = tol0 + g * d. Then on the full map,
# apply tol(d) with d = distance beyond the full-anchor library-RT range, and report the DIA-NN reference coverage
# (reference only) and the share of library precursors whose window widens (cost, same for targets and decoys).
print('\n== Part 2: distance-dependent RT tolerance (g fitted on held-out end anchors)')
for raw in ('raw01', 'raw02'):
  an = pd.read_csv(CAL[raw] + '/anchors.tsv', sep='\t').query('used == 1').sort_values('library_rt')
  x, y = an.library_rt.values, an.run_rt.values; n = len(x); h = int(0.15 * n)
  test = np.r_[np.arange(h), np.arange(n - h, n)]; fit = np.arange(h, n - h)
  M, (lo, hi) = maps(x[fit], y[fit]); sig = json.load(open(CAL[raw] + '/cal.json'))['rt']['sigma_s']; tol0 = 4 * sig
  dt = np.maximum(lo - x[test], x[test] - hi); e = np.abs(y[test] - M['L3'](x[test]))
  # g: smallest slope so that 95 % of held-out end anchors fall within tol0 + g d
  g = np.quantile(np.maximum(e - tol0, 0) / np.maximum(dt, 1e-3), 0.95)
  Mall, (alo, ahi) = maps(x, y)
  r = ref[ref.Run.str.endswith('_' + raw[-2:]) & ~ref.entrap].drop_duplicates('precursor_id'); rl = r.precursor_id.map(lib).values; ro = r.RT.values * 60
  d = np.maximum(np.maximum(alo - rl, rl - ahi), 0); res = np.abs(ro - Mall['L3'](rl)); out = d > 0
  dl = np.maximum(np.maximum(alo - lib.values, lib.values - ahi), 0)
  for name, tol in (('fixed 4 sigma', np.full_like(d, tol0)), (f'4 sigma + {g:.0f} s/unit x distance', tol0 + g * d)):
    tl = tol0 + (g * dl if 'distance' in name else 0 * dl)
    print(f'{raw} {name:36s}: DIA-NN ref outside anchors within tol {np.mean(res[out] <= tol[out]):.3f} (n {out.sum()}), all ref {np.mean(res <= tol):.3f}; '
          f'library mean half-window {tl.mean():.0f} s, max {tl.max():.0f} s')
