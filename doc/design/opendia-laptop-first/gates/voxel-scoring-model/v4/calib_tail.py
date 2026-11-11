"""PLAN_V6 per-fold decoy-tail calibration (label-free, no training). Usage: calib_tail.py <raw> <src arm> <out arm>.
Per fold k: decoy precursor winners d (scores of each precursor's best window, same tie rule as the evaluator).
s' = -log S_k(s) with S_k(s) = (#{d >= s} + 1) / (n_k + 1), linear in s between decoy scores, up to t_k (top-1 % decoy
score); above t_k: s' = -log S_k(t_k) + (s - t_k) / beta_k, beta_k = mean excess of the top-1 % decoys over t_k.
Applied to every window of the fold; windows at -1e30 (W rule) stay there."""
import sys, numpy as np, pandas as pd
RAW, SRC, OUT = sys.argv[1:4]
O = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2'
pop = pd.read_csv(f'{O}/{RAW}/population.tsv', sep='\t', usecols=['sample_id', 'precursor_id', 'is_decoy', 'fold'])
s = np.load(f'{O}/{RAW}/scores_{SRC}.npy').astype(np.float64)
win = pop.assign(s=s).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id')
out = np.full(len(s), -1e30, np.float64)
for k in sorted(pop.fold.unique()):
    d = np.sort(win[(win.fold == k) & (win.is_decoy == 1) & (win.s > -1e29)].s.values)  # ascending
    n = len(d); ntop = max(n // 100, 10); t = d[n - ntop]; beta = max(float((d[n - ntop:] - t).mean()), 1e-9)
    # survival at each decoy score: S(d_i) = (n - i + 1) / (n + 1) for ascending index i (count of decoys >= d_i incl. itself)
    surv = (n - np.arange(n) + 1) / (n + 1)
    xs = d[: n - ntop + 1]; ys = -np.log(surv[: n - ntop + 1])
    xs, idx = np.unique(xs, return_index=True); ys = ys[idx]  # unique scores keep the first (largest survival)
    m = (pop.fold.values == k) & (s > -1e29)
    v = s[m]
    low = np.interp(v, xs, ys)
    # below the lowest decoy: extend linearly with the slope of the first segment to stay strictly increasing
    if len(xs) > 1:
        sl = (ys[1] - ys[0]) / max(xs[1] - xs[0], 1e-12); low = np.where(v < xs[0], ys[0] + (v - xs[0]) * sl, low)
    hi = -np.log(surv[n - ntop]) + (v - t) / beta
    out[m] = np.where(v > t, hi, low)
    print(RAW, OUT, f'fold {k}: decoys {n}, t {t:.2f}, beta {beta:.2f}, -logS(t) {-np.log(surv[n - ntop]):.2f}', flush=True)
# Checks: within-fold order of finite scores unchanged; window choice unchanged.
for k in sorted(pop.fold.unique()):
    m = (pop.fold.values == k) & (s > -1e29); o1 = np.argsort(s[m], kind='stable'); assert np.all(np.diff(out[m][o1]) >= -1e-9), 'not monotone'
np.save(f'{O}/{RAW}/scores_{OUT}.npy', out.astype(np.float32))
w2 = pop.assign(s=out).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id')
print(RAW, OUT, 'window choice unchanged:', bool((w2.set_index('precursor_id').sample_id.sort_index().values == win.set_index('precursor_id').sample_id.sort_index().values).all()))
