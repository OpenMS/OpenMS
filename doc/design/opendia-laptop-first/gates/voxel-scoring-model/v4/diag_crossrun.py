"""D diagnostic (no training, no new runs): is cross-run evidence a separating signal for stage-D misses?
Label-free feature per window w of precursor p in run a: x_other(w) = the best D11FWT score of any window of the same
precursor p in the other run b whose RT, aligned to run a, is within 10 s of w's apex; -inf (no such window) otherwise.
Decoys get the same feature from their own windows in the other run, so the feature is symmetric for target-decoy
competition. RT alignment b -> a: binned medians of the apex RT of precursors that are TDC q <= 0.01 target winners in
both runs (decoy competition only, no entrapment labels).
Reading (entrapment and DIA-NN labels used here only to measure): AUC of the right windows of stage-D misses against
the top decoy winners (decoy winners above the 99th percentile of decoy winner scores), compared with the AUC of the
D11FWT score itself on the same sets."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2')
O = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2'
ARM = sys.argv[1] if len(sys.argv) > 1 else 'D11FWT'


def winners(pop, s):
  w = pop.assign(s=s).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id')
  best = w.sort_values(['s', 'is_decoy'], ascending=[False, False]).drop_duplicates('pair_key').sort_values('s', ascending=False)
  T = np.cumsum(1 - best.is_decoy.values); D = np.cumsum(best.is_decoy.values)
  best['q'] = np.minimum.accumulate(np.minimum(1, (D + 1) / np.maximum(T, 1))[::-1])[::-1]
  return w, best


def auc(a, b):
  a = np.asarray(a, float); b = np.asarray(b, float)
  allv = np.concatenate([a, b]); r = pd.Series(allv).rank().values
  return (r[:len(a)].sum() - len(a) * (len(a) + 1) / 2) / (len(a) * len(b))


pops, scores, W, Bt = {}, {}, {}, {}
for raw in ('raw01', 'raw02'):
  pops[raw] = pd.read_csv(f'{O}/{raw}/population.tsv', sep='\t')
  scores[raw] = np.load(f'{O}/{raw}/scores_{ARM}.npy').astype(np.float64)
  W[raw], Bt[raw] = winners(pops[raw], scores[raw])
# RT alignment from confident target winners in both runs.
conf = {r: Bt[r][(Bt[r].q <= 0.01) & (Bt[r].is_decoy == 0)].set_index('precursor_id').rt_apex for r in Bt}
both = conf['raw01'].to_frame('a').join(conf['raw02'].rename('b'), how='inner')
print(f'alignment anchors (q <= 0.01 target winners in both runs): {len(both)}; median rt02 - rt01 {np.median(both.b - both.a):.2f} s')
def binned_map(x, y, bins=25):
  """Binned-median knots of y over x (x sorted into equal-count bins); np.interp between knots."""
  o = np.argsort(x); x, y = x[o], y[o]; e = np.array_split(np.arange(len(x)), bins)
  kx = np.array([np.median(x[i]) for i in e]); ky = np.array([np.median(y[i]) for i in e])
  return lambda t: np.interp(t, kx, ky)


mapto = {('raw02', 'raw01'): binned_map(both.b.values, both.a.values), ('raw01', 'raw02'): binned_map(both.a.values, both.b.values)}
res = both.a - mapto[('raw02', 'raw01')](both.b.values)
print(f'alignment residual |rt01 - map(rt02)|: median {res.abs().median():.2f} s, q95 {res.abs().quantile(0.95):.2f} s')


def other_feature(a, b, tol=10.0):
  """For each window of run a: best score in run b of the same precursor within tol s after mapping b -> a."""
  pa, pb = pops[a].assign(row=np.arange(len(pops[a]))), pops[b].assign(s=scores[b])
  pb = pb[pb.s > -1e29]  # W rule: only each precursor's chosen window in run b is finite
  pb = pb.assign(rt_a=mapto[(b, a)](pb.rt_apex.values))[['precursor_id', 'rt_a', 's']]
  m = pa[['row', 'precursor_id', 'rt_apex']].merge(pb, on='precursor_id', how='inner')
  m = m[(m.rt_a - m.rt_apex).abs() <= tol]
  f = m.groupby('row').s.max()
  out = np.full(len(pa), -np.inf); out[f.index.values] = f.values
  return out


for a, b in (('raw01', 'raw02'), ('raw02', 'raw01')):
  x = other_feature(a, b)
  np.save(f'{O}/{a}/xrun_{ARM}.npy', x.astype(np.float32))
  pop = pops[a]; w = W[a].copy(); w['x'] = x[w.index.values]  # W keeps the population row index
  st = pd.read_csv(f'{O}/v3seq/misses_{ARM}_{a}.tsv', sep='\t')
  # right window of stage D: the winner (W rule) is the right window by construction of stage D
  D = w[w.precursor_id.isin(st[st.stage == 'D right window, below cut'].precursor_id)]
  E = w[w.precursor_id.isin(st[st.stage == 'E accepted'].precursor_id)]
  dec = w[(w.is_decoy == 1) & (w.s > -1e29)]; top = dec[dec.s >= dec.s.quantile(0.99)]
  fin = lambda v: np.where(np.isfinite(v), v, -50.0)
  print(f'\n== {a} (other run {b}), arm {ARM}')
  for name, g in (('E accepted', E), ('D below cut', D), ('top-1 % decoy winners', top), ('all decoy winners', dec)):
    print(f'{name:22s} n {len(g):6d}: has other-run window within 10 s {np.mean(np.isfinite(g.x)):.3f}; '
          f'other-run score median {np.median(fin(g.x.values)):.2f}')
  # own + other: sum of own score and the other-run score (both on the calibrated D11FWT scale, missing = 0 extra)
  comb = lambda g: g.s.values + np.where(np.isfinite(g.x.values), np.maximum(g.x.values, 0.0), 0.0)
  print(f'AUC D vs top-1 % decoys: own score {auc(D.s, top.s):.3f}; other-run score {auc(fin(D.x.values), fin(top.x.values)):.3f}; '
        f'own + max(other, 0) {auc(comb(D), comb(top)):.3f}')
