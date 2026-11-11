"""D diagnostic 2 (no training): why did cross-run features raise FDP_c? Entrapment labels are used only to measure.
Among D11FWT precursor winners, binned by own score (per-fold standardised as the evaluator), compare decoy, entrapment
and human winners on X presence (an other-run window of the same precursor within 10 s) and on the other run's best
D11FWT score there. If entrapment > decoy at matched own score, cross-run evidence is not target-decoy symmetric for
false targets, and decoy-based q-values cannot calibrate it.
Second part: the shared-signal hypothesis. For each winner, the number of other population windows (any precursor) in
the same run with apex within 2 s, precursor m/z within 1 Th, and at least 2 library fragment m/z within 10 ppm of the
winner's fragments. Reported by class and score bin."""
import sqlite3, numpy as np, pandas as pd
R = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; O = R + '/astral-voxel-pilot-20261110/v2'
c = sqlite3.connect('file:' + R + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c).set_index('precursor_id').e
pmz = pd.read_sql('select ID precursor_id, PRECURSOR_MZ mz from PRECURSOR', c).set_index('precursor_id').mz
tr = pd.read_sql('select m.PRECURSOR_ID pid, t.PRODUCT_MZ fmz from TRANSITION t join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID', c)
frag = {p: np.sort(g.fmz.values) for p, g in tr.groupby('pid')}
ARM = 'D11FWT'


def winners(pop, s):
  w = pop.assign(s=s.astype(np.float64), row=np.arange(len(pop))).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').copy()
  for f_, g in w[w.is_decoy == 1].groupby('fold'):
    med = g.s.median(); mad = 1.4826 * (g.s - med).abs().median(); w.loc[w.fold == f_, 's'] = (w.loc[w.fold == f_, 's'] - med) / mad
  return w


def shared(a, b, ppm=10.0):
  """Number of fragments of a within ppm of any fragment of b."""
  i = np.searchsorted(b, a); n = 0
  for x, k in zip(a, i):
    for j in (k - 1, k):
      if 0 <= j < len(b) and abs(b[j] - x) <= x * ppm * 1e-6: n += 1; break
  return n


for raw, other in (('raw01', 'raw02'), ('raw02', 'raw01')):
  pop = pd.read_csv(f'{O}/{raw}/population.tsv', sep='\t')
  w = winners(pop, np.load(f'{O}/{raw}/scores_{ARM}.npy'))
  x = np.load(f'{O}/{raw}/xrun_{ARM}.npy')  # other-run D11FWT score within 10 s (diag_crossrun.py), -inf if none
  w['x'] = x[w.row.values]
  w['cls'] = np.where(w.is_decoy == 1, 'decoy', np.where(w.precursor_id.map(ent) == 1, 'entrap', 'human'))
  dq = w[w.cls == 'decoy'].s.quantile([0.9, 0.99, 0.999]).values
  w['bin'] = pd.cut(w.s, [-np.inf, dq[0], dq[1], dq[2], np.inf], labels=['<q90', 'q90-q99', 'q99-q99.9', '>q99.9 decoy'])
  print(f'\n== {raw}: X presence and other-run score by own-score bin (bins = decoy-winner quantiles)')
  g = w.groupby(['bin', 'cls'], observed=True).agg(n=('s', 'size'), x_present=('x', lambda v: np.isfinite(v).mean()),
                                                     x_med=('x', lambda v: np.median(v[np.isfinite(v)]) if np.isfinite(v).any() else np.nan))
  print(g.round(3).to_string())
  # Shared-signal competitors for the winners in the top two bins (cost limits it to these).
  top = w[w.bin.isin(['q99-q99.9', '>q99.9 decoy'])].copy()
  allw = pop.assign(mz=pop.precursor_id.map(pmz).values)
  allw = allw.sort_values('rt_apex'); rts = allw.rt_apex.values; mzs = allw.mz.values; pids = allw.precursor_id.values
  cnt = []
  for p, rt in zip(top.precursor_id.values, top.rt_apex.values):
    lo, hi = np.searchsorted(rts, rt - 2.0), np.searchsorted(rts, rt + 2.0)
    m0 = pmz[p]; fa = frag[p]; n = 0; seen = set()
    for j in range(lo, hi):
      q = pids[j]
      if q == p or q in seen or abs(mzs[j] - m0) > 1.0: continue
      seen.add(q)
      if shared(fa, frag[q]) >= 2: n += 1
    cnt.append(n)
  top['competitors'] = cnt
  print(f'{raw}: shared-signal competitors (other precursors, apex +-2 s, m/z +-1 Th, >= 2 shared fragments within 10 ppm)')
  print(top.groupby(['bin', 'cls'], observed=True).competitors.agg(n='size', any_share=lambda v: (v > 0).mean(), mean='mean').round(3).to_string())
