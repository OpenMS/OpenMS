"""D diagnostic 3 (no training; entrapment and DIA-NN used only to measure).
(a) Entrapment construction: share of entrapment / decoy / human target precursors with another library precursor of a
    different class at the same charge and precursor m/z within 5 ppm (same composition, e.g. a shuffled peptide).
(b) Signal ownership: for each D11FWT precursor winner, 'explained' = another precursor whose winner apex is within 2 s,
    precursor m/z within 1 Th (proxy for the same isolation window), sharing >= 2 library fragments within 10 ppm, AND
    with a higher own score. Label-free by construction (uses only scores, library and RT). Reported by class and
    own-score bin, and for the stage-D misses and accepted (E) reference precursors."""
import sqlite3, numpy as np, pandas as pd
R = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; O = R + '/astral-voxel-pilot-20261110/v2'
c = sqlite3.connect('file:' + R + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c).set_index('precursor_id').e
pre = pd.read_sql('select ID precursor_id, PRECURSOR_MZ mz, CHARGE z, DECOY decoy from PRECURSOR', c).set_index('precursor_id')
pre['cls'] = np.where(pre.decoy == 1, 'decoy', np.where(pre.index.map(ent) == 1, 'entrap', 'human'))
# (a)
print('== (a) same charge, precursor m/z within 5 ppm of a precursor of another class')
s = pre.sort_values('mz'); mz = s.mz.values; z = s.z.values; cl = s.cls.values; hit = {k: np.zeros(len(s), bool) for k in ('decoy', 'entrap', 'human')}
lo = np.searchsorted(mz, mz * (1 - 5e-6)); hi = np.searchsorted(mz, mz * (1 + 5e-6), 'right')
for i in range(len(s)):
  for j in range(lo[i], hi[i]):
    if j != i and z[j] == z[i]: hit[cl[j]][i] = True
for k in ('human', 'entrap', 'decoy'):
  m = cl == k; print(f'{k:6s} n {m.sum():6d}: has same-mass human {hit["human"][m].mean():.3f}, entrap {hit["entrap"][m].mean():.3f}, decoy {hit["decoy"][m].mean():.3f}')
# (b)
tr = pd.read_sql('select m.PRECURSOR_ID pid, t.PRODUCT_MZ fmz from TRANSITION t join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID', c)
frag = {p: np.sort(g.fmz.values) for p, g in tr.groupby('pid')}


def nshared(a, b, ppm=10.0):
  i = np.searchsorted(b, a); n = 0
  for x, k in zip(a, i):
    if (k < len(b) and abs(b[k] - x) <= x * ppm * 1e-6) or (k > 0 and abs(b[k - 1] - x) <= x * ppm * 1e-6): n += 1
  return n


for raw in ('raw01', 'raw02'):
  pop = pd.read_csv(f'{O}/{raw}/population.tsv', sep='\t'); s_ = np.load(f'{O}/{raw}/scores_D11FWT.npy').astype(np.float64)
  w = pop.assign(s=s_).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').copy()
  for f_, g in w[w.is_decoy == 1].groupby('fold'):
    med = g.s.median(); mad = 1.4826 * (g.s - med).abs().median(); w.loc[w.fold == f_, 's'] = (w.loc[w.fold == f_, 's'] - med) / mad
  w = w[w.s > -1e20]; w['mz'] = w.precursor_id.map(pre.mz); w['cls'] = w.precursor_id.map(pre.cls)
  dq = w[w.cls == 'decoy'].s.quantile([0.99, 0.999]).values
  st = pd.read_csv(f'{O}/v3seq/misses_D11FWT_{raw}.tsv', sep='\t').set_index('precursor_id').stage
  w['stage'] = w.precursor_id.map(st)
  sel = w[(w.s >= dq[0]) | w.stage.isin(['D right window, below cut', 'E accepted'])]
  ws = w.sort_values('rt_apex'); rts, mzs, pids, ss = ws.rt_apex.values, ws.mz.values, ws.precursor_id.values, ws.s.values
  expl = []
  for p, rt, m0, sc in zip(sel.precursor_id.values, sel.rt_apex.values, sel.mz.values, sel.s.values):
    a, b = np.searchsorted(rts, rt - 2.0), np.searchsorted(rts, rt + 2.0); e = False
    for j in range(a, b):
      if ss[j] > sc and pids[j] != p and abs(mzs[j] - m0) <= 1.0 and nshared(frag[p], frag[pids[j]]) >= 2: e = True; break
    expl.append(e)
  sel = sel.assign(explained=expl)
  sel['bin'] = np.where(sel.s >= dq[1], '>q99.9 decoy', np.where(sel.s >= dq[0], 'q99-q99.9', '<q99'))
  print(f'\n== (b) {raw}: explained by a higher-scoring winner sharing >= 2 fragments (apex +-2 s, m/z +-1 Th)')
  print(sel.groupby(['bin', 'cls']).explained.agg(n='size', share='mean').round(3).to_string())
  print(sel[sel.stage.notna()].groupby('stage').explained.agg(n='size', share='mean').round(3).to_string())
