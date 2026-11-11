"""D diagnostic 4 (not a gate; existing D11FWT scores; entrapment used only to measure).
Pair-level match-between-runs style search-space selection, the cross-run use that keeps target-decoy competition
symmetric by construction: keep a target-decoy pair in run a only if its pair winner (best of target and decoy) in the
other run b reaches paired q <= Q_SEL. Then rerun the evaluator's TDC on run a for the kept pairs only and report N at
FDP_c and FDP_c at q <= 0.01, against no selection. A pair can pass through its decoy, so decoys are selected by the
same rule as targets."""
import sqlite3, numpy as np, pandas as pd
R = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; O = R + '/astral-voxel-pilot-20261110/v2'; RR = 148983 / 151239
c = sqlite3.connect('file:' + R + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c).set_index('precursor_id').e


def tdc(pop, s, keep=None):
  w = pop.assign(s=s.astype(np.float64)).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').copy()
  if keep is not None: w = w[w.pair_key.isin(keep)]
  for f_, g in w[w.is_decoy == 1].groupby('fold'):
    med = g.s.median(); mad = 1.4826 * (g.s - med).abs().median(); w.loc[w.fold == f_, 's'] = (w.loc[w.fold == f_, 's'] - med) / mad
  best = w.sort_values(['s', 'is_decoy'], ascending=[False, False]).drop_duplicates('pair_key').sort_values('s', ascending=False)
  T = np.cumsum(1 - best.is_decoy.values); D = np.cumsum(best.is_decoy.values)
  best['q'] = np.minimum.accumulate(np.minimum(1, (D + 1) / np.maximum(T, 1))[::-1])[::-1]
  return w, best


def n_at(t, f):
  e = t.sort_values('s', ascending=False).e.values; x = np.cumsum(e) * (1 + 1 / RR) / np.arange(1, len(e) + 1); ok = np.where(x <= f)[0]
  return int(ok.max() + 1) if len(ok) else 0


pops = {r: pd.read_csv(f'{O}/{r}/population.tsv', sep='\t') for r in ('raw01', 'raw02')}
sc = {r: np.load(f'{O}/{r}/scores_D11FWT.npy') for r in pops}
for a, b in (('raw01', 'raw02'), ('raw02', 'raw01')):
  _, bb = tdc(pops[b], sc[b])
  for qsel in (None, 0.5, 0.2, 0.05):
    keep = None if qsel is None else set(bb[bb.q <= qsel].pair_key)
    w, best = tdc(pops[a], sc[a], keep)
    t = w[w.is_decoy == 0].assign(e=lambda d: d.precursor_id.map(ent).values)
    tb = best[(best.is_decoy == 0) & (best.q <= 0.01)].assign(e=lambda d: d.precursor_id.map(ent).values)
    fq = tb.e.sum() * (1 + 1 / RR) / max(len(tb), 1)
    print(f'{a} selection by {b} pair q <= {qsel}: pairs {best.pair_key.nunique():6d} | N@FDP_c 1/2 % {n_at(t, 0.01)} / {n_at(t, 0.02)} | '
          f'q<=0.01 N {len(tb)} FDP_c {fq:.4f}')
