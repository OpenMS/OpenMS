"""Teacher-free RT fine-tune data: train = our own accepted IDs (tool stage-2, q<=0.01) with observed apex RT, split 85/15
for early stopping; eval = DIA-NN-confirmed precursors not in train (DIA-NN RT only used to score the predictor)."""
import sqlite3, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; D = B + '/astral-design-raw01-20261107'
cr = pd.read_csv(D + '/stage2-boundary-design/boundary/crossfit.tsv', sep='\t', usecols=['sample_id', 'precursor_id', 'is_decoy', 'precursor_best', 'q_value'])
s = pd.read_csv(D + '/pipeline/search.features.tsv', sep='\t', usecols=['candidate_row', 'rt_apex']).rename(columns={'candidate_row': 'sample_id'})
w = cr[cr.precursor_best == 1].merge(s, on='sample_id'); w['q'] = pd.to_numeric(w.q_value, errors='coerce')
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
pep = pd.read_sql('select pr.ID precursor_id, pr.CHARGE charge, pr.LIBRARY_RT lib_rt, p.MODIFIED_SEQUENCE modified_sequence from PRECURSOR pr join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=pr.ID join PEPTIDE p on p.ID=m.PEPTIDE_ID', c)
tr = w[(w.is_decoy == 0) & (w.q <= 0.01)][['precursor_id', 'rt_apex']].merge(pep, on='precursor_id')
tr = tr.drop_duplicates('modified_sequence')
rng = np.random.default_rng(0); val = rng.random(len(tr)) < 0.15
ref = pd.read_parquet(B + '/astral-diann-gap-20261108/ref_exact.parquet'); ref = ref[ref.Run.str.endswith('_01')]
ev = ref[['precursor_id']].assign(diann_rt=ref.RT.values * 60).merge(pep, on='precursor_id')
ev = ev[~ev.modified_sequence.isin(tr.modified_sequence)].drop_duplicates('modified_sequence')
cols = lambda d, rt: d.assign(retention_time=rt)[['modified_sequence', 'charge', 'retention_time']]
cols(tr[~val], tr.rt_apex[~val] / 60).to_csv('train.csv', index=False)
cols(tr[val], tr.rt_apex[val] / 60).to_csv('val.csv', index=False)
cols(ev, ev.diann_rt / 60).to_csv('eval.csv', index=False)
ev.to_csv('eval_meta.csv', index=False)
print('train', (~val).sum(), 'val', val.sum(), 'eval', len(ev))
