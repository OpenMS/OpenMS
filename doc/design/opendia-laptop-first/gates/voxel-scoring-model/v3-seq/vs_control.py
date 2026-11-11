"""PLAN_V3_SEQ rule 3: paired bootstrap of the arm's N@1 % / N@2 % against its control without the sequence input
(DS11 against D11, DVS against DV). Same metric and resampling as evaluate_v2.py (1000 resamples, seed 5)."""
import sqlite3, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; O = B + '/astral-voxel-pilot-20261110/v2'
R = 148983 / 151239
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c).set_index('precursor_id').e

def winners(pop, s):
    w = pop.assign(s=s).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').copy()
    for f_, g in w[w.is_decoy == 1].groupby('fold'):
        med = g.s.median(); mad = 1.4826 * (g.s - med).abs().median(); w.loc[w.fold == f_, 's'] = (w.loc[w.fold == f_, 's'] - med) / mad
    t = w[w.is_decoy == 0]
    return t.assign(e=t.precursor_id.map(ent).values).set_index('precursor_id')[['s', 'e']]

def n_at(s, e, t):
    o = np.argsort(-s, kind='stable'); e = e[o]; f = np.cumsum(e) * (1 + 1 / R) / np.arange(1, len(e) + 1); ok = np.where(f <= t)[0]
    return int(ok.max() + 1) if len(ok) else 0

rng = np.random.default_rng(5)
for raw in ('raw01', 'raw02'):
    pop = pd.read_csv(f'{O}/{raw}/population.tsv', sep='\t')
    for arm, ctl in (('DS11', 'D11'), ('DVS', 'DV')):
        try: a = winners(pop, np.load(f'{O}/{raw}/scores_{arm}.npy')); b = winners(pop, np.load(f'{O}/{raw}/scores_{ctl}.npy'))
        except FileNotFoundError: print(raw, arm, 'missing'); continue
        m = a.join(b.s.rename('b'), how='inner'); x, y, e = m.s.values, m.b.values, m.e.values.astype(int); d1, d2 = [], []
        for _ in range(1000):
            i = rng.integers(0, len(m), len(m)); d1.append(n_at(x[i], e[i], 0.01) - n_at(y[i], e[i], 0.01)); d2.append(n_at(x[i], e[i], 0.02) - n_at(y[i], e[i], 0.02))
        print(f'{raw} {arm} vs {ctl}: N@1% {n_at(x, e, 0.01)} vs {n_at(y, e, 0.01)} (90% {np.percentile(d1, 5):+.0f}..{np.percentile(d1, 95):+.0f}) | '
              f'N@2% {n_at(x, e, 0.02)} vs {n_at(y, e, 0.02)} (90% {np.percentile(d2, 5):+.0f}..{np.percentile(d2, 95):+.0f})', flush=True)
