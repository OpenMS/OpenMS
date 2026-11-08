"""Entrapment-matched N for raw 02 (same method as ../matched_fdp.py on raw 01): OSW+Percolator and the OpenDIA tool's
raw 02 stage-2 re-score (astral-raw02-stage2-check-20261108, an optimistic diagnostic). It also reports OSW's own
m_score <= 0.01 export size and FDP_c, as the calibration reference for both runs."""
import sqlite3, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
R = 148983 / 151239
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c)
O = pd.read_csv('/home/singjc/Documents/datasets/proteobench-astral-entrapment/100K_sampled_fasta/opendia_results/OpenDIA.results.tsv', sep='\t', usecols=['run_name', 'FullPeptideName', 'Charge', 'm_score', 'd_score', 'ProteinName'])

def curve(name, e):
    e = np.asarray(e); n = np.arange(1, len(e) + 1); fdp = np.cumsum(e) * (1 + 1 / R) / n
    out = [int(np.where(fdp <= t)[0].max() + 1) if (fdp <= t).any() else 0 for t in (0.005, 0.01, 0.015, 0.02)]
    print(f'{name:52s} N at FDP_c <= 0.5/1/1.5/2%: {out}')

for run in ('01', '02'):
    o = O[O.run_name.str.endswith('_' + run)].copy(); o['mkey'] = o.FullPeptideName + '/' + o.Charge.astype(str)
    o['e'] = o.ProteinName.apply(lambda s: int(all('_p_target' in p for p in str(s).split(';'))))
    o = o.sort_values('d_score', ascending=False).drop_duplicates('mkey')
    curve(f'raw {run} OSW+Percolator', o.e.values)
    a = o[o.m_score <= 0.01]
    print(f'   raw {run} OSW m_score<=0.01 export: N {len(a)}, entrapment {int(a.e.sum())}, FDP_c {a.e.sum() * (1 + 1 / R) / len(a):.4f}')
cr = pd.read_csv(B + '/astral-raw02-stage2-check-20261108/crossfit/crossfit.tsv', sep='\t')
cols = [c_ for c_ in ('precursor_best', 'in_rt_window') if c_ in cr.columns]
w = cr[(cr.is_decoy == 0)]
for c_ in cols:
    w = w[w[c_] == 1]
w = w.merge(ent, on='precursor_id').sort_values('stage2_score', ascending=False)
curve('raw 02 OpenDIA tool stage-2 re-score (filters: ' + ','.join(cols) + ')', w.e.values)
