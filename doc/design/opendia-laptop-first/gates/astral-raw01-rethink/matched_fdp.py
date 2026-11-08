"""Honest comparison at matched entrapment FDP (raw 01): OSW+Percolator (same library) vs our tool stage-2
(search+boundary) and the offline MLP ensemble. Entrapment = all protein accessions contain _p_target; FDP_c combined."""
import sqlite3, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; D = B + '/astral-design-raw01-20261107'
R = 148983 / 151239
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c)
lib = pd.read_parquet(B + '/astral-diann-gap-20261108/lib_targets.parquet'); lib['mkey'] = lib.modseq + '/' + lib.charge.astype(str)
o = pd.read_csv('/home/singjc/Documents/datasets/proteobench-astral-entrapment/100K_sampled_fasta/opendia_results/OpenDIA.results.tsv', sep='\t', usecols=['run_name', 'FullPeptideName', 'Charge', 'm_score', 'd_score', 'ProteinName'])
o = o[o.run_name.str.endswith('_01')]; o['mkey'] = o.FullPeptideName + '/' + o.Charge.astype(str)
o['e'] = o.ProteinName.apply(lambda s: int(all('_p_target' in p for p in str(s).split(';'))))
o = o.sort_values('d_score', ascending=False).drop_duplicates('mkey')
def curve(name, score_sorted_e):
    e = np.asarray(score_sorted_e); n = np.arange(1, len(e) + 1); ne = np.cumsum(e); fdp = ne * (1 + 1 / R) / n
    out = []
    for t in (0.005, 0.01, 0.015, 0.02, 0.033):
        ok = np.where(fdp <= t)[0]; out.append(int(ok.max() + 1) if len(ok) else 0)
    print(f'{name:38s} N at FDP_c <= 0.5% / 1% / 1.5% / 2% / 3.3%: {out}')
curve('OSW+Percolator (m_score<=0.01 export)', o.e.values)
cr = pd.read_csv(D + '/stage2-boundary-design/boundary/crossfit.tsv', sep='\t', usecols=['precursor_id', 'is_decoy', 'precursor_best', 'stage2_score', 'in_rt_window'])
w = cr[(cr.precursor_best == 1) & (cr.is_decoy == 0) & (cr.in_rt_window == 1)].merge(ent, on='precursor_id').sort_values('stage2_score', ascending=False)
curve('ours: tool stage-2 search+boundary', w.e.values)
print('OSW exported', len(o), 'entrap', int(o.e.sum()))
