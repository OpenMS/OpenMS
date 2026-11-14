"""Descriptive (evaluation labels): which precursor pairs share amino-acid composition and charge within the scored
population (raw01)? Class pairs human/entrapment/decoy. Explains what features_J removes relative to features_I."""
import sqlite3, numpy as np, pandas as pd
D = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/gap-20261114/armE_raw01'
pop = pd.read_csv(f'{D}/population.tsv', sep='\t', usecols=['precursor_id']).drop_duplicates()
c = sqlite3.connect(f'file:{D}/library.pqp?mode=ro', uri=True)
s = pd.read_sql("""select p.ID pid, p.CHARGE z, p.DECOY d, pe.UNMODIFIED_SEQUENCE seq from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m
  on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID""", c).drop_duplicates('pid')
ent = pd.read_sql("""select pp.PRECURSOR_ID pid, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c).set_index('pid').e
s = s[s.pid.isin(pop.precursor_id)].copy()
s['cls'] = np.where(s.d == 1, 'decoy', np.where(s.pid.map(ent) == 1, 'entrap', 'human'))
s['key'] = s.seq.map(lambda x: ''.join(sorted(x))) + '/' + s.z.astype(str)
g = s.groupby('key')
multi = g.filter(lambda x: len(x) > 1)
pairs = multi.merge(multi, on='key'); pairs = pairs[pairs.pid_x < pairs.pid_y]
print('precursors in population', len(s), '| in a same-composition group', len(multi))
print(pairs.apply(lambda r: '-'.join(sorted((r.cls_x, r.cls_y))), axis=1).value_counts().to_string())
print('share of each class with a same-composition partner:', (multi.cls.value_counts() / s.cls.value_counts()).round(3).to_dict())
