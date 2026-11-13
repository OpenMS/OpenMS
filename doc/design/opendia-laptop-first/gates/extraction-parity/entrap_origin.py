"""Are entrapment peptides composition-preserving permutations of human peptides (same precursor mass, shared
termini)? And are decoys mass-shifted from their target? Library only, no run data."""
import sqlite3, pandas as pd, numpy as np
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
pep = pd.read_sql("""select pe.ID, pe.UNMODIFIED_SEQUENCE seq, pe.DECOY decoy,
  min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e from PEPTIDE pe
  join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pe.ID""", c)
t = pep[pep.decoy == 0].drop_duplicates('seq').copy()
t['comp'] = t.seq.map(lambda s: ''.join(sorted(s)))
hum = t[t.e == 0]; ent = t[t.e == 1]
hc = set(hum.comp); hs = set(hum.seq)
print('human peptides', len(hum), 'entrapment peptides', len(ent))
print('entrapment with identical composition to a human peptide:', round(ent.comp.isin(hc).mean(), 4))
# Same composition, same C-terminal residue, same first residue
hkey = set(zip(hum.comp, hum.seq.str[-1], hum.seq.str[0]))
print('  ... and same first and last residue:', round(np.mean([(a, b[-1], b[0]) in hkey for a, b in zip(ent.comp, ent.seq)]), 4))
# Random baseline: human peptides sharing composition with another human peptide
cnt = hum.comp.value_counts()
print('human peptides whose composition is shared by another human peptide:', round((hum.comp.map(cnt) > 1).mean(), 4))
# Example pairs
m = ent[ent.comp.isin(hc)].head(5).merge(hum[['comp', 'seq']], on='comp', suffixes=('_entrap', '_human'))
print(m[['seq_entrap', 'seq_human']].to_string(index=False))
# Precursor m/z offset between each target precursor and its decoy
pre = pd.read_sql('select ID, TRAML_ID traml, PRECURSOR_MZ mz, CHARGE z, DECOY decoy, LIBRARY_RT rt from PRECURSOR', c)
tg = pre[pre.decoy == 0].set_index('traml'); dc = pre[pre.decoy == 1].assign(traml=lambda d: d.traml.str.replace('^DECOY_', '', regex=True)).set_index('traml')
j = tg.join(dc, rsuffix='_d', how='inner')
dm = (j.mz_d - j.mz) * j.z
print('decoy - target precursor mass (Da), quantiles 5/25/50/75/95:', np.round(np.quantile(dm, [0.05, 0.25, 0.5, 0.75, 0.95]), 3))
print('share of decoys with |mass shift| < 0.01 Da:', round(float((dm.abs() < 0.01).mean()), 4), ' RT identical:', round(float((j.rt == j.rt_d).mean()), 4))
