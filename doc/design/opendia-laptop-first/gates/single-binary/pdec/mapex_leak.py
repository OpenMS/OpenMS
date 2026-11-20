# Diagnostic only (entrapment labels for evaluation): do the m_* features separate decoys from entrapment targets?
import sqlite3, sys, glob, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
L = '/home/singjc/.jcode/scratch/dia2/ctest_cache/c/hela_entrap_codon.predicted_decoys.pqp'
c = sqlite3.connect(f'file:{L}?mode=ro', uri=True)
e = pd.read_sql("""select p.ID precursor_id, p.DECOY d, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID
  join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c).set_index('precursor_id')
f = pd.read_csv(glob.glob('mapex/*.features.tsv')[0], sep='\t')
pid = [x for x in f.columns if 'precursor' in x][0]
f['e'] = f[pid].map(e.e).values; f['d'] = f[pid].map(e.d).values
cols = [x for x in f.columns if x.startswith('m_') or x.startswith('x_im')] + ['score']
for name, sub in (('all windows', f), ('top 20 % score', f[f.score >= f.score.quantile(0.8)])):
    false = sub[(sub.d == 1) | (sub.e == 1)]
    y = (false.d == 1).astype(int)
    print(name, 'decoys', int(y.sum()), 'entrapment', int((1 - y).sum()))
    for col in cols:
        print(f'  {col:22s} AUC(decoy vs entrapment) {roc_auc_score(y, false[col]):.3f}')
