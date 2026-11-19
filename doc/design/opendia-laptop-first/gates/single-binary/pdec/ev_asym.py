"""PRECHECK_EVIDENCE asymmetry check (as PRECHECK_PGQ): for each feature, AUC of entrapment-target best windows against
decoy best windows among precursors whose pair has a best q <= 0.05. A fair feature has AUC ~0.5 (entrapments are false
targets, like decoys). Also AUC of HeLa targets vs decoys (separation) for the same rows. Labels for evaluation only."""
import sqlite3, glob, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
W = '/home/singjc/.jcode/scratch/dia2/gap7/ev/D_new_s0'
L = '/home/singjc/.jcode/scratch/dia2/ctest_cache/c/hela_entrap_codon.predicted_decoys.pqp'
F = pd.read_csv(glob.glob(f'{W}/*.features.tsv')[0], sep='\t')
c = sqlite3.connect(f'file:{L}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID pid, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c).set_index('pid').e
B = F[F.score > -1e29].copy()
B['e'] = pq.reindex(B.precursor_id).values
B = B[B.q_value <= 0.05]
meta = ['candidate_row', 'precursor_id', 'is_decoy', 'fold', 'rt_apex', 'score', 'q_value', 'e']
X = [x for x in B.columns if x not in meta]
ent = B[(B.is_decoy == 0) & (B.e == 1)]; dec = B[B.is_decoy == 1]; hel = B[(B.is_decoy == 0) & (B.e == 0)]
print(f'kept windows q<=0.05: entrapment {len(ent)}, decoy {len(dec)}, HeLa {len(hel)}')
rows = []
for x in X + ['score']:
    a = roc_auc_score(np.r_[np.ones(len(ent)), np.zeros(len(dec))], np.r_[ent[x].fillna(0), dec[x].fillna(0)])
    h = roc_auc_score(np.r_[np.ones(len(hel)), np.zeros(len(dec))], np.r_[hel[x].fillna(0), dec[x].fillna(0)])
    rows.append(dict(feature=x, ent_vs_dec=round(a, 3), hela_vs_dec=round(h, 3), flag='*' if abs(a - 0.5) > 0.05 else ''))
print(pd.DataFrame(rows).to_string(index=False))
