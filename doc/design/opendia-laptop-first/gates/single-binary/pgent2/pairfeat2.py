"""Decompose the high-score asymmetry. Among E-vs-own-decoy pairs (best q<=0.05), the winner is decided by co-elution
(f12_coelution AUC 0.97) - expected, since co-elution IS signal. The question is whether a false target gets co-eluting
signal more often than its decoy. Compare the distribution of the strongest evidence features over ALL entrapment vs ALL
decoy winners (not only high-scoring): fraction with f12_support >= 3 and f12_coelution > 0.5."""
import sqlite3, numpy as np, pandas as pd
c=sqlite3.connect('file:hela_ath2.pqp?mode=ro',uri=True)
ent=pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
 from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by 1""",c).set_index('precursor_id').e
w=pd.read_csv('win_D_cur_s0.tsv.gz',sep='\t'); w=w[w.score>-1e29]
w['e']=w.precursor_id.map(ent).fillna(0).values*(w.is_decoy==0)
E=w[w.e==1]; D=w[w.is_decoy==1]
for thr in (2,3,4,5,6):
    fe=(E.f12_support>=thr).mean(); fd=(D.f12_support>=thr).mean()
    print(f'f12_support >= {thr}: entrapment {fe:.5f} decoy {fd:.5f} ratio {fe/max(fd,1e-9):.2f}')
for thr in (0.5,0.7,0.8,0.9):
    fe=((E.f12_support>=3)&(E.f12_coelution>thr)).mean(); fd=((D.f12_support>=3)&(D.f12_coelution>thr)).mean()
    print(f'support>=3 & coelution>{thr}: entrapment {fe:.5f} decoy {fd:.5f} ratio {fe/max(fd,1e-9):.2f}')
