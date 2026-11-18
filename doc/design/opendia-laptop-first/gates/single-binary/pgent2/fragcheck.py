"""Why do false targets get more co-eluting fragment hits than their pseudo-reverse decoys? Compare the library fragment
m/z sets: (1) fragment m/z distribution of entrapments vs decoys; (2) how many of a precursor's 12 fragments fall on
common human fragment m/z (a pool built from the HeLa target fragments, binned at 10 ppm). A false target from a real
proteome shares sub-sequences (dipeptides, b2/y2 etc.) with real peptides more often than a reversed sequence?"""
import sqlite3, numpy as np, pandas as pd
c=sqlite3.connect('file:hela_ath2.pqp?mode=ro',uri=True)
T=pd.read_sql("""select tpm.PRECURSOR_ID pid, t.PRODUCT_MZ mz, t.TYPE ty, t.ORDINAL ord, t.CHARGE fz, p.DECOY d
 from TRANSITION t join TRANSITION_PRECURSOR_MAPPING tpm on tpm.TRANSITION_ID=t.ID join PRECURSOR p on p.ID=tpm.PRECURSOR_ID""",c)
ent=pd.read_sql("""select pp.PRECURSOR_ID pid, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
 from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by 1""",c).set_index('pid').e
T['e']=T.pid.map(ent).fillna(0).astype(int)
T['grp']=np.where(T.d==1,'decoy',np.where(T.e==1,'entrap','human'))
print(T.groupby('grp').agg(n=('mz','size'),mz_med=('mz','median'),ord_med=('ord','median'),y_share=('ty',lambda s:(s=='y').mean()),z2=('fz',lambda s:(s==2).mean())).round(3))
hum=T[T.grp=='human'].mz.values
bins=np.unique(np.round(np.log(hum)/1e-5).astype(np.int64))
def onhuman(m):
    k=np.round(np.log(m)/1e-5).astype(np.int64)
    return np.isin(k,bins)|np.isin(k-1,bins)|np.isin(k+1,bins)
for g in ('entrap','decoy'):
    s=T[T.grp==g]; s=s.assign(h=onhuman(s.mz.values))
    per=s.groupby('pid').h.sum()
    print(g, 'fragment on a human-target fragment m/z (+-10 ppm):', round(s.h.mean(),3), '| precursors with >=3 such fragments:', round((per>=3).mean(),3))
