"""Pseudo-reverse decoys keep the target's library intensities but in which fragment order? Check for a pair: are the
decoy's 6 search fragments the 'same ordinals' as the target's (OpenSwathDecoyGenerator maps each fragment to the reversed
sequence's ion of the same type and ordinal) and does the decoy keep the target's intensity per slot? And for entrapments,
DIALibGen picks its top-6 by predicted intensity for its own sequence. So an entrapment's fragments are its most intense
predicted ions; a decoy's are the target's top ions mapped to reversed-sequence ions, which are NOT the decoy's own most
intense. If the entrapment's chosen ions are intrinsically more likely to be observed (e.g. y ions of high ordinal),
the entrapment gets more fragment hits when random signal is present."""
import sqlite3, numpy as np, pandas as pd
c=sqlite3.connect('file:hela_ath2.pqp?mode=ro',uri=True)
T=pd.read_sql("""select tpm.PRECURSOR_ID pid, t.TYPE ty, t.ORDINAL ord, t.LIBRARY_INTENSITY li, t.PRODUCT_MZ mz, p.DECOY d, p.PRECURSOR_MZ pmz, p.CHARGE z
 from TRANSITION t join TRANSITION_PRECURSOR_MAPPING tpm on tpm.TRANSITION_ID=t.ID join PRECURSOR p on p.ID=tpm.PRECURSOR_ID""",c)
ent=pd.read_sql("""select pp.PRECURSOR_ID pid, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
 from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by 1""",c).set_index('pid').e
T['e']=T.pid.map(ent).fillna(0).astype(int); T['grp']=np.where(T.d==1,'decoy',np.where(T.e==1,'entrap','human'))
T['above_prec']=(T.mz>T.pmz).astype(int)
T=T.sort_values(['pid','li'],ascending=[True,False]); T['rk']=T.groupby('pid').cumcount()
top=T[T.rk<6]
print(top.groupby('grp').agg(y=('ty',lambda s:(s=='y').mean()), ord_mean=('ord','mean'), above_prec=('above_prec','mean'), mz_med=('mz','median')).round(3))
