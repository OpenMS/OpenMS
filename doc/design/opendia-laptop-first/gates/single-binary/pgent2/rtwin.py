"""Do false targets get more co-eluting evidence because of WHERE they search? A decoy searches at its target's predicted
RT (copied); an entrapment at its own predicted RT. Check the library RT distribution: entrapment vs human vs decoy, and
the density of real HeLa signal there (human targets accepted at q<=0.01, observed apex)."""
import sqlite3, numpy as np, pandas as pd
c=sqlite3.connect('file:hela_ath2.pqp?mode=ro',uri=True)
P=pd.read_sql("select ID pid, DECOY d, LIBRARY_RT lrt, LIBRARY_DRIFT_TIME im, PRECURSOR_MZ mz, CHARGE z from PRECURSOR",c)
ent=pd.read_sql("""select pp.PRECURSOR_ID pid, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
 from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by 1""",c).set_index('pid').e
P['e']=P.pid.map(ent).fillna(0).astype(int); P['grp']=np.where(P.d==1,'decoy',np.where(P.e==1,'entrap','human'))
w=pd.read_csv('win_D_cur_s0.tsv.gz',sep='\t',usecols=['precursor_id','is_decoy','q_value','rt_apex'])
acc=w[(w.is_decoy==0)&(w.q_value<=0.01)].merge(P,left_on='precursor_id',right_on='pid')
acc=acc[acc.grp=='human']
# density of accepted human peptides per library-RT decile and per (library RT, IM, m/z) cell
edges=np.quantile(P.lrt,np.linspace(0,1,11))
P['rtbin']=np.clip(np.searchsorted(edges,P.lrt)-1,0,9)
dens=pd.cut(acc.lrt,edges,labels=False,include_lowest=True).value_counts(normalize=True).sort_index()
print('share of each group per library-RT decile (and accepted human density):')
print(pd.concat([P.groupby('grp').rtbin.value_counts(normalize=True).unstack(0).round(3), dens.rename('accepted_human').round(3)],axis=1))
for g in ('human','entrap','decoy'):
    s=P[P.grp==g]; print(g,'IM median',round(s.im.median(),3),'m/z median',round(s.mz.median(),1),'z2 share',round((s.z==2).mean(),3))
