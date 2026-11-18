"""Sequence-statistics hypothesis: natural sequences (human targets, Arabidopsis entrapments) share local residue
statistics near the C-terminus; pseudo-reverse decoys put the target's N-terminal residue in front of K/R. Real co-eluting
signal comes from natural (human) peptides, so low-ordinal y ions of a natural false target hit real fragment m/z more
often. Check: residue at position n-1 (before C-term K/R) distributions, and the share of y1-y3 among top-6 fragments."""
import sqlite3, numpy as np, pandas as pd
c=sqlite3.connect('file:hela_ath2.pqp?mode=ro',uri=True)
pe=pd.read_sql("""select distinct pe.UNMODIFIED_SEQUENCE s, pe.DECOY d, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
 from PEPTIDE pe left join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID left join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by pe.ID""",c)
pe['grp']=np.where(pe.d==1,'decoy',np.where(pe.e==1,'entrap','human'))
pe['pre']=pe.s.str[-2]; pe['pre2']=pe.s.str[-3:-1]
f=pe.groupby('grp').pre.value_counts(normalize=True).unstack(0).fillna(0)
tv=lambda a,b: 0.5*np.abs(f[a]-f[b]).sum()
print('residue before C-term: total variation distance human-entrap', round(tv('human','entrap'),3), 'human-decoy', round(tv('human','decoy'),3))
g=pe.groupby('grp').pre2.value_counts(normalize=True).unstack(0).fillna(0)
tv2=lambda a,b: 0.5*np.abs(g[a]-g[b]).sum()
print('dipeptide before C-term: TV human-entrap', round(tv2('human','entrap'),3), 'human-decoy', round(tv2('human','decoy'),3))
T=pd.read_sql("""select tpm.PRECURSOR_ID pid, t.TYPE ty, t.ORDINAL ord, t.LIBRARY_INTENSITY li, p.DECOY d
 from TRANSITION t join TRANSITION_PRECURSOR_MAPPING tpm on tpm.TRANSITION_ID=t.ID join PRECURSOR p on p.ID=tpm.PRECURSOR_ID""",c)
ent=pd.read_sql("""select pp.PRECURSOR_ID pid, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
 from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by 1""",c).set_index('pid').e
T['grp']=np.where(T.d==1,'decoy',np.where(T.pid.map(ent).fillna(0)==1,'entrap','human'))
T=T.sort_values(['pid','li'],ascending=[True,False]); T['rk']=T.groupby('pid').cumcount(); top=T[T.rk<6]
print(top.assign(lowy=(top.ty=='y')&(top.ord<=3)).groupby('grp').lowy.mean().round(3).rename('share of top-6 that are y1-y3'))
