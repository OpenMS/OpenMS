"""PRECHECK_TDC (diagnostic on ODP D and H outputs): at the tool's own cut q<=c, count accepted target winners T,
accepted decoy winners D, entrapments E among T, and accepted decoys whose pair target is an entrapment (D_E) or HeLa (D_H).
Paired TDC is unbiased for entrapments iff E ~ D_E. Also recompute alternative q estimators on the same ranking."""
import sqlite3, glob, sys, numpy as np, pandas as pd
c=sqlite3.connect('file:cut/hela_r5.predicted_decoys.pqp?mode=ro',uri=True)
pq=pd.read_sql("""select p.ID precursor_id, p.TRAML_ID t, p.DECOY d, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
 from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID
 join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""",c)
pq['key']=np.where(pq.d==1,pq.t.str.replace('^DECOY_','',regex=True),pq.t)
tgt_e=pq[pq.d==0].set_index('key').e
pq['pe']=pq.key.map(tgt_e)          # entrapment status of the pair (via its target)
info=pq.set_index('precursor_id')[['d','pe']]
for R in ('D','H'):
  for arm in ('cur','new'):
    agg={c_:dict(T=0,D=0,E=0,DE=0,DH=0) for c_ in (0.005,0.01,0.02)}
    for f in sorted(glob.glob(f'runs/{R}_{arm}_s*/*.gz')):
        t=pd.read_csv(f,sep='\t'); t=t[t.score>-1e29].sort_values('score',ascending=False).drop_duplicates('precursor_id')
        t=t.join(info,on='precursor_id')
        for cut,a in agg.items():
            x=t[t.q_value<=cut]; tt=x[x.is_decoy==0]; dd=x[x.is_decoy==1]
            a['T']+=len(tt); a['D']+=len(dd); a['E']+=int(tt.pe.sum()); a['DE']+=int(dd.pe.sum()); a['DH']+=int((dd.pe==0).sum())
    for cut,a in agg.items():
        print(f"{R} {arm} q<={cut}: T {a['T']} D {a['D']} (D/T {a['D']/a['T']:.4f})  E {a['E']} (E/T {a['E']/a['T']:.4f})  D_E {a['DE']} D_H {a['DH']}  E/D_E {a['E']/max(a['DE'],1):.3f}")
