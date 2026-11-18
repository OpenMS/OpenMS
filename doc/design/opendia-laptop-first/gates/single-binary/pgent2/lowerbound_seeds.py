"""PGENT2 (homolog-free Arabidopsis, r=1.56): per seed, lower bound N_E/N and combined at q<=0.01; and N_E/N at q<=0.005, 0.02."""
import sqlite3, glob, numpy as np, pandas as pd
c=sqlite3.connect('file:hela_ath2.pqp?mode=ro',uri=True)
ent=pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
 from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by 1""",c).set_index('precursor_id').e
r=1.5602; k=1+1/r
for arm in ('cur','new'):
    for qq in (0.005,0.01,0.02):
        L=[];C=[]
        for f in sorted(glob.glob(f'runs/D_{arm}_s*/*.gz')):
            t=pd.read_csv(f,sep='\t'); t=t[(t.score>-1e29)&(t.is_decoy==0)&(t.q_value<=qq)]
            E=t.precursor_id.map(ent).fillna(0).sum(); N=len(t); L.append(E/N); C.append(E*k/N)
        print(f'{arm} q<={qq}: lower {np.mean(L):.4f} ({min(L):.4f}-{max(L):.4f})  combined {np.mean(C):.4f}')
