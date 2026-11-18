"""Paired TDC assumption check (D, PGENT2 library, seeds 0-4): for entrapment target precursors and their own paired decoy,
P(entrapment beats its decoy) among pairs whose better member reaches a given run q. Must be 0.5 if decoys model false
targets. Also the same for human-target pairs (expected > 0.5, true targets present)."""
import sqlite3, glob, numpy as np, pandas as pd
c=sqlite3.connect('file:hela_ath2.pqp?mode=ro',uri=True)
P=pd.read_sql("""select p.ID precursor_id, p.TRAML_ID t, p.DECOY d, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
 from PRECURSOR p left join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID left join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID
 left join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""",c)
P['key']=np.where(P.d==1,P.t.str.replace('^DECOY_','',regex=True),P.t)
tg=P[P.d==0].set_index('key'); dc=P[P.d==1].set_index('key')
pairs=tg[['precursor_id','e']].join(dc[['precursor_id']],rsuffix='_dec',how='inner')
for arm in ('cur','new'):
    rows=[]
    for f in sorted(glob.glob(f'runs/D_{arm}_s*/*.gz')):
        t=pd.read_csv(f,sep='\t'); t=t[t.score>-1e29].sort_values('score',ascending=False).drop_duplicates('precursor_id').set_index('precursor_id')
        x=pairs.copy(); x['st']=t.score.reindex(x.precursor_id).values; x['sd']=t.score.reindex(x.precursor_id_dec).values
        x['qt']=t.q_value.reindex(x.precursor_id).values; x['qd']=t.q_value.reindex(x.precursor_id_dec).values
        x=x.dropna(subset=['st','sd']); x['best_q']=np.minimum(x.qt,x.qd); x['win']=(x.st>x.sd).astype(int)
        rows.append(x)
    x=pd.concat(rows)
    for qq in (0.01,0.05,0.2,1.0):
        b=x[x.best_q<=qq]
        for lab,g in (('entrapment',b[b.e==1]),('human',b[b.e==0])):
            n=len(g); k=g.win.sum()
            se=np.sqrt(0.25/n) if n else 0
            print(f'{arm} best q<={qq:<4}: {lab:10s} pairs {n:7d}  P(target beats its decoy) {k/max(n,1):.3f} (+-{1.64*se:.3f})')
