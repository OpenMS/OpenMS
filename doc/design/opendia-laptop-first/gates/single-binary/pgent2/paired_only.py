"""Predicted-decoy library: 32 % of targets have no decoy partner (no competition, so they win by default). Recompute TDC
q and the bounds restricted to pairs that have both members, for the predicted-decoy and the pseudo-reverse (PGENT2) libs."""
import sqlite3, glob, numpy as np, pandas as pd
def run(lib, files, label):
    c=sqlite3.connect(f'file:{lib}?mode=ro',uri=True)
    P=pd.read_sql("""select p.ID precursor_id, p.TRAML_ID t, p.DECOY d, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
     from PRECURSOR p left join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID left join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID
     left join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""",c)
    P['key']=np.where(P.d==1,P.t.str.replace('^DECOY_','',regex=True),P.t)
    paired=set(P[P.d==1].key)&set(P[P.d==0].key)
    P['paired']=P.key.isin(paired)
    tg=P[(P.d==0)&P.paired]; r=(tg.e==1).sum()/(tg.e==0).sum(); k=1+1/r
    for f in files:
        t=pd.read_csv(f,sep='\t'); t=t[t.score>-1e29].sort_values('score',ascending=False).drop_duplicates('precursor_id').merge(P,on='precursor_id')
        un=t[(t.d==0)&~t.paired]
        a=t[t.paired]
        best=a.sort_values(['score','d'],ascending=[False,False]).drop_duplicates('key')
        D=np.cumsum(best.d.values); T=np.cumsum(1-best.d.values); q=np.minimum.accumulate(np.minimum(1,(D+1)/np.maximum(T,1))[::-1])[::-1]
        acc=best[(q<=0.01)&(best.d==0)]; N=len(acc); E=int(acc.e.sum())
        accu=t[(t.d==0)&(t.q_value<=0.01)]; Nu=int((~accu.paired).sum()); Eu=int(accu[~accu.paired].e.sum())
        print(f'{label} {f.split("/")[-2]}: unpaired targets in library {(~P[P.d==0].paired).sum()}, accepted unpaired {Nu} (entrap {Eu}) | '
              f'paired-only TDC q<=0.01: N {N} lower {E/N:.4f} combined {E*k/N:.4f} (r {r:.2f})')
run('hela_pdec.pqp', sorted(glob.glob('runs/D_*_s0/*.gz')), 'pdec')
run('../entlib3/hela_ath2.pqp', sorted(glob.glob('../entlib3/runs/D_*_s0/*.gz')), 'pseudo-rev')
