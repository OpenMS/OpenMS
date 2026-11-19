"""Diagnostic (labels for evaluation only): do accepted entrapment winners share fragments with a co-eluting identified HeLa
target more than accepted decoy winners do? H and D, CUR, seed 0. For each accepted winner (q<=0.01): among identified HeLa
targets (q<=0.01, non-entrapment) with |rt_apex diff| <= 10 s and |precursor m/z diff| <= 12.5, the maximum number of the
winner's fragments within 20 ppm of that neighbour's fragments."""
import sqlite3, glob, numpy as np, pandas as pd, bisect
c=sqlite3.connect('file:cut/hela_r5.predicted_decoys.pqp?mode=ro',uri=True)
pq=pd.read_sql("""select p.ID precursor_id, p.TRAML_ID t, p.DECOY d, p.PRECURSOR_MZ mz, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
 from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID
 join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""",c).set_index('precursor_id')
pq['key']=np.where(pq.d==1,pq.t.str.replace('^DECOY_','',regex=True),pq.t)
te=pq[pq.d==0].set_index('key').e; pq['pe']=pq.key.map(te)
def run(f):
    t=pd.read_csv(f,sep='\t'); t=t[t.score>-1e29].sort_values('score',ascending=False).drop_duplicates('precursor_id')
    t=t[t.q_value<=0.01].join(pq[['d','mz','e','pe']],on='precursor_id')
    ids=t.precursor_id.tolist()
    fr=pd.read_sql(f"select tm.PRECURSOR_ID p, tr.PRODUCT_MZ m from TRANSITION_PRECURSOR_MAPPING tm join TRANSITION tr on tr.ID=tm.TRANSITION_ID where tm.PRECURSOR_ID in ({','.join(map(str,ids))})",c)
    frag=fr.groupby('p').m.apply(lambda s: np.sort(s.values)).to_dict()
    hela=t[(t.d==0)&(t.e==0)].sort_values('rt_apex'); hr=hela.rt_apex.values; hm=hela.mz.values; hp=hela.precursor_id.values
    def shared(row):
        a=frag.get(row.precursor_id); 
        if a is None: return 0
        lo=bisect.bisect_left(hr,row.rt_apex-10); hi=bisect.bisect_right(hr,row.rt_apex+10); best=0
        for k in range(lo,hi):
            if hp[k]==row.precursor_id or abs(hm[k]-row.mz)>12.5: continue
            b=frag.get(hp[k]); 
            if b is None: continue
            idx=np.searchsorted(b,a); n=0
            for v,i in zip(a,idx):
                for j in (i-1,i):
                    if 0<=j<len(b) and abs(b[j]-v)<=v*20e-6: n+=1; break
            best=max(best,n)
        return best
    out={}
    for lab,sub in (('entrapment target',t[(t.d==0)&(t.e==1)]),('decoy, entrapment pair',t[(t.d==1)&(t.pe==1)]),('decoy, HeLa pair',t[(t.d==1)&(t.pe==0)]),('HeLa target (sample 1000)',t[(t.d==0)&(t.e==0)].sample(1000,random_state=0))):
        s=sub.apply(shared,axis=1); out[lab]=(len(s), s.mean(), (s>=2).mean(), (s>=3).mean())
    return out
for f in sorted(glob.glob('runs/[DH]_cur_s0/*.gz')):
    print(f.split('/')[1])
    for k,(n,m,g2,g3) in run(f).items(): print(f'  {k:28s} n {n:5d}  mean shared {m:.2f}  >=2 {g2:.3f}  >=3 {g3:.3f}')
