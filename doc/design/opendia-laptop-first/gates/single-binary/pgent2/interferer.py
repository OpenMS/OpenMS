"""Shared-ladder test (D s0 CUR, PGENT2 lib). For each winner (entrapment or decoy) with best q <= 0.05 band and a random
sample of all winners: max number of its 12 library fragments (+-15 ppm) shared with ONE accepted human precursor
(q<=0.01) whose observed apex is within 10 s and whose precursor m/z falls in the same 25 Th isolation window (diaPASEF
windows are ~25 Th). Compare entrapment vs decoy."""
import sqlite3, numpy as np, pandas as pd, bisect
c=sqlite3.connect('file:hela_ath2.pqp?mode=ro',uri=True)
P=pd.read_sql("select ID pid, DECOY d, PRECURSOR_MZ pmz from PRECURSOR",c).set_index('pid')
ent=pd.read_sql("""select pp.PRECURSOR_ID pid, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
 from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by 1""",c).set_index('pid').e
T=pd.read_sql("select tpm.PRECURSOR_ID pid, t.PRODUCT_MZ mz from TRANSITION t join TRANSITION_PRECURSOR_MAPPING tpm on tpm.TRANSITION_ID=t.ID",c)
frags=T.groupby('pid').mz.apply(lambda s: np.sort(s.values)).to_dict()
w=pd.read_csv('win_D_cur_s0.tsv.gz',sep='\t',usecols=['precursor_id','is_decoy','q_value','score','rt_apex']); w=w[w.score>-1e29]
w['e']=w.precursor_id.map(ent).fillna(0).values*(w.is_decoy==0); w['pmz']=w.precursor_id.map(P.pmz)
hum=w[(w.is_decoy==0)&(w.e==0)&(w.q_value<=0.01)].sort_values('rt_apex')
hrt=hum.rt_apex.values; hid=hum.precursor_id.values; hmz=hum.pmz.values
def shared(pid, rt, pmz):
    f=frags.get(pid); 
    if f is None: return 0
    lo=bisect.bisect_left(hrt,rt-10); hi=bisect.bisect_right(hrt,rt+10); best=0
    for j in range(lo,hi):
        if abs(hmz[j]-pmz)>12.5 or hid[j]==pid: continue
        g=frags.get(hid[j]); 
        if g is None: continue
        idx=np.searchsorted(g,f); k=0
        for a,i in zip(f,idx):
            for ii in (i-1,i):
                if 0<=ii<len(g) and abs(g[ii]-a)<=a*15e-6: k+=1; break
        best=max(best,k)
    return best
cut=w[(w.is_decoy==0)&(w.q_value<=0.05)].score.min()
for name,sel in (('high band (score >= q0.05 cut)', w.score>=cut), ('random winners', w.index.isin(w.sample(20000,random_state=0).index))):
    s=w[sel]
    for lab,g in (('entrapment',s[s.e==1]),('decoy',s[s.is_decoy==1])):
        g=g.sample(min(len(g),3000),random_state=1)
        k=np.array([shared(p,r,m) for p,r,m in zip(g.precursor_id,g.rt_apex,g.pmz)])
        print(f'{name:32s} {lab:10s} n {len(k):5d} | mean shared {k.mean():.2f} | >=3 shared {np.mean(k>=3):.3f} | >=4 {np.mean(k>=4):.3f}')
