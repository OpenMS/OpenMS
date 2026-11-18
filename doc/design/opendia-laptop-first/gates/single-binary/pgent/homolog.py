"""b2: are accepted Arabidopsis entrapments near-homologs of human peptides? For each entrapment peptide: min Hamming
distance (I=L) to any human tryptic peptide of the same length (0-2 missed cleavages, 7-44), and the longest suffix
shared with a human peptide of the same length (y-ion overlap). Accepted (q <= 0.01) vs a random sample of all
entrapment target winners (null)."""
import sqlite3, glob, numpy as np, pandas as pd, collections, random
exec(open('make_fasta.py').read().split("c=sqlite3")[0])
HUM='/home/singjc/Documents/github/testing_auto_irt/diaPASEF_evosep_published/uniprotkb_Human_AND_reviewed_true_AND_m_2026_02_26.fasta'
norm=lambda p:p.replace('I','L')
hum=collections.defaultdict(set)
for s in read(HUM):
    for p in digest(s): hum[len(p)].add(norm(p))
A={ln:np.frombuffer(''.join(sorted(v)).encode(),dtype=np.uint8).reshape(-1,ln) for ln,v in hum.items()}
def near(p):
    q=np.frombuffer(norm(p).encode(),dtype=np.uint8); M=A.get(len(p))
    if M is None: return 99,0
    neq=(M!=q); ham=int(neq.sum(1).min())
    # longest shared suffix
    rev=neq[:,::-1]; first=np.where(rev.any(1),rev.argmax(1),len(p)); return ham,int(first.max())
c=sqlite3.connect('file:hela_ath.pqp?mode=ro',uri=True)
P=pd.read_sql("""select p.ID precursor_id, pe.UNMODIFIED_SEQUENCE seq, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
 from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
 join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID where p.DECOY=0 group by p.ID""",c).set_index('precursor_id')
for arm in ('cur','new'):
    t=pd.read_csv(glob.glob(f'runs/D_{arm}_s0/*.gz')[0],sep='\t')
    w=t[(t.score>-1e29)&(t.is_decoy==0)].sort_values('score',ascending=False).drop_duplicates('precursor_id').join(P,on='precursor_id')
    ew=w[w.e==1]
    acc=ew[ew.q_value<=0.01]; top=ew[ew.q_value<=0.005]; nul=ew.sample(min(1000,len(ew)),random_state=0)
    for name,g in (('accepted q<=0.01',acc),('  of which q<=0.005',top),('null: random entrapment winners',nul)):
        r=np.array([near(s) for s in g.seq.unique()])
        print(f'{arm} {name:32s} n {len(r):4d} | Hamming<=1 {np.mean(r[:,0]<=1):.2f} <=2 {np.mean(r[:,0]<=2):.2f} | shared y-suffix>=6 {np.mean(r[:,1]>=6):.2f}')
