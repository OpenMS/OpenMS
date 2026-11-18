"""Wen et al. 2025 bounds at q <= 0.01, per run/arm (pooled over seeds) for the shuffled 1:1 entrapment library (v1):
lower bound N_E/N (can show failure), combined N_E(1+1/r)/N (upper bound), and the paired upper bound at peptide level
(r=1 pairs from the FASTA): FDP = (N_E + N_{E>=s>T} + 2 N_{E>T>=s}) / (N_T + N_E), peptide score = best precursor score,
peptide accepted if its best precursor has q <= 0.01 (s = lowest accepted score)."""
import sqlite3, glob, numpy as np, pandas as pd
c=sqlite3.connect('file:hela_entrap_dec.pqp?mode=ro',uri=True)
m=pd.read_sql("""select p.ID precursor_id, pe.UNMODIFIED_SEQUENCE seq, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
 from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING pm on pm.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=pm.PEPTIDE_ID
 join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID where p.DECOY=0 group by p.ID""",c).set_index('precursor_id')
pair={}; prev=None
for l in open('hela_entrap.fasta'):
    if l.startswith('>'): name=l[1:].split('|')[1]
    else:
        s=l.strip()
        if name.endswith('_p_target'): pair[s]=prev
        else: prev=s
r=1.1225; k=1+1/r
for run in ('D','H'):
  for arm in ('cur','new'):
    lb=[];cb=[];pb=[]; E_=N_=0
    for f in sorted(glob.glob(f'runs/runs/{run}_{arm}_s*/*.gz')):
        t=pd.read_csv(f,sep='\t'); t=t[(t.score>-1e29)&(t.is_decoy==0)].join(m,on='precursor_id')
        acc=t[t.q_value<=0.01]; N=len(acc); E=int(acc.e.sum()); E_+=E; N_+=N
        s=acc.score.min()
        pep=t.groupby('seq').agg(score=('score','max'),e=('e','max'))
        Tp=pep[pep.e==0].score; Ep=pep[pep.e==1]
        acc_p=pep[pep.score>=s]; NT=int((acc_p.e==0).sum()); NE=int((acc_p.e==1).sum())
        par=Ep.index.map(pair); ts=Tp.reindex(par).fillna(-np.inf).values; es=Ep.score.values
        a=int(((es>=s)&(ts<s)).sum()); b=int(((es>ts)&(ts>=s)).sum())
        pb.append((NE+a+2*b)/(NT+NE)); lb.append(E/N); cb.append(E*k/N)
    print(f'{run} {arm}: precursor lower bound {np.mean(lb):.4f}  combined {np.mean(cb):.4f} | peptide paired upper bound {np.mean(pb):.4f}  (mean over {len(lb)} seeds)')

print('--- peptide level, same acceptance set: lower / combined / paired')
for run in ('D','H'):
  for arm in ('cur','new'):
    L=[];C=[];Pp=[]
    for f in sorted(glob.glob(f'runs/runs/{run}_{arm}_s*/*.gz')):
        t=pd.read_csv(f,sep='\t'); t=t[(t.score>-1e29)&(t.is_decoy==0)].join(m,on='precursor_id')
        s=t[t.q_value<=0.01].score.min()
        pep=t.groupby('seq').agg(score=('score','max'),e=('e','max'))
        Tp=pep[pep.e==0].score; Ep=pep[pep.e==1]
        acc_p=pep[pep.score>=s]; NT=int((acc_p.e==0).sum()); NE=int((acc_p.e==1).sum())
        par=Ep.index.map(pair); ts=Tp.reindex(par).fillna(-np.inf).values; es=Ep.score.values
        a=int(((es>=s)&(ts<s)).sum()); b=int(((es>ts)&(ts>=s)).sum())
        L.append(NE/(NT+NE)); C.append(2*NE/(NT+NE)); Pp.append((NE+a+2*b)/(NT+NE))
    print(f'{run} {arm}: lower {np.mean(L):.4f}  combined(r=1) {np.mean(C):.4f}  paired {np.mean(Pp):.4f}')
