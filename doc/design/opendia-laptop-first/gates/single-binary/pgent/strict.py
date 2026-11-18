"""b1 strict entrapment set: drop entrapment precursors whose peptide is within Hamming 2 (I=L) of a same-length human
tryptic peptide, from the accepted set AND from r (so the estimator stays consistent). FDP_c at q<=0.01 and rho."""
import sqlite3, glob, numpy as np, pandas as pd
exec(open('homolog.py').read().split("for arm in")[0])
seqs=P[P.e==1].seq.unique()
far={s for s in seqs if near(s)[0]>2}
P['strict']=np.where(P.e==1, P.seq.isin(far).astype(int), 0)
P['drop']=(P.e==1)&(P.strict==0)
r=(P.strict==1).sum()/((P.e==0).sum()); k=1+1/r
print(f'entrapment peptides {len(seqs)}, kept as strict (Hamming > 2) {len(far)} ({len(far)/len(seqs):.2f}); r_strict {r:.4f}')
for arm in ('cur','new'):
    t=pd.read_csv(glob.glob(f'runs/D_{arm}_s0/*.gz')[0],sep='\t')
    w=t[t.score>-1e29].sort_values('score',ascending=False).drop_duplicates('precursor_id').join(P[['strict','drop']],on='precursor_id')
    w=w[~w['drop'].fillna(False).astype(bool)]
    acc=w[(w.is_decoy==0)&(w.q_value<=0.01)]; N=len(acc); E=int(acc.strict.fillna(0).sum())
    tt=w[w.is_decoy==0].sort_values(['q_value','score'],ascending=[True,False]); e=tt.strict.fillna(0).values
    f=np.cumsum(e)*k/np.arange(1,len(e)+1); ok=np.where(f<=0.01)[0]; n1=int(ok.max()+1) if len(ok) else 0
    print(f'{arm}: q<=0.01 N {N}, strict entrapments {E}, FDP_c {E*k/N:.4f} (90 % {__import__("scipy.stats").stats.chi2.ppf(.05,2*E)/2*k/N:.4f}-{__import__("scipy.stats").stats.chi2.ppf(.95,2*E+2)/2*k/N:.4f}); N@FDP_c 1 % {n1}')
