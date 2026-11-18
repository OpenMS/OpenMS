"""q4: entrapment counts behind FDP_c at q <= 0.01, and a Poisson 90 % interval, CUR vs NEW."""
import sqlite3, pandas as pd, numpy as np
from scipy.stats import chi2
L='/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/gap-20261114/armE_raw01/library.pqp'
c=sqlite3.connect(f'file:{L}?mode=ro',uri=True)
ent=pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""",c).set_index('precursor_id').e
tg=pd.read_sql('select ID precursor_id from PRECURSOR where DECOY=0',c).precursor_id
R=(ent.loc[tg]==1).sum()/(ent.loc[tg]==0).sum(); k=1+1/R
print(f'R={R:.4f}, FDP_c = entrapments x {k:.2f} / N')
tot={}
for raw in ('01','02','03'):
    for arm,root in (('cur','ms1' if raw!='03' else 'pg'),('new','pg')):
        Es=[];Ns=[]
        for s in range(5):
            t=pd.read_csv(f'{root}/raw{raw}_s{s}_{arm}/p.tsv.gz',sep='\t'); t=t[(t.is_decoy==0)&(t.q_value<=0.01)]
            Es.append(int(t.precursor_id.map(ent).fillna(0).sum())); Ns.append(len(t))
        E=sum(Es); N=sum(Ns); lo=chi2.ppf(.05,2*E)/2*k/N; hi=chi2.ppf(.95,2*E+2)/2*k/N
        tot.setdefault(arm,[0,0]); tot[arm][0]+=E; tot[arm][1]+=N
        print(f'raw{raw} {arm}: entrapments per seed {Es}, pooled FDP_c {E*k/N:.4f} (90 % {lo:.4f}-{hi:.4f})')
for arm,(E,N) in tot.items():
    print(f'all raws {arm}: pooled FDP_c {E*k/N:.4f} (90 % {chi2.ppf(.05,2*E)/2*k/N:.4f}-{chi2.ppf(.95,2*E+2)/2*k/N:.4f}), E={E}, N={N}')
