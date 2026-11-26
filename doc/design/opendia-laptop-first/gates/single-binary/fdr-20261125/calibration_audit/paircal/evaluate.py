"""Evaluation ONLY. This module is never imported by the label-free selector.
Usage: evaluate.py library.pqp selected.tsv out.tsv
Selected file columns: dataset, seed, rule, precursor_id, score, q_value, decoy.
"""
import sys
import sqlite3
import numpy as np
import pandas as pd
from scipy.stats import chi2

library, selected, out = sys.argv[1:]
with sqlite3.connect(f'file:{library}?mode=ro', uri=True) as c:
    labels = pd.read_sql("SELECT p.ID precursor_id, p.DECOY decoy, pe.UNMODIFIED_SEQUENCE seq, MIN(CASE WHEN INSTR(pr.PROTEIN_ACCESSION, '_p_target')>0 THEN 1 ELSE 0 END) e FROM PRECURSOR p JOIN PRECURSOR_PEPTIDE_MAPPING m ON m.PRECURSOR_ID=p.ID JOIN PEPTIDE pe ON pe.ID=m.PEPTIDE_ID JOIN PEPTIDE_PROTEIN_MAPPING pp ON pp.PEPTIDE_ID=m.PEPTIDE_ID JOIN PROTEIN pr ON pr.ID=pp.PROTEIN_ID GROUP BY p.ID", c)
audit = pd.read_csv('/home/singjc/.jcode/scratch/fdr-20261125/label_audit_hela.tsv', sep='\t')
exact = set(audit.seq)
targets = labels[labels.decoy == 0].copy()
targets['corrected'] = targets.e.where(~targets.seq.isin(exact), 0)
sel = pd.read_csv(selected, sep='\t', usecols=['precursor_id','score','decoy','q_value','dataset','seed','rule'], dtype={'precursor_id':'int32','score':'float64','decoy':'uint8','q_value':'float64','dataset':'category','seed':'uint8','rule':'category'})
results = []
for policy in ['e', 'corrected']:
    ratio = targets[policy].sum() / (len(targets)-targets[policy].sum())
    multiplier = 1+1/ratio
    emap = targets.set_index('precursor_id')[policy]
    for (dataset, seed, rule), g in sel.groupby(['dataset','seed','rule'], observed=True):
        t = g[g.decoy==0].sort_values(['q_value','score','precursor_id'], ascending=[True,False,True]).copy()
        t['e'] = t.precursor_id.map(emap)
        assert not t.e.isna().any()
        acc = t[t.q_value <= .01]
        n, e = len(acc), int(acc.e.sum())
        realized = np.cumsum(t.e.to_numpy())*multiplier/np.arange(1,len(t)+1)
        ok = np.flatnonzero(realized <= .01)
        results.append(dict(policy=policy,dataset=dataset,seed=seed,rule=rule,N=n,E=e,combined=multiplier*e/max(n,1),N1=int(ok[-1]+1) if len(ok) else 0,multiplier=multiplier))
r = pd.DataFrame(results)
r.to_csv(out, sep='\t', index=False)
print(r.to_string(index=False))
print('\nPOOLED')
for (policy,dataset,rule),g in r.groupby(['policy','dataset','rule']):
    n,e = g.N.sum(), g.E.sum()
    k = g.multiplier.iloc[0]
    lo = chi2.ppf(.05,2*e)/2 if e else 0
    hi = chi2.ppf(.95,2*e+2)/2
    print(f'{policy} {dataset} {rule}: N={n} E={e} combined={k*e/max(n,1):.6f} 90%CI=[{k*lo/max(n,1):.6f},{k*hi/max(n,1):.6f}] median_N1={g.N1.median():.0f} per_seed_N1={g.N1.tolist()}')
