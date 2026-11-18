"""Which evidence separates our rejected right peaks (R) from decoys, in OpenSWATH's own feature scores?
R = SVM q<=0.01 IDs whose right peak OpenDIA scores below its cut. Compare OpenSWATH sub-scores on the same precursor's
SVM feature for R vs SVM-scored decoys at matched SVM d-score, and vs IDs OpenDIA accepts (A).
Also: our 24 features for R vs score-matched decoys (where our scorer sees no difference by construction)."""
import sqlite3, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
R = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/openswath-hela-published'
B = '/home/singjc/.jcode/scratch/dia2'
fs = pd.read_csv(f'{B}/img/funnel_imgrid_svm.tsv', sep='\t')
c = sqlite3.connect(f'file:{R}/svm.osw?mode=ro', uri=True)
cols = [r[1] for r in c.execute('pragma table_info(FEATURE_MS2)') if r[1].startswith('VAR_')]
c1 = [r[1] for r in c.execute('pragma table_info(FEATURE_MS1)') if r[1].startswith('VAR_')]
q = f"""select f.PRECURSOR_ID pid, p.DECOY decoy, s.SCORE d, s.QVALUE q, {','.join('m.'+x for x in cols)}, {','.join('m1.'+x+' as MS1_'+x for x in c1)}
        from SCORE_MS2 s join FEATURE f on f.ID=s.FEATURE_ID join PRECURSOR p on p.ID=f.PRECURSOR_ID
        join FEATURE_MS2 m on m.FEATURE_ID=f.ID left join FEATURE_MS1 m1 on m1.FEATURE_ID=f.ID where s.RANK=1"""
X = pd.read_sql(q, c)
allc = cols + ['MS1_' + x for x in c1]
Rset = set(fs[fs.g == 'right peak below cut'].pid); Aset = set(fs[fs.g == 'accepted'].pid)
Rx = X[(X.decoy == 0) & X.pid.isin(Rset)]; Ax = X[(X.decoy == 0) & X.pid.isin(Aset)]
dec = X[X.decoy == 1].sort_values('d'); ds = dec.d.values
# Decoys at the SVM score of R (score-matched in the reference's own score): shows what SVM features carry beyond our score.
print(f'R {len(Rx)}, A {len(Ax)}, SVM decoys {len(dec)}; SVM d-score R median {Rx.d.median():.2f}, A {Ax.d.median():.2f}, decoys {dec.d.median():.2f}')
y = np.r_[np.ones(len(Rx)), np.zeros(len(dec))]
rows = []
for k in allc:
    a = np.r_[Rx[k].fillna(0).values, dec[k].fillna(0).values]
    if np.nanstd(a) == 0: continue
    au = roc_auc_score(y, a); rows.append((k, round(max(au, 1 - au), 3), '+' if au >= 0.5 else '-'))
print('OpenSWATH sub-score AUC, R vs all SVM decoys (top 20):')
print(pd.DataFrame(rows, columns=['score', 'AUC', 'dir']).sort_values('AUC', ascending=False).head(20).to_string(index=False))
