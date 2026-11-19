"""PRECHECK_SCORER36 addendum: TPR at fixed decoy FPR per abundance quartile (tool score), r~7 D s0 default.
Positives: best tool-score window within 6 s of DIA-NN RT per DIA-NN non-entrapment q<=0.01 precursor.
Negatives: best tool-score window of each decoy precursor. Diagnostic only."""
import numpy as np, pandas as pd
G = '/home/singjc/.jcode/scratch/dia2/gap7'
F = pd.read_parquet(f'{G}/D36_feat.parquet', columns=['precursor_id', 'is_decoy', 'rt_apex', 'score'])
d = pd.read_parquet('/home/singjc/.jcode/scratch/dia2/entlib12/D/report.parquet', columns=['Modified.Sequence', 'Precursor.Charge', 'Q.Value', 'RT', 'Precursor.Quantity'])
import sqlite3
c = sqlite3.connect('file:/home/singjc/.jcode/scratch/dia2/entlib10/cut/hela_r5.predicted_decoys.pqp?mode=ro', uri=True)
tg = pd.read_sql("""select p.ID pid, p.CHARGE z, pe.MODIFIED_SEQUENCE ms, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID where p.DECOY=0 group by p.ID""", c)
tg = tg[tg.e == 0].assign(sk=tg.ms + '/' + tg.z.astype(str)).drop_duplicates('sk').set_index('sk')
d = d[d['Q.Value'] <= 0.01]
d['pid'] = tg.pid.reindex(d['Modified.Sequence'] + '/' + d['Precursor.Charge'].astype(str)).values
d = d.dropna(subset=['pid']).astype({'pid': int}).drop_duplicates('pid').set_index('pid')
s = np.where(F.score.values > -1e29, F.score.values, -np.inf)
F = F.assign(s=s)
neg = F[F.is_decoy == 1].groupby('precursor_id').s.max().values
neg = neg[np.isfinite(neg)]
T = F[(F.is_decoy == 0) & F.precursor_id.isin(d.index)].copy()
T = T[np.abs(T.rt_apex.values - d.RT.reindex(T.precursor_id).values * 60) <= 6]
P = T.groupby('precursor_id').s.max()
q = pd.qcut(np.log10(d['Precursor.Quantity'].reindex(P.index).clip(lower=1)), 4, labels=['Q1', 'Q2', 'Q3', 'Q4'])
print(f'decoy precursors {len(neg)}, DIA-NN-confirmed windows {len(P)} (of {len(d)} DIA-NN IDs)')
for fpr in (0.0001, 0.001, 0.01, 0.05):
    thr = np.quantile(neg, 1 - fpr)
    print(f'decoy FPR {fpr:5.3f}: TPR', ' '.join(f'{k} {(P[q == k] > thr).mean():.3f}' for k in ['Q1', 'Q2', 'Q3', 'Q4']), f'all {(P > thr).mean():.3f}')
print('DIA-NN IDs with no window within 6 s:', len(d) - len(P))
miss = d.index.difference(P.index)
anyw = F[(F.is_decoy == 0) & F.precursor_id.isin(miss)].precursor_id.unique()
qa = pd.qcut(np.log10(d['Precursor.Quantity'].clip(lower=1)), 4, labels=['Q1', 'Q2', 'Q3', 'Q4'])
print(f'  of these: no window at all {len(miss) - len(anyw)}, only windows >6 s away {len(anyw)}')
print('  miss share per DIA-NN quartile:', ' '.join(f'{k} {qa.reindex(miss).eq(k).sum() / qa.eq(k).sum():.3f}' for k in ['Q1', 'Q2', 'Q3', 'Q4']))
thr = np.quantile(neg, 1 - 1e-4)
print(f'expected IDs = TPR(1e-4) x confirmed = {(P > thr).sum()} (tool N@FDP_c 1 % 11,258)')
