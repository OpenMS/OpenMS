"""PRECHECK_VALLEY: r~1.56 D s0, default (ev/D_new_s0, features) vs -valley_fraction 0.8 -max_peaks 8 (ev/v08m8).
For DIA-NN non-entrapment q<=0.01 IDs: right candidate present (apex within 6 s), accepted at the tool's q<=0.01
with a right window. Evaluation labels only."""
import glob, sqlite3, numpy as np, pandas as pd
G = '/home/singjc/.jcode/scratch/dia2/gap7/ev'
L = '/home/singjc/.jcode/scratch/dia2/ctest_cache/c/hela_entrap_codon.predicted_decoys.pqp'
DN = '/home/singjc/.jcode/scratch/dia2/entlib12/D/report.parquet'
c = sqlite3.connect(f'file:{L}?mode=ro', uri=True)
tg = pd.read_sql("""select p.ID pid, p.CHARGE z, pe.MODIFIED_SEQUENCE ms, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
 from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
 join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID where p.DECOY=0 group by p.ID""", c)
lab = tg[tg.e == 0].assign(k=tg.ms + '/' + tg.z.astype(str)).drop_duplicates('k').set_index('k')
d = pd.read_parquet(DN, columns=['Modified.Sequence', 'Precursor.Charge', 'Q.Value', 'RT', 'Precursor.Quantity'])
d = d[d['Q.Value'] <= 0.01]; k = d['Modified.Sequence'] + '/' + d['Precursor.Charge'].astype(str)
d = d[k.isin(lab.index).values]; d['pid'] = lab.pid.reindex(k[k.isin(lab.index)]).values; d = d.drop_duplicates('pid').set_index('pid')
d['qq'] = pd.qcut(np.log10(d['Precursor.Quantity'].clip(lower=1)), 4, labels=['Q1', 'Q2', 'Q3', 'Q4'])
out = {}
for arm in ('D_new_s0', 'v08m8'):
    F = pd.read_csv(glob.glob(f'{G}/{arm}/*.features.tsv')[0], sep='\t', usecols=['precursor_id', 'is_decoy', 'rt_apex', 'score', 'q_value'])
    T = F[(F.is_decoy == 0) & F.precursor_id.isin(d.index)].copy()
    T['right'] = np.abs(T.rt_apex - d.RT.reindex(T.precursor_id).values * 60) <= 6
    T['kept'] = T.score > -1e29
    g = T.groupby('precursor_id')
    s = pd.DataFrame({'any_right': g.right.any(), 'acc_right': g.apply(lambda x: bool((x.kept & x.right & (x.q_value <= 0.01)).any()))}).reindex(d.index).fillna(False)
    out[arm] = s
    print(f'{arm}: windows {len(F)}; DIA-NN IDs {len(d)}: right candidate {s.any_right.mean():.3f}, accepted with right window {s.acc_right.mean():.3f}')
    print('   right candidate by quartile', s.any_right.groupby(d.qq).mean().round(3).to_dict())
a, b = out['D_new_s0'], out['v08m8']
gain = ~a.any_right & b.any_right
print(f'gained a right candidate: {gain.sum()}; of these accepted (right window) {b.acc_right[gain].mean():.3f}; '
      f'by quartile {gain.groupby(d.qq).sum().to_dict()}')
lost = a.acc_right & ~b.acc_right; new = ~a.acc_right & b.acc_right
print(f'accepted-right: lost {lost.sum()}, gained {new.sum()} (of which via a newly present right candidate {(new & gain).sum()})')
