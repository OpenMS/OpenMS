"""PRECHECK_IM2D (diagnostic): RT x mobility picking (-im_peaks 8) against the current default, r~1.56 HeLa D seed 0.
Funnel of DIA-NN non-entrapment q<=0.01 IDs (evaluation only) for both runs, plus how often the added band candidates
are right / kept / accepted. Usage: precheck_im2d.py <default features dir> <im features dir>"""
import glob, sys, sqlite3, numpy as np, pandas as pd
A, B = sys.argv[1:3]
L = '/home/singjc/.jcode/scratch/dia2/ctest_cache/c/hela_entrap_codon.predicted_decoys.pqp'
DN = '/home/singjc/.jcode/scratch/dia2/entlib12/D/report.parquet'
c = sqlite3.connect(f'file:{L}?mode=ro', uri=True)
tg = pd.read_sql("""select p.ID pid, p.CHARGE z, pe.MODIFIED_SEQUENCE ms, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
 from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
 join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID where p.DECOY=0 group by p.ID""", c)
tg['key'] = tg.ms + '/' + tg.z.astype(str); lab = tg[tg.e == 0].drop_duplicates('key').set_index('key')
d = pd.read_parquet(DN, columns=['Modified.Sequence', 'Precursor.Charge', 'Q.Value', 'RT', 'Precursor.Quantity'])
d = d[d['Q.Value'] <= 0.01]; d['key'] = d['Modified.Sequence'] + '/' + d['Precursor.Charge'].astype(str)
d = d[d.key.isin(lab.index)].drop_duplicates('key'); d['pid'] = lab.pid.reindex(d.key).values; d = d.set_index('pid')
qq = pd.qcut(np.log10(d['Precursor.Quantity'].clip(lower=1)), 4, labels=['Q1', 'Q2', 'Q3', 'Q4'])

def funnel(path, name):
    F = pd.read_csv(glob.glob(f'{path}/*.features.tsv')[0], sep='\t', usecols=['candidate_row', 'precursor_id', 'is_decoy', 'rt_apex', 'score', 'q_value'])
    T = F[(F.is_decoy == 0) & F.precursor_id.isin(d.index)].copy()
    T['right'] = np.abs(T.rt_apex - (d.RT * 60).reindex(T.precursor_id).values) <= 6
    T['kept'] = T.score > -1e29
    g = T.groupby('precursor_id')
    s = pd.DataFrame({'any_right': g.right.any(), 'kept_right': g.apply(lambda x: bool((x.kept & x.right).any())),
                      'q': g.apply(lambda x: x.loc[x.kept, 'q_value'].min() if x.kept.any() else 1.0)}).reindex(d.index)
    s['any_right'] = s.any_right.fillna(False).astype(bool); s['kept_right'] = s.kept_right.fillna(False).astype(bool); s['q'] = s.q.fillna(1.0)
    acc = s.q <= 0.01
    print(f'== {name}: windows {len(F)} (targets {int((F.is_decoy == 0).sum())}, decoys {int(F.is_decoy.sum())})')
    print(f'  DIA-NN IDs with a right candidate {s.any_right.mean():.3f}; right kept {s.kept_right.mean():.3f}; accepted (right) {(acc & s.kept_right).mean():.3f}; accepted any {acc.mean():.3f}')
    print('  right candidate by quartile:', s.any_right.groupby(qq, observed=True).mean().round(3).to_dict())
    print('  accepted by quartile:      ', acc.groupby(qq, observed=True).mean().round(3).to_dict())
    return F, s

Fa, sa = funnel(A, 'default')
Fb, sb = funnel(B, '-im_peaks 8')
new_right = sb.any_right & ~sa.any_right
print(f'\nDIA-NN IDs that gain a right candidate: {new_right.sum()} ({new_right.mean():.3f}); of these kept {sb.kept_right[new_right].mean():.2f}, accepted {(sb.q[new_right] <= 0.01).mean():.2f}')
print(f'DIA-NN IDs that lose a right candidate: {(sa.any_right & ~sb.any_right).sum()}')
print(f'Accepted DIA-NN IDs: {int((sa.q <= 0.01).sum())} -> {int((sb.q <= 0.01).sum())}')
