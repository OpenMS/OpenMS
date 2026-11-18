"""HeLa remaining-gap funnel on the current code (963543c7f0, seed 0, D 2737.d, features written).
References: DIA-NN 2.2 and OpenSWATH + PyProphet SVM (same library), q <= 0.01.
For each reference ID we do not report: not scored, wrong peak (our best candidate > tol from the reference RT),
right peak below our cut. Also: is the reference's peak among our candidates at all, and abundance/score of each group."""
import sqlite3, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
H = '/home/singjc/Documents/github/testing_auto_irt/diaPASEF_evosep_published'
B = '/home/singjc/.jcode/scratch/dia2'
R = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/openswath-hela-published'
RUN = f'{B}/img/pre_imgrid'
lib = sqlite3.connect(f'file:{H}/200623_hela_evosep_fractions_library_ptypic_decoy.pqp?mode=ro', uri=True)
pre = pd.read_sql("select p.ID pid, p.CHARGE z, p.DECOY decoy, p.LIBRARY_RT lrt, pe.MODIFIED_SEQUENCE ms from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID", lib)
pre['key'] = pre.ms + '/' + pre.z.astype(str)
tg = pre[pre.decoy == 0].set_index('key')
F = pd.read_csv(f'{RUN}/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.features.tsv', sep='\t')
win = F.sort_values('score', ascending=False).drop_duplicates('precursor_id').set_index('precursor_id')
acc = set(pd.read_csv(f'{RUN}/precursors.tsv', sep='\t').query('q_value <= 0.01').precursor_id)
cut = win.loc[list(acc)].score.min()
print(f'OpenDIA: {len(acc)} at q<=0.01, score cut {cut:.3f}; scored precursors {len(win)} ({int((win.is_decoy==0).sum())} targets)')

dn = pd.read_parquet(f'{B}/diann_hela/diann_report.parquet', columns=['Modified.Sequence', 'Precursor.Charge', 'Q.Value', 'RT', 'Precursor.Quantity'])
dn['key'] = dn['Modified.Sequence'] + '/' + dn['Precursor.Charge'].astype(str)
dn = dn[(dn['Q.Value'] <= 0.01) & dn.key.isin(tg.index)].drop_duplicates('key')
dn = pd.DataFrame({'pid': tg.loc[dn.key, 'pid'].values, 'rt': dn.RT.values * 60, 'qty': dn['Precursor.Quantity'].values})
c = sqlite3.connect(f'file:{R}/svm.osw?mode=ro', uri=True)
sv = pd.read_sql("""select f.PRECURSOR_ID pid, f.EXP_RT rt, s.QVALUE q from SCORE_MS2 s join FEATURE f on f.ID=s.FEATURE_ID
                    join PRECURSOR p on p.ID=f.PRECURSOR_ID where p.DECOY=0 and s.RANK=1""", c)
sv = sv[sv.q <= 0.01].drop_duplicates('pid')
TOL = 3.43
def funnel(name, ref):
    ref = ref.copy()
    ref['scored'] = ref.pid.isin(win.index)
    ref['our_rt'] = ref.pid.map(win.rt_apex)
    ref['acc'] = ref.pid.isin(acc)
    m = F[F.precursor_id.isin(ref.pid)].merge(ref[['pid', 'rt']], left_on='precursor_id', right_on='pid')
    m['hit'] = (m.rt_apex - m.rt).abs() <= TOL
    has = set(m[m.hit].precursor_id)
    rs = m[m.hit].sort_values('score', ascending=False).drop_duplicates('precursor_id').set_index('precursor_id').score
    ref['right_cand'] = ref.pid.isin(has)
    ref['right_score'] = ref.pid.map(rs)
    ref['g'] = np.select([ref.acc, ~ref.scored, (ref.our_rt - ref.rt).abs() <= TOL, ref.right_cand],
                         ['accepted', 'not scored', 'right peak below cut', 'wrong peak, right one is a candidate'], 'wrong peak, right one never a candidate')
    print(f'\n== {name}: {len(ref)} IDs in the library')
    print(ref.g.value_counts().to_string())
    rb = ref[ref.g == 'right peak below cut'].right_score
    print(f'  right peak below cut: score quantiles {rb.quantile([.25,.5,.75,.9]).round(2).to_dict()} (cut {cut:.2f})')
    wc = ref[ref.g == 'wrong peak, right one is a candidate']
    print(f'  wrong peak with right candidate: right-candidate score median {wc.right_score.median():.2f}, would pass cut {(wc.right_score >= cut).mean():.3f}')
    return ref
fd = funnel('DIA-NN q<=0.01', dn)
fs = funnel('OpenSWATH SVM q<=0.01', sv)
both = set(dn.pid) & set(sv.pid)
print(f'\nDIA-NN and SVM both: {len(both)}; OpenDIA has {len(both & acc)} ({len(both & acc)/len(both):.3f})')
gb = fd[fd.pid.isin(both)].g.value_counts()
print('of those, where OpenDIA loses them:\n' + gb.to_string())
fd['qp'] = np.log10(fd.qty.clip(lower=1)).rank(pct=True)
print('\nDIA-NN abundance percentile median by group:\n' + fd.groupby('g').qp.median().round(2).to_string())
fd.to_csv(f'{B}/img/funnel_imgrid_diann.tsv', sep='\t', index=False); fs.to_csv(f'{B}/img/funnel_imgrid_svm.tsv', sep='\t', index=False)
