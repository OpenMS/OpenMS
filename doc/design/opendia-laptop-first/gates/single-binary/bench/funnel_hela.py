"""HeLa funnel against DIA-NN (diagnostic). For DIA-NN's q <= 0.01 precursors that OpenDIA does not report:
not in the population (never scored), scored on the wrong peak (|our apex - DIA-NN RT| > tol), or right peak but below
our q <= 0.01 cut. Also: RT error of our right-peak winners, decoy/target score picture, and DIA-NN RT vs library RT
(how wide the RT window must be). Uses the run-level precursors.tsv (all scored windows, best per precursor)."""
import sqlite3, numpy as np, pandas as pd, json
H = '/home/singjc/Documents/github/testing_auto_irt/diaPASEF_evosep_published'
B = '/home/singjc/.jcode/scratch/e2e/hela'
lib = sqlite3.connect(f'file:{H}/200623_hela_evosep_fractions_library_ptypic_decoy.pqp?mode=ro', uri=True)
pre = pd.read_sql("select p.ID pid, p.CHARGE z, p.DECOY decoy, p.LIBRARY_RT lrt, pe.MODIFIED_SEQUENCE ms from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID", lib)
pre['key'] = pre.ms + '/' + pre.z.astype(str)
tg = pre[pre.decoy == 0].set_index('key')
run = pd.read_csv(f'{B}/opendia/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.precursors.tsv', sep='\t')
print('run rows', len(run), 'targets', int((run.is_decoy == 0).sum()), 'decoys', int((run.is_decoy == 1).sum()))
best = run.sort_values('score', ascending=False).drop_duplicates('precursor_id').set_index('precursor_id')
dn = pd.read_parquet(f'{B}/diann_report.parquet', columns=['Modified.Sequence', 'Precursor.Charge', 'Q.Value', 'RT', 'RT.Start', 'RT.Stop', 'Precursor.Quantity', 'IM'])
dn['key'] = dn['Modified.Sequence'] + '/' + dn['Precursor.Charge'].astype(str)
dn = dn[(dn['Q.Value'] <= 0.01) & dn.key.isin(tg.index)].drop_duplicates('key')
dn['pid'] = tg.loc[dn.key, 'pid'].values
dn['lrt'] = tg.loc[dn.key, 'lrt'].values
cal = json.load(open(f'{B}/opendia/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.calibration.json'))
print('calibration keys', list(cal.keys())[:20])
dn['in_pop'] = dn.pid.isin(best.index)
dn['our_rt'] = dn.pid.map(best.rt_apex)
dn['our_q'] = dn.pid.map(best.q_value)
dn['dt'] = dn.our_rt - dn.RT * 60
fw = ((dn['RT.Stop'] - dn['RT.Start']) * 60).median()
tol = max(3.0, fw / 2)
dn['grp'] = np.where(~dn.in_pop, 'not_scored', np.where(dn.dt.abs() > tol, 'wrong_peak', np.where(dn.our_q <= 0.01, 'accepted', 'right_below_cut')))
print(f'DIA-NN q<=0.01 in library: {len(dn)}; median DIA-NN peak width {fw:.2f} s, right-peak tol {tol:.2f} s')
print(dn.grp.value_counts().to_string())
dn['lq'] = np.log10(dn['Precursor.Quantity'].clip(lower=1)); dn['qp'] = dn.lq.rank(pct=True)
print(dn.groupby('grp').qp.median().round(3).to_string())
print('wrong-peak |dt| quantiles (s):', dn[dn.grp == 'wrong_peak'].dt.abs().quantile([.25, .5, .75, .9]).round(1).to_dict())
# DIA-NN RT against library RT mapped by our calibration: residual spread tells the window needed.
lr = np.polyfit(dn.lrt, dn.RT * 60, 3); res = dn.RT * 60 - np.polyval(lr, dn.lrt)
print(f'DIA-NN RT vs library RT (cubic): residual p50 {np.median(np.abs(res)):.1f} s, p95 {np.quantile(np.abs(res), .95):.1f} s, p99 {np.quantile(np.abs(res), .99):.1f} s')
print('run span s', run.rt_apex.min(), run.rt_apex.max())
# Not scored: why? candidate RT window vs DIA-NN RT
ns = dn[dn.grp == 'not_scored']
print('not scored: residual |res| quantiles', np.abs(res[dn.grp == 'not_scored']).quantile([.5, .9]).round(1).to_dict())
dn.to_csv(f'{B}/funnel_hela.tsv', sep='\t', index=False)
