"""PLAN_RTCENTRE evaluation: N, agreement with DIA-NN, calibrated-map accuracy, for runs under a directory.
Usage: eval_rtc.py <runs dir> <prefix D|H> <diann report.parquet>"""
import sys, json, glob, sqlite3, numpy as np, pandas as pd
R, P, REF = sys.argv[1], sys.argv[2], sys.argv[3]
H = '/home/singjc/Documents/github/testing_auto_irt/diaPASEF_evosep_published'
lib = sqlite3.connect(f'file:{H}/200623_hela_evosep_fractions_library_ptypic_decoy.pqp?mode=ro', uri=True)
pre = pd.read_sql("select p.ID pid, p.CHARGE z, p.DECOY decoy, p.LIBRARY_RT lrt, pe.MODIFIED_SEQUENCE ms from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID", lib)
pre['key'] = pre.ms + '/' + pre.z.astype(str)
key_of = dict(zip(pre.pid, pre.key)); tg = pre[pre.decoy == 0].set_index('key')
dn = pd.read_parquet(REF, columns=['Modified.Sequence', 'Precursor.Charge', 'Q.Value', 'RT'])
dn['key'] = dn['Modified.Sequence'] + '/' + dn['Precursor.Charge'].astype(str)
dn = dn[dn['Q.Value'] <= 0.01].drop_duplicates('key')
D = set(dn.key); dl = dn[dn.key.isin(tg.index)].copy(); dl['lrt'] = tg.loc[dl.key, 'lrt'].values
rows = []
for d in sorted(glob.glob(f'{R}/{P}_*_s*')):
    if not d.split('/')[-1].count('_s'): continue
    import os
    if not os.path.isdir(d): continue
    arm, seed = d.split('/')[-1].split('_')[1], d.split('_s')[-1]
    t = pd.read_csv(f'{d}/precursors.tsv', sep='\t')
    O = set(t[t.q_value <= 0.01].precursor_id.map(key_of))
    cal = json.load(open(glob.glob(f'{d}/*.calibration.json')[0]))['rt']
    pred = np.interp(dl.lrt, cal['knots_library_rt'], cal['knots_run_rt_s'])
    err = np.abs(dl.RT.values * 60 - pred)
    rows.append(dict(arm=arm, seed=int(seed), N=len(O), agree=len(O & D) / len(O), map_p95=np.quantile(err, .95), map_p50=np.median(err),
                     early_err=np.median((dl.RT.values * 60 - pred)[dl.lrt.values < 12])))
df = pd.DataFrame(rows).sort_values(['arm', 'seed'])
print(df.round(4).to_string(index=False))
print(df.groupby('arm').agg(N_med=('N', 'median'), N_min=('N', 'min'), N_max=('N', 'max'), agree_med=('agree', 'median'),
                            map_p95=('map_p95', 'median'), early=('early_err', 'median')).round(4).to_string())
