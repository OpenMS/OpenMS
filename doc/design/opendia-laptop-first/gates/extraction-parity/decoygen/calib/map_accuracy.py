"""Evaluation only: accuracy of calibration RT maps against DIA-NN 2.3.1's run RTs (ProteoBench reference report, Q <= 0.01,
human targets). Each DIA-NN precursor is matched to the PAIRLIB target library by stripped sequence + charge (unmodified
precursors only); the map predicts run RT from library RT. Reports |predicted - DIA-NN RT| quantiles and the share within
the map's own 4-sigma tolerance. Usage: map_accuracy.py <cal.json> ... """
import json, sys, sqlite3, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
ref = pd.read_parquet(B.replace('github/OpenMS-r1-benchmark-20261107', 'datasets/proteobench-astral-entrapment/proteobench_diann_v231_res_20260717_123332/input_file.parquet'),
                      columns=['Run', 'Stripped.Sequence', 'Modified.Sequence', 'Precursor.Charge', 'RT', 'Q.Value', 'Decoy'])
ref = ref[(ref['Q.Value'] <= 0.01) & (ref.Decoy == 0) & (ref['Modified.Sequence'] == ref['Stripped.Sequence'].str.replace('C', 'C(UniMod:4)'))]
c = sqlite3.connect(f'file:{B}/pairlib-20261113/lib/library.pqp?mode=ro', uri=True)
lib = pd.read_sql("""select pe.UNMODIFIED_SEQUENCE seq, pe.MODIFIED_SEQUENCE ms, p.CHARGE z, p.LIBRARY_RT lrt from PRECURSOR p
  join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID where p.DECOY=0""", c)
lib = lib[~lib.ms.str.contains('UniMod:35')].drop_duplicates(['seq', 'z'])
for path in sys.argv[1:]:
    raw = '01' if 'raw01' in path else '02'
    r = ref[ref.Run.str.endswith('_' + raw)].merge(lib, left_on=['Stripped.Sequence', 'Precursor.Charge'], right_on=['seq', 'z'])
    j = json.load(open(path))['rt']
    pred = np.interp(r.lrt, j['knots_library_rt'], j['knots_run_rt_s'])
    e = np.abs(pred - r.RT.values * 60.0)
    tol = 4 * j['sigma_s']
    print(f'{path:45s} n {len(e)} |err| p50 {np.median(e):6.2f} p95 {np.quantile(e, .95):6.2f} p99 {np.quantile(e, .99):7.2f} s | '
          f'within own 4 sigma ({tol:5.1f} s) {np.mean(e <= tol):.4f}')
