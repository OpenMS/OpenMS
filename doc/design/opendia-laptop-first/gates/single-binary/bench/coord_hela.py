"""HeLa: for DIA-NN IDs, compare DIA-NN RT / IM with our calibrated library coordinates.
RT map: isotonic knots from calibration.json (linear interpolation, end slopes). IM: intercept + slope * library IM."""
import json, sqlite3, numpy as np, pandas as pd
B = '/home/singjc/.jcode/scratch/e2e/hela'
H = '/home/singjc/Documents/github/testing_auto_irt/diaPASEF_evosep_published'
c = json.load(open(f'{B}/opendia/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.calibration.json'))
kx, ky = np.array(c['rt']['knots_library_rt']), np.array(c['rt']['knots_run_rt_s'])
def rtmap(x):
    y = np.interp(x, kx, ky)
    y = np.where(x < kx[0], ky[0] + (x - kx[0]) * c['rt']['low_slope_s'], y)
    return np.where(x > kx[-1], ky[-1] + (x - kx[-1]) * c['rt']['high_slope_s'], y)
f = pd.read_csv(f'{B}/funnel_hela.tsv', sep='\t')
lib = sqlite3.connect(f'file:{H}/200623_hela_evosep_fractions_library_ptypic_decoy.pqp?mode=ro', uri=True)
im = pd.read_sql("select ID pid, LIBRARY_DRIFT_TIME lim from PRECURSOR", lib).set_index('pid').lim
f['pred_rt'] = rtmap(f.lrt.values)
f['drt'] = f.RT * 60 - f.pred_rt
f['lim'] = f.pid.map(im)
f['cim'] = c['im']['intercept_vssc'] + c['im']['slope'] * f.lim
f['dim_cal'] = f.IM - f.cim
f['dim_raw'] = f.IM - f.lim
F = pd.read_csv(f'{B}/feat/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.features.tsv', sep='\t', usecols=['precursor_id', 'rt_apex'])
cand = F.merge(f[['pid', 'RT']], left_on='precursor_id', right_on='pid')
present = cand[(cand.rt_apex - cand.RT * 60).abs() <= 3.43].precursor_id.unique()
f['never'] = (f.grp == 'wrong_peak') & ~f.pid.isin(present)
tol_rt, tol_im = c['rt']['tolerance_half_width_s'], c['im']['tolerance_half_width_vssc']
for g, d in [('accepted', f[f.grp == 'accepted']), ('right_below_cut', f[f.grp == 'right_below_cut']), ('wrong: in candidates', f[(f.grp == 'wrong_peak') & ~f.never]), ('wrong: never a candidate', f[f.never])]:
    print(f'{g:26s} n {len(d):5d} | |dRT| p50 {d.drt.abs().median():5.2f} p95 {d.drt.abs().quantile(.95):5.2f} out of RT window {(d.drt.abs() > tol_rt).mean():.3f}'
          f' | dIM cal p50 {d.dim_cal.median():+.4f} |dIM| p95 {d.dim_cal.abs().quantile(.95):.4f} out of IM tol {(d.dim_cal.abs() > tol_im).mean():.3f} | dIM uncal p50 {d.dim_raw.median():+.4f} out of +-0.03 uncal {(d.dim_raw.abs() > 0.03).mean():.3f}')
