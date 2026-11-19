"""PRECHECK_IMPICK part 2: why is the right peak not a candidate for 15.8 % of DIA-NN IDs (r~1.56 D s0, default)?
Checks, with DIA-NN RT and 1/K0 as evaluation truth: (a) DIA-NN RT outside the calibrated RT window, (b) DIA-NN 1/K0
outside the calibrated 1/K0 window (calibration JSON of the run), (c) candidate cap (5 peaks) reached, (d) other."""
import glob, json, sqlite3, numpy as np, pandas as pd
W = '/home/singjc/.jcode/scratch/dia2/gap7/ev/D_new_s0'
L = '/home/singjc/.jcode/scratch/dia2/ctest_cache/c/hela_entrap_codon.predicted_decoys.pqp'
DN = '/home/singjc/.jcode/scratch/dia2/entlib12/D/report.parquet'
F = pd.read_csv(glob.glob(f'{W}/*.features.tsv')[0], sep='\t', usecols=['precursor_id', 'is_decoy', 'rt_apex', 'abs_library_rt_error', 'score', 'q_value'])
J = json.load(open(glob.glob(f'{W}/*.calibration.json')[0]))
rt_knots_x, rt_knots_y = J['rt']['knots_library_rt'], J['rt']['knots_run_rt_s']
tol_rt = J['rt']['tolerance_half_width_s']
im_a, im_b, tol_im = J['im']['intercept_vssc'], J['im']['slope'], J['im']['tolerance_half_width_vssc']
print(f'RT tolerance +-{tol_rt:.1f} s, 1/K0 tolerance +-{tol_im:.4f}, 1/K0 map {im_a:.4f} + {im_b:.4f} x')
c = sqlite3.connect(f'file:{L}?mode=ro', uri=True)
tg = pd.read_sql("""select p.ID pid, p.CHARGE z, p.LIBRARY_RT lrt, p.LIBRARY_DRIFT_TIME lim, pe.MODIFIED_SEQUENCE ms,
 min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
 from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
 join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID where p.DECOY=0 group by p.ID""", c)
tg['key'] = tg.ms + '/' + tg.z.astype(str); lab = tg[tg.e == 0].drop_duplicates('key').set_index('key')
d = pd.read_parquet(DN, columns=['Modified.Sequence', 'Precursor.Charge', 'Q.Value', 'RT', 'IM', 'Precursor.Quantity'])
d = d[d['Q.Value'] <= 0.01]; d['key'] = d['Modified.Sequence'] + '/' + d['Precursor.Charge'].astype(str)
d = d[d.key.isin(lab.index)].drop_duplicates('key')
d = d.join(lab[['pid', 'lrt', 'lim']], on='key').set_index('pid')
# Monotone piecewise-linear map (the JSON knots; ends extrapolated with the end slopes).
x, y = np.array(rt_knots_x), np.array(rt_knots_y)
def rtmap(v):
    out = np.interp(v, x, y)
    lo, hi = v < x[0], v > x[-1]
    out[lo] = y[0] + (v[lo] - x[0]) * J['rt']['low_slope_s']; out[hi] = y[-1] + (v[hi] - x[-1]) * J['rt']['high_slope_s']
    return out
d['mapped_rt'] = rtmap(d.lrt.values); d['rt_err'] = d.RT * 60 - d.mapped_rt
d['mapped_im'] = im_a + im_b * d.lim; d['im_err'] = d.IM - d.mapped_im
T = F[(F.is_decoy == 0) & F.precursor_id.isin(d.index)].copy()
T['right'] = np.abs(T.rt_apex - (d.RT * 60).reindex(T.precursor_id).values) <= 6
g = T.groupby('precursor_id')
d['n'] = g.size(); d['any_right'] = g.right.any()
nc = d[~d.any_right]
print(f'DIA-NN IDs {len(d)}; right peak not a candidate {len(nc)} ({len(nc) / len(d):.3f})')
a = np.abs(nc.rt_err) > tol_rt; b = np.abs(nc.im_err) > tol_im; cap = nc.n >= 5
print(f'  (a) DIA-NN RT outside the RT window          {a.mean():.3f}')
print(f'  (b) DIA-NN 1/K0 outside the 1/K0 window      {b.mean():.3f}  (all DIA-NN IDs: {(np.abs(d.im_err) > tol_im).mean():.3f})')
print(f'  (c) inside both windows, 5-candidate cap hit {(~a & ~b & cap).mean():.3f}')
print(f'  (d) inside both windows, below the cap       {(~a & ~b & ~cap).mean():.3f}')
print(f'  |1/K0 error| p50 / p95: all DIA-NN IDs {np.median(np.abs(d.im_err)):.4f} / {np.quantile(np.abs(d.im_err), 0.95):.4f}; '
      f'not-a-candidate {np.median(np.abs(nc.im_err)):.4f} / {np.quantile(np.abs(nc.im_err), 0.95):.4f}')
print(f'  |RT error| p50 / p95 s: all {np.median(np.abs(d.rt_err)):.1f} / {np.quantile(np.abs(d.rt_err), 0.95):.1f}; '
      f'not-a-candidate {np.median(np.abs(nc.rt_err)):.1f} / {np.quantile(np.abs(nc.rt_err), 0.95):.1f}')
nc_q = pd.qcut(np.log10(d['Precursor.Quantity'].clip(lower=1)), 4, labels=['Q1', 'Q2', 'Q3', 'Q4'])
print('  not-a-candidate share by DIA-NN quantity quartile:', (~d.any_right).groupby(nc_q, observed=True).mean().round(3).to_dict())
