"""Which candidate evidence separates the right-peak-but-rejected DIA-NN IDs (R) from decoy winners at the same score,
without also separating entrapment winners from those decoys? Evaluation only: labels pick the groups, never a method.
Matching: each R / entrapment winner (score below the cut) gets the 3 decoy winners nearest in D11FDWT score.
Reported per feature: AUC(R vs matched decoys) and AUC(entrapment vs matched decoys); useful = first well above 0.5
while the second stays near 0.5. Features (all label-free, from the arm-E arrays, central bins 10:22 = apex +-5.6 s):
  ms1_corr        Pearson r of MS1 M+0 with the summed fragment profile
  ms1_iso_corr    Pearson r of M+1 with M+0
  ms1_logint      log10(1 + central M+0 area)
  ppm_mean_abs    |intensity-weighted mean fragment ppm| (calibrated), over fragments with signal
  ppm_frag_sd     SD across fragments of their weighted mean ppm (>= 2 fragments with signal)
  frag_logint     log10(1 + central summed fragment area)
  frag_nz3        fragments with >= 3 nonzero central bins
  top_share       largest fragment's share of the central area (one-fragment dominance)
  others_corr     mean r of each fragment with the sum of the other fragments, fragments with signal only"""
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
D = B + '/gap-20261114/armE_raw01'
exec(open(B + '/gap-20261114/funnel.py').read().split("# Why stage 1 loses")[0].split('print(\'\\nFunnel')[0])
wr = w.set_index('precursor_id')
R = wr.loc[sorted(near_is_winner - accepted)]
E = w[(w.is_decoy == 0) & w.precursor_id.map(ent).eq(1) & (w.s < cut)]
Dw = w[w.is_decoy == 1].sort_values('s')
ds = Dw.s.values
def match(g, k=3):
    i = np.searchsorted(ds, g.s.values)
    idx = np.clip(np.r_[[i + o for o in range(-k // 2 + 1, k // 2 + 2)]].T, 0, len(ds) - 1)[:, :k].ravel()
    return Dw.row.values[idx]
grp = {'R': R.row.values, 'R_dec': match(R), 'E': E.row.values, 'E_dec': match(E)}
allrows = np.unique(np.concatenate(list(grp.values())))
XIC = np.load(f'{D}/xic.npy', mmap_mode='r'); MS1 = np.load(f'{D}/ms1.npy', mmap_mode='r'); PPM = np.load(f'{D}/ppm.npy', mmap_mode='r')
x = np.asarray(XIC[allrows], np.float64)[:, :, 10:22]; m = np.asarray(MS1[allrows], np.float64)[:, :, 10:22]
p = np.asarray(PPM[allrows], np.float64)[:, :, 10:22]
def corr(a, b):
    a = a - a.mean(-1, keepdims=True); b = b - b.mean(-1, keepdims=True)
    return (a * b).sum(-1) / np.sqrt((a * a).sum(-1) * (b * b).sum(-1) + 1e-12)
tot = x.sum(1); area = x.sum(2); has = area > 0
fppm = np.where(has, p.sum(2) / np.maximum(area, 1e-12), np.nan)
wmean = np.nansum(fppm * area, 1) / np.maximum(area.sum(1), 1e-12)
F = pd.DataFrame({
    'ms1_corr': corr(m[:, 0], tot), 'ms1_iso_corr': corr(m[:, 1], m[:, 0]), 'ms1_logint': np.log10(1 + m[:, 0].sum(1)),
    'ppm_mean_abs': np.abs(wmean), 'ppm_frag_sd': np.where(has.sum(1) >= 2, np.nanstd(fppm, 1), np.nan),
    'frag_logint': np.log10(1 + tot.sum(1)), 'frag_nz3': ((x > 0).sum(2) >= 3).sum(1),
    'top_share': area.max(1) / np.maximum(area.sum(1), 1e-12),
    'others_corr': np.nanmean(np.where(has, corr(x, tot[:, None, :] - x), np.nan), 1)}, index=allrows)
def auc(a, b, f):
    va, vb = F.loc[a, f].values, F.loc[b, f].values
    va, vb = va[np.isfinite(va)], vb[np.isfinite(vb)]
    return roc_auc_score(np.r_[np.ones(len(va)), np.zeros(len(vb))], np.r_[va, vb])
print(f'R {len(R)} (score p50 {R.s.median():.1f}), entrapment winners below cut {len(E)} (score p50 {E.s.median():.1f}); 3 score-matched decoys each')
print(f'{"feature":14s} {"AUC R>dec":>10s} {"AUC E>dec":>10s} {"R p50":>8s} {"Rdec p50":>9s}')
for f in F.columns:
    print(f'{f:14s} {auc(grp["R"], grp["R_dec"], f):10.3f} {auc(grp["E"], grp["E_dec"], f):10.3f} {np.nanmedian(F.loc[grp["R"], f]):8.3f} {np.nanmedian(F.loc[grp["R_dec"], f]):9.3f}')
# Same, but E restricted to the score range of R (E below cut is dominated by low scores)
Es = E[E.s >= R.s.quantile(0.1)]
print(f'\nE restricted to R score range (>= R p10 {R.s.quantile(0.1):.1f}): n {len(Es)}')
ge = match(Es)
for f in F.columns:
    print(f'  {f:14s} AUC E>dec {auc(Es.row.values, ge, f):.3f}')
