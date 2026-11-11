"""Stage-D deep dive, step 1: per-window evidence from the raw extraction arrays (xic.npy, ppm.npy, ms1.npy, lib.npy,
xic12/lib12/orig12), descriptive only. DIA-NN and entrapment labels are used only to define the groups.
Groups (D11FWT, FDP_c 1 % cut, misses_D11FWT_<raw>.tsv):
  E   accepted reference precursors (their winning window)
  D   stage-D misses (winner = right window, below the cut)
  T1  top-1 % decoy winners;  T01 top-0.1 % decoy winners
Per window, central bins (apex +- 4.7 s = bins 12..19 of 32 over +- 15 s):
  n_sig        fragments (of 6) with any central signal
  n_peak       fragments whose central maximum is the profile maximum over the 32 bins (apex inside the core)
  log_area     log10 summed central fragment intensity
  frac_core    central share of the summed 32-bin fragment intensity (low = signal is a slope or a neighbour peak)
  corr_med     median pairwise Pearson correlation of the fragment profiles over the 32 bins (fragments with signal)
  top1_share   share of the most intense fragment in the central area (high = one fragment carries the window)
  lib_cos      cosine sqrt(central area) vs sqrt(library intensity)
  rank_agree   Spearman of central fragment areas against library intensities
  ppm_abs      intensity-weighted |mean ppm - 0| over central bins, median over fragments with signal (export targets are
               already calibrated, so 0 is the expected value)
  ms1_m0       log10 central M+0 area;  ms1_corr  Pearson of M+0 profile with summed fragments; iso_ratio M+1/M+0 area
  bins_nz      central bins with any fragment signal (of 8; low = sparse sampling)
Also per window: missing fragments (no signal anywhere) and, for 12-fragment arrays, n_sig12 and extra (non-original)
fragments with central signal."""
import sys, numpy as np, pandas as pd
O = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2'
C0, C1 = 12, 20


def pearson_rows(a):
  a = a - a.mean(1, keepdims=True); n = np.sqrt((a * a).sum(1)); ok = n > 0
  a = a[ok] / n[ok, None]; c = a @ a.T; iu = np.triu_indices(len(a), 1)
  return np.median(c[iu]) if len(iu[0]) else np.nan


def feats(rows, A):
  rows = np.sort(rows); xic = np.asarray(A['xic'][rows], np.float64); ppm = np.asarray(A['ppm'][rows], np.float64)
  ms1 = np.asarray(A['ms1'][rows], np.float64); lib = np.asarray(A['lib'][rows], np.float64)
  x12 = np.asarray(A['xic12'][rows], np.float64); o12 = np.asarray(A['orig12'][rows]) > 0
  core = xic[:, :, C0:C1]; area = core.sum(2); tot = xic.sum((1, 2)); s = xic.sum(1)
  out = pd.DataFrame(index=rows)
  out['n_sig'] = (area > 0).sum(1)
  out['n_missing'] = (xic.max(2) == 0).sum(1)
  out['n_peak'] = ((core.max(2) >= xic.max(2)) & (area > 0)).sum(1)
  out['log_area'] = np.log10(area.sum(1) + 1)
  out['frac_core'] = area.sum(1) / np.maximum(tot, 1e-12)
  out['corr_med'] = [pearson_rows(x[x.max(1) > 0]) for x in xic]
  out['top1_share'] = area.max(1) / np.maximum(area.sum(1), 1e-12)
  sa, sl = np.sqrt(area), np.sqrt(lib)
  out['lib_cos'] = (sa * sl).sum(1) / np.maximum(np.linalg.norm(sa, axis=1) * np.linalg.norm(sl, axis=1), 1e-12)
  out['rank_agree'] = [pd.Series(a).rank().corr(pd.Series(l).rank()) if (a > 0).sum() >= 3 else np.nan for a, l in zip(area, lib)]
  mp = np.where(core > 0, ppm[:, :, C0:C1] / np.maximum(core, 1e-12), np.nan)
  wm = np.nansum(mp * core, 2) / np.maximum(area, 1e-12)
  out['ppm_abs'] = [np.median(np.abs(w[a > 0])) if (a > 0).any() else np.nan for w, a in zip(wm, area)]
  m0 = ms1[:, 0]; out['ms1_m0'] = np.log10(m0[:, C0:C1].sum(1) + 1)
  out['ms1_corr'] = [np.corrcoef(a, b)[0, 1] if a.std() > 0 and b.std() > 0 else np.nan for a, b in zip(m0, s)]
  out['iso_ratio'] = ms1[:, 1, C0:C1].sum(1) / np.maximum(m0[:, C0:C1].sum(1), 1e-12)
  out['bins_nz'] = (core.sum(1) > 0).sum(1)
  a12 = x12[:, :, C0:C1].sum(2); out['n_sig12'] = (a12 > 0).sum(1); out['extra_sig'] = ((a12 > 0) & ~o12).sum(1)
  return out


for raw in (sys.argv[1:] or ['raw01', 'raw02']) if __name__ == '__main__' else []:
  A = {k: np.load(f'{O}/{raw}/{k}.npy', mmap_mode='r') for k in ('xic', 'ppm', 'ms1', 'lib', 'xic12', 'orig12')}
  pop = pd.read_csv(f'{O}/{raw}/population.tsv', sep='\t'); s = np.load(f'{O}/{raw}/scores_D11FWT.npy').astype(np.float64)
  w = pop.assign(s=s, row=np.arange(len(pop))).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id')
  w = w[w.s > -1e20]
  st = pd.read_csv(f'{O}/v3seq/misses_D11FWT_{raw}.tsv', sep='\t').set_index('precursor_id').stage
  dec = w[w.is_decoy == 1].sort_values('s', ascending=False)
  groups = {'E accepted': w[w.precursor_id.map(st) == 'E accepted'].row.values,
            'D below cut': w[w.precursor_id.map(st) == 'D right window, below cut'].row.values,
            'T1 top-1 % decoys': dec.row.values[:len(dec) // 100], 'T01 top-0.1 % decoys': dec.row.values[:len(dec) // 1000]}
  tabs = {k: feats(v, A) for k, v in groups.items()}
  allf = pd.concat({k: t for k, t in tabs.items()}, names=['group', 'row'])
  allf.to_csv(f'{O}/v4/dd_feats_{raw}.tsv', sep='\t')
  print(f'\n== {raw}: medians (n per group: ' + ', '.join(f'{k} {len(v)}' for k, v in groups.items()) + ')')
  print(allf.groupby(level='group').median().T.round(3).to_string())
