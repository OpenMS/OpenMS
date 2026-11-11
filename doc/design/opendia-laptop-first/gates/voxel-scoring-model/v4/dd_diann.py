"""Stage-D deep dive, step 3: does DIA-NN measure the same signal as our winning window? (descriptive; DIA-NN report is
reference only). Join our D11FWT winners (E accepted and D below cut) with DIA-NN's report row for the same precursor/run:
  rt_diff        our apex - DIA-NN RT (s), and DIA-NN FWHM / RT.Start..Stop
  qty            DIA-NN Precursor.Quantity (log10) against our central fragment area (log10, dd_feats)
  ms1            DIA-NN Ms1.Area (log10) against our central M+0 area
  evidence       DIA-NN Evidence, Mass.Evidence, Ms1.Profile.Corr, Empirical.Quality, Quantity.Quality
  best_fr        DIA-NN Best.Fr.Mz: is it one of our 6 library fragments (within 10 ppm)? Which rank in our library?
  pep            DIA-NN PEP, Q.Value
Then: among D, which of these explains the low score (by comparison with E at matched DIA-NN quantity)."""
import sqlite3, numpy as np, pandas as pd, pyarrow.parquet as pq
R = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; O = R + '/astral-voxel-pilot-20261110/v2'
DN = '/home/singjc/Documents/datasets/proteobench-astral-entrapment/proteobench_diann_v231_res_20260717_123332/input_file.parquet'
cols = ['Run', 'Stripped.Sequence', 'Modified.Sequence', 'Precursor.Charge', 'RT', 'RT.Start', 'RT.Stop', 'FWHM', 'Precursor.Quantity', 'Ms1.Area',
        'Ms1.Apex.Mz.Delta', 'Evidence', 'Mass.Evidence', 'Ms1.Profile.Corr', 'Empirical.Quality', 'Quantity.Quality', 'Q.Value', 'PEP',
        'Best.Fr.Mz', 'Best.Fr.Mz.Delta', 'Decoy']
dn = pq.read_table(DN, columns=cols).to_pandas(); dn = dn[dn.Decoy == 0]
dn['key'] = dn['Stripped.Sequence'] + '/' + dn['Precursor.Charge'].astype(str)
c = sqlite3.connect('file:' + R + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
tr = pd.read_sql('select m.PRECURSOR_ID pid, t.PRODUCT_MZ fmz, t.LIBRARY_INTENSITY li from TRANSITION t join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID', c)
tr['rank'] = tr.groupby('pid').li.rank(ascending=False, method='first')
fr = {p: g[['fmz', 'rank']].values for p, g in tr.groupby('pid')}
for raw in ('raw01', 'raw02'):
  pop = pd.read_csv(f'{O}/{raw}/population.tsv', sep='\t'); s = np.load(f'{O}/{raw}/scores_D11FWT.npy').astype(np.float64)
  w = pop.assign(s=s, row=np.arange(len(pop))).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id')
  st = pd.read_csv(f'{O}/v3seq/misses_D11FWT_{raw}.tsv', sep='\t')
  st = st[st.stage.isin(['E accepted', 'D right window, below cut'])].merge(w[['precursor_id', 'rt_apex', 'row', 's']], on='precursor_id')
  f = pd.read_csv(f'{O}/v4/dd_feats_{raw}.tsv', sep='\t'); f = f[f.group.isin(['E accepted', 'D below cut'])].drop_duplicates('row').set_index('row')
  st = st.join(f[['log_area', 'ms1_m0', 'lib_cos', 'corr_med', 'n_peak', 'frac_core']], on='row')
  d = dn[dn.Run.str.endswith('_' + raw[-2:])].drop_duplicates('key').set_index('key')
  st = st.join(d, on='key', rsuffix='_dn')
  st['rt_diff'] = st.rt_apex - st.RT * 60; st['fwhm_s'] = st.FWHM * 60
  st['lq'] = np.log10(st['Precursor.Quantity']); st['lms1'] = np.log10(st['Ms1.Area'].clip(lower=1))
  def best_rank(row):
    a = fr.get(row.precursor_id)
    if a is None or not np.isfinite(row['Best.Fr.Mz']): return np.nan
    d_ = np.abs(a[:, 0] - row['Best.Fr.Mz']) / row['Best.Fr.Mz'] * 1e6; i = d_.argmin()
    return a[i, 1] if d_[i] <= 10 else 0  # 0 = DIA-NN's best fragment is not one of our 6
  st['best_fr_rank'] = st.apply(best_rank, axis=1)
  st['grp'] = np.where(st.stage == 'E accepted', 'E', 'D')
  show = ['rt_diff', 'fwhm_s', 'lq', 'log_area', 'lms1', 'ms1_m0', 'Evidence', 'Mass.Evidence', 'Ms1.Profile.Corr', 'Empirical.Quality',
          'Quantity.Quality', 'PEP', 'Q.Value', 'lib_cos', 'corr_med']
  print(f'\n== {raw}: medians, E (accepted) vs D (below cut); our window joined with the DIA-NN report row')
  print(st.groupby('grp')[show].median().T.round(3).to_string())
  print(f'{raw}: |rt_diff| quantiles (s) by group'); print(st.groupby('grp').rt_diff.apply(lambda v: v.abs().quantile([0.5, 0.9]).round(1).tolist()).to_string())
  print(f'{raw}: DIA-NN best fragment among our 6 library fragments (0 = not in our library; 1 = our top library fragment)')
  print(pd.crosstab(st.grp, st.best_fr_rank, normalize='index').round(3).to_string())
  # Matched on DIA-NN quantity: within quantity deciles, the D share and how log_area compares (is our signal lower than
  # DIA-NN's quantity implies?)
  st['qbin'] = pd.qcut(st.lq, 10, labels=False, duplicates='drop')
  g = st.groupby('qbin').agg(n=('grp', 'size'), D_share=('grp', lambda v: (v == 'D').mean()), lq=('lq', 'median'))
  for k in ('log_area', 'lib_cos', 'corr_med', 'Evidence'):
    g[f'{k}_E'] = st[st.grp == 'E'].groupby('qbin')[k].median(); g[f'{k}_D'] = st[st.grp == 'D'].groupby('qbin')[k].median()
  print(f'{raw}: by DIA-NN quantity decile'); print(g.round(3).to_string())
  st.to_csv(f'{O}/v4/dd_diann_{raw}.tsv', sep='\t', index=False)
