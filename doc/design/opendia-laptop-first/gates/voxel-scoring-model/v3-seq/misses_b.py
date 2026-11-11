"""D11 misses, part b (descriptive): (1) wrong-window group C: all 11 features + score of the wrong winner against the right
window, and BASE's window choice for the same precursors; (2) group D against top decoy winners: per-feature AUC; (3) DIA-NN
q of each stage; (4) RT error of the right window against the calibration sigma."""
import sqlite3, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; O = B + '/astral-voxel-pilot-20261110/v2'
FEATS = ['transition_support', 'log_matched_scans', 'log_peak_width', 'abs_library_rt_error', 'boundary_log_fwhm', 'boundary_log_width',
         'boundary_coelution', 'boundary_library_correlation', 'boundary_apex_agreement', 'byion_b', 'byion_y']
BASEF = {'raw01': B + '/astral-design-raw01-20261107/stage2-byion-design/byion/features.tsv', 'raw02': B + '/astral-raw02-byion-design-20261109/byion/features.tsv'}
dq = pd.read_parquet(B + '/astral-diann-gap-20261108/diann_runs_1pct.parquet')
for raw, tag in (('raw01', '_01'), ('raw02', '_02')):
    pop = pd.read_csv(f'{O}/{raw}/population.tsv', sep='\t'); H = pd.DataFrame(np.load(f'{O}/{raw}/heuristics.npy'), columns=FEATS)
    s = np.load(f'{O}/{raw}/scores_D11.npy')
    base = pd.read_csv(BASEF[raw], sep='\t', usecols=['sample_id', 'stage2_score']).set_index('sample_id').stage2_score.loc[pop.sample_id].values
    st = pd.read_csv(f'{O}/v3seq/misses_D11_{raw}.tsv', sep='\t')
    q = dq[dq.Run.str.endswith(tag)].drop_duplicates('key').set_index('key')['Q.Value']; st['dq'] = st.key.map(q)
    print(f'\n=== {raw}: DIA-NN q by stage (median, share q<=1e-3)')
    print(st.groupby('stage').dq.agg(['count', 'median', lambda x: (x <= 1e-3).mean()]).round(5).to_string())
    w = pop.assign(s=s, b=base, row=np.arange(len(pop))).merge(st[['precursor_id', 'rt_s', 'stage']], on='precursor_id')
    w['right'] = (w.rt_apex - w.rt_s).abs() <= 10
    w['d11_rank'] = w.groupby('precursor_id').s.rank(ascending=False, method='first')
    w['base_rank'] = w.groupby('precursor_id').b.rank(ascending=False, method='first')
    C = w[w.stage == 'C wrong window won']
    cr = C[C.right].sort_values('d11_rank').drop_duplicates('precursor_id').set_index('precursor_id')
    cw = C[C.d11_rank == 1].set_index('precursor_id').loc[cr.index]
    print(f'C: n {len(cr)}; BASE ranks a right window first for {int((C[C.right].groupby("precursor_id").base_rank.min() == 1).sum())}')
    print(f'C: windows per precursor median {C.groupby("precursor_id").size().median():.0f}; right window D11 rank median {cr.d11_rank.median():.0f}')
    print(f'C: D11 score wrong winner minus right window: median {(cw.s.values - cr.s.values).__array__().mean():.2f} (mean), share with z-sum higher in wrong:')
    hw, hr = H.iloc[cw.row.values].reset_index(drop=True), H.iloc[cr.row.values].reset_index(drop=True)
    print('   feature                       wrong-winner median   right-window median   share wrong > right')
    for f in FEATS: print(f'   {f:30s} {hw[f].median():12.3f} {hr[f].median():20.3f} {(hw[f] > hr[f]).mean():18.2f}')
    # D vs top decoy winners (window rows): per-feature AUC (D > decoy)
    win = pop.assign(s=s, row=np.arange(len(pop))).sort_values('s', ascending=False).drop_duplicates('precursor_id')
    dw = win[win.is_decoy == 1]; top = dw.row.values[:len(dw) // 100]
    D = w[(w.stage == 'D right window, below cut') & (w.d11_rank == 1)].row.values
    E = w[(w.stage == 'E accepted') & (w.d11_rank == 1)].row.values
    print(f'D ({len(D)}) vs top-1% decoy winners ({len(top)}) and accepted E: per-feature AUC(D > top decoy), AUC(E > D)')
    for f in FEATS:
        a = mannwhitneyu(H[f].values[D], H[f].values[top]).statistic / (len(D) * len(top)); b = mannwhitneyu(H[f].values[E], H[f].values[D]).statistic / (len(E) * len(D))
        print(f'   {f:30s} {a:.3f}  {b:.3f}')
