"""PRECHECK_NOCAND part 2: picker misses inside the search window, from the search's own trace matrix.
Mini library of sampled DIA-NN IDs (600 picker misses, 150 cap misses, 300 found), OpenDIA dev tool search with the run's
calibration (same picker as the workflow search), -traces_out. For each precursor: consensus trace (sum of fragments)
and support near DIA-NN's RT vs the precursor's best picked peak. Diagnostic only."""
import numpy as np, pandas as pd
O = '/home/singjc/.jcode/scratch/dia2/gap7/nocand'
S = pd.read_csv(f'{O}/sample.tsv', sep='\t', index_col=0)
P = pd.read_csv(f'{O}/peaks.tsv', sep='\t')
T = pd.read_csv(f'{O}/traces.tsv', sep='\t')
F = pd.read_parquet('/home/singjc/.jcode/scratch/dia2/gap7/D36_feat.parquet', columns=['precursor_id', 'is_decoy', 'rt_apex'])
F = F[F.precursor_id.isin(S.index) & (F.is_decoy == 0)]
# Reproduction: does the dev-tool search return the same apexes as the workflow?
a = P.groupby('precursor_id').rt_apex.apply(lambda x: tuple(sorted(np.round(x, 3))))
b = F.groupby('precursor_id').rt_apex.apply(lambda x: tuple(sorted(np.round(x, 3))))
both = a.index.intersection(b.index)
print(f'apex sets identical to the workflow for {np.mean([a[i] == b[i] for i in both]):.3f} of {len(both)} precursors '
      f'(workflow-only {len(b.index.difference(a.index))}, dev-only {len(a.index.difference(b.index))})')
P['drt'] = np.abs(P.rt_apex - S.rt_s.reindex(P.precursor_id).values)
S['right_dev'] = P.groupby('precursor_id').drt.min().reindex(S.index) <= 6
print('dev search right-peak share by category:', S.groupby('cat').right_dev.mean().round(3).to_dict())

rows = []
for pid, t in T.groupby('precursor_id', sort=False):
    m = t.pivot_table(index='rt', columns='slot', values='intensity', aggfunc='max').fillna(0.0)
    rt = m.index.values; X = m.values; cons = X.sum(1)
    sm = np.convolve(np.r_[cons[:1], cons, cons[-1:]], [0.25, 0.5, 0.25], 'valid')
    tr = S.rt_s[pid]
    near = np.abs(rt - tr) <= 6
    if not near.any(): rows.append((pid, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 0)); continue
    i = np.where(near)[0][np.argmax(sm[near])]
    # local maximum of smoothed consensus at DIA-NN RT? support = fragments with signal at apex +- 1
    lm = (i > 0 and sm[i] >= sm[i - 1]) and (i < len(sm) - 1 and sm[i] >= sm[i + 1])
    sup = int((X[max(i - 1, 0):i + 2] > 0).any(0).sum())
    nz_near = (X[near] > 0).mean()
    # best picked peak height (consensus at its apex)
    pk = P[P.precursor_id == pid]
    best_h = 0.0
    for ra in pk.rt_apex.values:
        j = np.argmin(np.abs(rt - ra)); best_h = max(best_h, cons[j])
    # max within the whole window
    rows.append((pid, cons[i], sm[i], best_h, cons.max(), float(lm), sup, X.shape[1], nz_near))
R = pd.DataFrame(rows, columns=['pid', 'h_dn', 'sm_dn', 'h_best_pick', 'h_max', 'local_max', 'support', 'frags', 'nz_near'][:len(rows[0])]).set_index('pid')
R = R.join(S[['cat', 'qq', 'n']])
R['rel_to_pick'] = R.h_dn / R.h_best_pick.replace(0, np.nan)
R['rel_to_max'] = R.h_dn / R.h_max.replace(0, np.nan)
pd.set_option('display.width', 200)
print('\nPer category (median unless noted):')
print(R.groupby('cat').agg(n=('h_dn', 'size'), h_dn=('h_dn', 'median'), rel_to_pick=('rel_to_pick', 'median'), rel_to_max=('rel_to_max', 'median'),
                           zero_at_dn=('h_dn', lambda x: np.mean(x == 0)), local_max=('local_max', 'mean'), support=('support', 'median'),
                           support_le2=('support', lambda x: np.mean(x <= 2)), nz_near=('nz_near', 'median')).round(3).to_string())
M = R[R.cat == 'picker']
print('\npicker misses, support at DIA-NN RT (fragments with signal at apex +-1 scan):', M.support.value_counts().sort_index().to_dict())
print('picker misses, consensus at DIA-NN RT relative to window max: share <0.1', round(np.mean(M.rel_to_max < 0.1), 3),
      ' 0.1-0.5', round(np.mean((M.rel_to_max >= 0.1) & (M.rel_to_max < 0.5)), 3), ' >=0.5', round(np.mean(M.rel_to_max >= 0.5), 3))
print('picker misses that are a smoothed local max at DIA-NN RT but not picked:', round(M.local_max.mean(), 3))
R.to_csv(f'{O}/traces_summary.tsv', sep='\t')

# Merge test: is DIA-NN's RT inside the [rt_start, rt_end] of a picked peak whose apex is > 6 s away?
P['inside'] = (S.rt_s.reindex(P.precursor_id).values >= P.rt_start) & (S.rt_s.reindex(P.precursor_id).values <= P.rt_end)
ins = P.groupby('precursor_id').inside.any().reindex(S.index).fillna(False)
S['inside_pick'] = ins
print('\nDIA-NN RT inside a picked peak region (apex elsewhere):', S.groupby('cat').inside_pick.mean().round(3).to_dict())
Q = P[P.inside].copy(); Q['cat'] = S.cat.reindex(Q.precursor_id).values; Q = Q[Q.cat != 'found']
Q['half_span'] = (Q.rt_end - Q.rt_start) / 2
print('merged-into peaks (misses): apex distance to DIA-NN RT p50 %.1f s, peak width p50 %.1f s (found peaks p50 %.1f s)' % (
    Q.drt.median(), (Q.rt_end - Q.rt_start).median(), (P[S.cat.reindex(P.precursor_id).values == 'found'].pipe(lambda x: x.rt_end - x.rt_start)).median()))

# Valley test against the nearest picked apex, on the [1 2 1]-smoothed consensus (picker: separate only if valley <= 0.5 x lower max).
vr = {}
for pid, t in T[T.precursor_id.isin(S.index[S.cat != 'found'])].groupby('precursor_id', sort=False):
    m = t.pivot_table(index='rt', columns='slot', values='intensity', aggfunc='max').fillna(0.0)
    rt = m.index.values; cons = m.values.sum(1)
    sm = np.convolve(np.r_[cons[:1], cons, cons[-1:]], [0.25, 0.5, 0.25], 'valid')
    near = np.abs(rt - S.rt_s[pid]) <= 6
    if not near.any(): continue
    i = np.where(near)[0][np.argmax(sm[near])]
    pk = P[P.precursor_id == pid]
    if not len(pk): vr[pid] = (np.nan, np.nan, np.inf); continue
    j = np.argmin(np.abs(rt - pk.rt_apex.values[np.argmin(np.abs(pk.rt_apex.values - rt[i]))]))
    lo_, hi_ = sorted((i, j)); valley = sm[lo_:hi_ + 1].min()
    vr[pid] = (valley / max(min(sm[i], sm[j]), 1e-9), sm[i] / max(sm[j], 1e-9), abs(hi_ - lo_))
V = pd.DataFrame.from_dict(vr, orient='index', columns=['valley_ratio', 'height_ratio', 'scans_apart'])
V['cat'] = S.cat.reindex(V.index)
for cat, v in V.groupby('cat'):
    print(f'{cat}: valley/lower max > 0.5 (merged by valley rule) {np.mean(v.valley_ratio > 0.5):.3f}; '
          f'scans to nearest pick p50 {v.scans_apart.median():.0f}; DIA-NN max / nearest pick max p50 {v.height_ratio.median():.2f}')
