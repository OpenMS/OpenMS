"""Post-gate design question (diagnostic, not a gate): how many training IDs does the RT correction need?
Train on the top-N OpenDIA IDs by score (a proxy for calibration anchors: 140 on D, 441 on H), cross-fitted as before."""
import os, sys, importlib.util
sys.argv = ['x']
src = open('/nfs/research/juan/singj/opendia/imfeat-20261120/precheck_rtim.py').read()
head = src.split('out = {}')[0]
exec(head)
for run in ('D', 'H'):
    rd = f'{C}/entlib-mt-20261120/runs/{run}_mt_s0'
    cal = json.load(open(glob.glob(f'{rd}/*.calibration.json')[0]))
    P = pd.read_csv([f for f in glob.glob(f'{rd}/*.precursors.tsv') if not f.endswith('/precursors.tsv')][0], sep='\t', usecols=['precursor_id', 'is_decoy', 'rt_apex', 'score', 'q_value'])
    P = P[P.score > -1e29].sort_values('score', ascending=False).drop_duplicates('precursor_id').set_index('precursor_id')
    ours = P[(P.is_decoy == 0) & (P.q_value <= 0.01)]
    d = pd.read_parquet(f'{C}/entlib-diann5-20261119/{run}/report.parquet', columns=['Modified.Sequence', 'Precursor.Charge', 'Q.Value', 'Decoy', 'RT'])
    d = d[(d['Q.Value'] <= 0.01) & (d.Decoy == 0)].copy(); d['pid'] = sk.reindex(d['Modified.Sequence'] + '/' + d['Precursor.Charge'].astype(str)).values
    d = d.dropna(subset=['pid']); d['pid'] = d.pid.astype(np.int64); d = d[tg.e.reindex(d.pid).values == 0].drop_duplicates('pid').set_index('pid')
    ids = sorted(set(ours.index) | set(d.index)); T = tg.loc[ids].copy()
    T['ours'] = T.index.isin(ours.index); T['dnn'] = T.index.isin(d.index)
    T['rt_obs'] = ours.rt_apex.reindex(T.index).values; T['rt_dnn'] = d.RT.reindex(T.index).values * 60; T['rt_map'] = apply_map(cal, T.lrt.values)
    T['rank'] = ours.score.rank(ascending=False).reindex(T.index).values
    X = featurise(T); folds = fold_of(T.ms.values); ev = T.dnn.values; mi = ev & ~T.ours.values
    base = T.rt_dnn.values - T.rt_map.values
    print(f'\n{run}: map p95 all {q(base[ev], .95):.1f}, missed {q(base[mi], .95):.1f} s')
    for n in (140, 441, 1000, 3000, 10000, int(T.ours.sum())):
        trm = T.ours.values & (T['rank'].values <= n)
        y = np.where(trm, T.rt_obs.values - T.rt_map.values, 0.0)
        pr, al = crossfit(X, y, trm, folds); r = base - pr
        print(f'  top {n:6d} IDs (train {trm.sum():6d}): p95 all {q(r[ev], .95):6.1f}  missed {q(r[mi], .95):6.1f} s  ({(q(r[mi], .95) / q(base[mi], .95) - 1) * 100:+.1f} %), alpha {np.median(al):.0f}')
