"""D11 misses, part e (descriptive): score saturation (floor) and window choice. For DIA-NN reference precursors: does the
D11 winner match the DIA-NN RT, against the BASE winner; and are the competing windows in the score floor."""
import numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; O = B + '/astral-voxel-pilot-20261110/v2'
BASEF = {'raw01': B + '/astral-design-raw01-20261107/stage2-byion-design/byion/features.tsv', 'raw02': B + '/astral-raw02-byion-design-20261109/byion/features.tsv'}
for raw in ('raw01', 'raw02'):
    pop = pd.read_csv(f'{O}/{raw}/population.tsv', sep='\t'); s = np.load(f'{O}/{raw}/scores_D11.npy')
    base = pd.read_csv(BASEF[raw], sep='\t', usecols=['sample_id', 'stage2_score']).set_index('sample_id').stage2_score.loc[pop.sample_id].values
    for k in range(3):
        m = pop.fold.values == k; q = np.quantile(s[m], [0.01, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999])
        print(f'{raw} fold {k} D11 window score quantiles 1/25/50/75/90/99/99.9 %: ' + ' '.join(f'{v:.1f}' for v in q),
              f'| share within 1 of the fold minimum: {(s[m] < s[m].min() + 1).mean():.2f}, distinct values: {len(np.unique(np.round(s[m], 3)))}')
    st = pd.read_csv(f'{O}/v3seq/misses_D11_{raw}.tsv', sep='\t')
    w = pop.assign(s=s, b=base).merge(st[['precursor_id', 'rt_s', 'stage']], on='precursor_id'); w['right'] = (w.rt_apex - w.rt_s).abs() <= 10
    has = w.groupby('precursor_id').right.any(); w = w[w.precursor_id.isin(has[has].index)]
    dw = w.sort_values('s', ascending=False).drop_duplicates('precursor_id').set_index('precursor_id').right
    bw = w.sort_values('b', ascending=False).drop_duplicates('precursor_id').set_index('precursor_id').right.loc[dw.index]
    print(f'{raw}: precursors with a right window {len(dw)}; winner is right: D11 {int(dw.sum())}, BASE {int(bw.sum())}; '
          f'BASE right & D11 wrong {int((bw & ~dw).sum())}, D11 right & BASE wrong {int((dw & ~bw).sum())}')
    C = w[w.stage == 'C wrong window won']
    fl = {k: s[pop.fold.values == k].min() for k in range(3)}
    C = C.assign(floor=C.fold.map(fl)); C['at_floor'] = C.s < C.floor + 5
    g = C.groupby('precursor_id')
    print(f'{raw} C: right window at floor (within 5 of the fold min) {int(C[C.right].groupby("precursor_id").at_floor.all().sum())} of {g.ngroups}; '
          f'wrong winner at floor {int(C.sort_values("s", ascending=False).drop_duplicates("precursor_id").at_floor.sum())}')
