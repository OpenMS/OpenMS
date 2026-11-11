"""Graded follow-up to sparse_evidence.py: below-cut precursor winners (D11FW), bucketed by k = fragments with
correlation >= 0.8 to the other fragments and by MS1 co-elution (ms1_corr >= 0.6). Per bucket: target and decoy counts,
excess T - D (decoy-count estimate of true precursors), and DIA-NN stage-D precursors in it (reference only)."""
import numpy as np, pandas as pd, sys
sys.argv = ['x']; exec(open(__file__.replace('sparse_grid', 'sparse_evidence')).read().split('for raw, tag in')[0])
for raw in ('raw01', 'raw02'):
    pop = pd.read_csv(f'{O}/{raw}/population.tsv', sep='\t'); s = np.load(f'{O}/{raw}/scores_{ARM}.npy')
    win = pop.assign(s=s, row=np.arange(len(pop))).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').copy()
    for f_, g in win[win.is_decoy == 1].groupby('fold'):
        med = g.s.median(); mad = 1.4826 * (g.s - med).abs().median(); win.loc[win.fold == f_, 's'] = (win.loc[win.fold == f_, 's'] - med) / mad
    st = pd.read_csv(f'{O}/v3seq/misses_{ARM}_{raw}.tsv', sep='\t').set_index('precursor_id').stage
    win['stage'] = win.precursor_id.map(st); cut = win[(win.is_decoy == 0) & (win.stage == 'E accepted')].s.min()
    b = win[win.s < cut].copy()
    X = np.load(f'{O}/{raw}/xic12.npy', mmap_mode='r'); LI = np.load(f'{O}/{raw}/lib12.npy', mmap_mode='r'); M = np.load(f'{O}/{raw}/ms1feat.npy', mmap_mode='r')
    rows = b.row.values; o = np.argsort(rows); sr = rows[o]; inv = np.argsort(o)
    k = np.zeros(len(rows), int)
    for a in range(0, len(sr), 50000):
        x = np.asarray(X[sr[a:a + 50000]], np.float64); li = np.asarray(LI[sr[a:a + 50000]])
        oth = x.sum(1, keepdims=True) - x; xc = x - x.mean(2, keepdims=True); oc = oth - oth.mean(2, keepdims=True)
        r = (xc * oc).sum(2) / np.sqrt((xc * xc).sum(2) * (oc * oc).sum(2) + 1e-12)
        r = np.where((li > 0) & (x[:, :, 12:20].max(2) > 0), r, 0); k[a:a + 50000] = (r >= 0.8).sum(1)
    b['k'] = np.minimum(k[inv], 5); b['ms1'] = np.asarray(M[sr, 1])[inv] >= 0.6
    b['isD'] = (b.is_decoy == 0) & (b.stage == 'D right window, below cut')
    g = b.groupby(['k', 'ms1']).agg(T=('is_decoy', lambda v: int((v == 0).sum())), Dc=('is_decoy', 'sum'), dianD=('isD', 'sum'))
    g['excess'] = g['T'] - g['Dc']; g['T/D'] = (g['T'] / g['Dc'].clip(lower=1)).round(2)
    print(f'\n== {raw} {ARM} below the 1 % cut ({len(b)} winners); k = fragments with r >= 0.8 (5 = 5+), ms1 = M+0 co-elution r >= 0.6')
    print(g.to_string())
