"""Diagnostic: fragment mass accuracy as sparse-evidence. From xic.npy / ppm.npy (6 original fragments; ppm = sum
intensity x (ppm - calibrated offset) per bin over +- 17.5 ppm). Per window, central bins (apex +- 4.7 s), fragments with signal:
  ppm_abs_med   median |intensity-weighted mean ppm| over fragments (random matches spread over +- 17.5)
  n_tight       fragments with |mean ppm| <= 3
  tight_share   intensity share of fragments with |mean ppm| <= 3
AUC of stage-D right windows (DIA-NN reference) against the top 1 % decoy winners and against below-cut decoy winners."""
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
O = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2'; ARM = 'D11FW'
def feats(raw, rows):
    X = np.load(f'{O}/{raw}/xic.npy', mmap_mode='r'); P = np.load(f'{O}/{raw}/ppm.npy', mmap_mode='r')
    o = np.argsort(rows); inv = np.argsort(o); x = np.asarray(X[rows[o]], np.float64)[inv][:, :, 12:20]; p = np.asarray(P[rows[o]], np.float64)[inv][:, :, 12:20]
    a = x.sum(2); m = np.where(a > 0, p.sum(2) / np.maximum(a, 1e-12), np.nan); sig = a > 0
    tight = sig & (np.abs(np.nan_to_num(m, nan=99)) <= 3)
    return pd.DataFrame({'ppm_abs_med': np.nan_to_num(np.nanmedian(np.abs(m), 1), nan=17.5), 'n_tight': tight.sum(1),
                         'tight_share': (a * tight).sum(1) / np.maximum(a.sum(1), 1e-12)})
for raw in ('raw01', 'raw02'):
    pop = pd.read_csv(f'{O}/{raw}/population.tsv', sep='\t'); s = np.load(f'{O}/{raw}/scores_{ARM}.npy')
    win = pop.assign(s=s, row=np.arange(len(pop))).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').copy()
    for f_, g in win[win.is_decoy == 1].groupby('fold'):
        med = g.s.median(); mad = 1.4826 * (g.s - med).abs().median(); win.loc[win.fold == f_, 's'] = (win.loc[win.fold == f_, 's'] - med) / mad
    st = pd.read_csv(f'{O}/v3seq/misses_{ARM}_{raw}.tsv', sep='\t').set_index('precursor_id').stage; win['stage'] = win.precursor_id.map(st)
    cut = win[(win.is_decoy == 0) & (win.stage == 'E accepted')].s.min()
    D = win[(win.is_decoy == 0) & (win.stage == 'D right window, below cut')]; E = win[(win.is_decoy == 0) & (win.stage == 'E accepted')]
    dw = win[win.is_decoy == 1].sort_values('s', ascending=False); top = dw.head(len(dw) // 100); bel = dw[dw.s < cut].sample(20000, random_state=0)
    fD, fE, fT, fB = (feats(raw, g.row.values) for g in (D, E, top, bel))
    print(f'\n== {raw} {ARM}: medians E / D / top-1% decoys | AUC D vs top-1% decoys (higher = D looks more real), D vs below-cut decoys')
    for c, sign in (('ppm_abs_med', -1), ('n_tight', 1), ('tight_share', 1)):
        a1 = mannwhitneyu(sign * fD[c], sign * fT[c]).statistic / (len(fD) * len(fT)); a2 = mannwhitneyu(sign * fD[c], sign * fB[c]).statistic / (len(fD) * len(fB))
        print(f'{c:12s} {fE[c].median():6.2f} / {fD[c].median():6.2f} / {fT[c].median():6.2f} | {a1:.3f}  {a2:.3f}')
