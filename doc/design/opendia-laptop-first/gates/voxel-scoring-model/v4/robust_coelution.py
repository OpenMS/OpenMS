"""Diagnostic: interference-robust, normalised co-elution against the current sum-of-others measure (D11FW winners).
Per window: each fragment's XIC scaled to unit maximum; consensus = per-bin median over fragments with signal;
r_med = correlation of each fragment with the consensus. Features:
  sum_top3   mean top-3 correlation to the sum of the other fragments (current style)
  med_top3   mean top-3 correlation to the median consensus
  med_k      fragments with r_med >= 0.8
  expl_li    predicted intensity explained: sum(lib12 of fragments with r_med >= 0.8) / sum(lib12)  (missing weak fragments cost little)
AUC of stage-D right windows (DIA-NN reference) against the top 1 % decoy winners, and against all below-cut decoy winners."""
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
O = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2'; ARM = 'D11FW'
def corr(a, b):
    a = a - a.mean(-1, keepdims=True); b = b - b.mean(-1, keepdims=True)
    return (a * b).sum(-1) / np.sqrt((a * a).sum(-1) * (b * b).sum(-1) + 1e-12)
def feats(raw, rows):
    X = np.load(f'{O}/{raw}/xic12.npy', mmap_mode='r'); LI = np.load(f'{O}/{raw}/lib12.npy', mmap_mode='r')
    o = np.argsort(rows); inv = np.argsort(o); x = np.asarray(X[rows[o]], np.float64)[inv]; li = np.asarray(LI[rows[o]], np.float64)[inv]
    sig = (li > 0) & (x[:, :, 12:20].max(2) > 0)
    oth = x.sum(1, keepdims=True) - x; rs = np.where(sig, corr(x, oth), 0)
    xn = x / np.maximum(x.max(2, keepdims=True), 1e-12); xn = np.where(sig[:, :, None], xn, np.nan)
    cons = np.nan_to_num(np.nanmedian(xn, 1))[:, None, :]
    rm = np.where(sig, corr(np.nan_to_num(xn), np.broadcast_to(cons, x.shape)), 0)
    top = lambda r: -np.sort(-r, 1)[:, :3].mean(1)
    return pd.DataFrame({'sum_top3': top(rs), 'med_top3': top(rm), 'med_k': (rm >= 0.8).sum(1),
                         'expl_li': (li * (rm >= 0.8)).sum(1) / np.maximum(li.sum(1), 1e-12)})
for raw in ('raw01', 'raw02'):
    pop = pd.read_csv(f'{O}/{raw}/population.tsv', sep='\t'); s = np.load(f'{O}/{raw}/scores_{ARM}.npy')
    win = pop.assign(s=s, row=np.arange(len(pop))).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').copy()
    for f_, g in win[win.is_decoy == 1].groupby('fold'):
        med = g.s.median(); mad = 1.4826 * (g.s - med).abs().median(); win.loc[win.fold == f_, 's'] = (win.loc[win.fold == f_, 's'] - med) / mad
    st = pd.read_csv(f'{O}/v3seq/misses_{ARM}_{raw}.tsv', sep='\t').set_index('precursor_id').stage; win['stage'] = win.precursor_id.map(st)
    cut = win[(win.is_decoy == 0) & (win.stage == 'E accepted')].s.min()
    D = win[(win.is_decoy == 0) & (win.stage == 'D right window, below cut')]; E = win[(win.is_decoy == 0) & (win.stage == 'E accepted')]
    dw = win[win.is_decoy == 1].sort_values('s', ascending=False); top = dw.head(len(dw) // 100); bel = dw[dw.s < cut].sample(20000, random_state=0)
    tb = win[(win.is_decoy == 0) & (win.s < cut)].sample(20000, random_state=0)
    fD, fE, fT, fB, fU = (feats(raw, g.row.values) for g in (D, E, top, bel, tb))
    print(f'\n== {raw} {ARM}: medians E / D / top-1% decoys | AUC D vs top-1% decoys, D vs below-cut decoys | AUC below-cut targets vs below-cut decoys')
    for c in fD.columns:
        a1 = mannwhitneyu(fD[c], fT[c]).statistic / (len(fD) * len(fT)); a2 = mannwhitneyu(fD[c], fB[c]).statistic / (len(fD) * len(fB))
        a3 = mannwhitneyu(fU[c], fB[c]).statistic / (len(fU) * len(fB))
        print(f'{c:9s} {fE[c].median():6.2f} / {fD[c].median():6.2f} / {fT[c].median():6.2f} | {a1:.3f}  {a2:.3f} | {a3:.3f}')
