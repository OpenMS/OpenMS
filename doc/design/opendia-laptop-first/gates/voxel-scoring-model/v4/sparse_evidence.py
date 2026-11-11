"""Diagnostic (no training, no gate): do the stage-D misses (DIA-NN 1 % precursors whose right window scores below the
FDP_c 1 % cut; reference only) carry a 'sparse but strong' pattern (a few strongly co-eluting fragments plus co-eluting
MS1), and how often do decoy winners carry the same pattern? Arm D11FW. Per window, from xic12 / lib12 / ms1feat:
  n_strong   fragments with central signal and correlation >= 0.9 to the sum of the other fragments
  top3_corr  mean of the 3 highest fragment correlations
  wtd_corr   predicted-intensity-weighted mean fragment correlation (absent fragments count as 0)
  ms1_corr   M+0 vs fragment-sum correlation; ms1_iso isotope-envelope cosine
Pattern S: top3_corr >= 0.9 and ms1_corr >= 0.8 and ms1_iso >= 0.9."""
import numpy as np, pandas as pd
O = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2'
ARM = 'D11FW'; C0, C1 = 12, 20
ref_all = pd.read_parquet('/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-diann-gap-20261108/ref_exact.parquet')

def win_feats(raw, rows):
    X = np.load(f'{O}/{raw}/xic12.npy', mmap_mode='r'); LI = np.load(f'{O}/{raw}/lib12.npy', mmap_mode='r'); M = np.load(f'{O}/{raw}/ms1feat.npy', mmap_mode='r')
    o = np.argsort(rows); inv = np.argsort(o); sr = rows[o]
    x = np.asarray(X[sr], np.float64)[inv]; li = np.asarray(LI[sr], np.float64)[inv]; m = np.asarray(M[sr], np.float64)[inv]
    present = li > 0; tot = x.sum(1, keepdims=True); oth = tot - x
    xc = x - x.mean(2, keepdims=True); oc = oth - oth.mean(2, keepdims=True)
    r = np.where(present, (xc * oc).sum(2) / np.sqrt((xc * xc).sum(2) * (oc * oc).sum(2) + 1e-12), 0)
    sig = x[:, :, C0:C1].max(2) > 0; r = np.where(sig, r, 0)
    top3 = -np.sort(-r, 1)[:, :3].mean(1)
    return pd.DataFrame({'n_strong': (r >= 0.9).sum(1), 'top3_corr': top3, 'wtd_corr': (r * li).sum(1) / np.maximum(li.sum(1), 1e-12),
                         'ms1_corr': m[:, 1], 'ms1_iso': m[:, 2], 'ms1_present': m[:, 0] > 0})

for raw, tag in (('raw01', '_01'), ('raw02', '_02')):
    pop = pd.read_csv(f'{O}/{raw}/population.tsv', sep='\t'); s = np.load(f'{O}/{raw}/scores_{ARM}.npy')
    win = pop.assign(s=s, row=np.arange(len(pop))).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').copy()
    for f_, g in win[win.is_decoy == 1].groupby('fold'):
        med = g.s.median(); mad = 1.4826 * (g.s - med).abs().median(); win.loc[win.fold == f_, 's'] = (win.loc[win.fold == f_, 's'] - med) / mad
    st = pd.read_csv(f'{O}/v3seq/misses_{ARM}_{raw}.tsv', sep='\t').set_index('precursor_id').stage
    tw = win[win.is_decoy == 0].copy(); tw['stage'] = tw.precursor_id.map(st)
    dw = win[win.is_decoy == 1].sort_values('s', ascending=False)
    E = tw[tw.stage == 'E accepted']; cut = E.s.min()
    D = tw[tw.stage == 'D right window, below cut']  # the winner is the right window for stage D
    groups = {'E accepted': E, 'D below cut': D, 'targets not in DIA-NN, below cut': tw[tw.stage.isna() & (tw.s < cut)].sample(20000, random_state=0),
              'decoy winners, top 1 %': dw.head(len(dw) // 100), 'decoy winners, all': dw.sample(20000, random_state=0)}
    print(f'\n== {raw} {ARM}: cut (normalised) {cut:.1f}')
    print(f"{'group':36s} {'n':>6s} {'n_strong':>8s} {'top3':>6s} {'wtd':>6s} {'ms1 pres':>8s} {'ms1 corr':>8s} {'ms1 iso':>7s} {'S rate':>7s} {'S & ts<=4':>9s}")
    H = np.load(f'{O}/{raw}/heuristics.npy', mmap_mode='r')
    for name, g in groups.items():
        f = win_feats(raw, g.row.values); ts = np.asarray(H[np.sort(g.row.values), 0])[np.argsort(np.argsort(g.row.values))]
        S = (f.top3_corr >= 0.9) & (f.ms1_corr >= 0.8) & (f.ms1_iso >= 0.9)
        print(f"{name:36s} {len(g):6d} {f.n_strong.median():8.0f} {f.top3_corr.median():6.2f} {f.wtd_corr.median():6.2f} {f.ms1_present.mean():8.2f} "
              f"{f.ms1_corr[f.ms1_present].median():8.2f} {f.ms1_iso[f.ms1_present].median():7.2f} {S.mean():7.3f} {(S & (ts <= 4)).mean():9.3f}")
    # D precursors with the sparse-strong pattern: how far below the cut, and how many decoys above that score show S
    fD = win_feats(raw, D.row.values); SD = ((fD.top3_corr >= 0.9) & (fD.ms1_corr >= 0.8) & (fD.ms1_iso >= 0.9)).values
    fdec = win_feats(raw, dw.row.values); Sdec = ((fdec.top3_corr >= 0.9) & (fdec.ms1_corr >= 0.8) & (fdec.ms1_iso >= 0.9)).values
    print(f'   D with pattern S: {SD.sum()} of {len(D)}; their score/cut median {np.median(D.s.values[SD] / cut):.2f}. '
          f'Decoy winners with S: {Sdec.sum()} of {len(dw)}; target winners with S (all): ', end='')
    ft = win_feats(raw, tw.row.values); St = ((ft.top3_corr >= 0.9) & (ft.ms1_corr >= 0.8) & (ft.ms1_iso >= 0.9)).values
    print(f'{St.sum()} of {len(tw)} -> S-only target/decoy ratio {St.sum() / max(Sdec.sum(), 1):.2f} (all winners {len(tw) / len(dw):.2f})')
    # Within S: targets minus decoys gives a decoy-count estimate of true S precursors below the cut
    below = tw.s.values < cut; dbel = dw.s.values < cut
    print(f'   below the cut: S targets {int((St & below).sum())}, S decoys {int((Sdec & dbel).sum())} -> excess {int((St & below).sum() - (Sdec & dbel).sum())}')
