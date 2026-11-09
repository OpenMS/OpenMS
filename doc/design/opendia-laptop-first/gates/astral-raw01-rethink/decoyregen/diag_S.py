"""Why family S collapses on raw 02: semi-supervised iteration-0 positive count and score finiteness per family."""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/../gap2')
import oswlike_score as S
full = {k: getattr(S, k).copy() for k in ('X0', 'prec', 'dec', 'fold', 'pair', 'entr', 'right')}
fam = S.m.family.values
for F in ('O', 'P', 'S'):
    keep = (fam == 'T') | (fam == F)
    for k, v in full.items(): setattr(S, k, v[keep])
    X = S.X0; Z = np.nan_to_num((X - np.nanmean(X, 0)) / (np.nanstd(X, 0) + 1e-9)); main = Z[:, S.NAMES.index('xcorr_shape_w')]
    tr = np.where(S.fold != 0)[0]
    w = S.winners(tr, main); kept, q = S.tdc_q(w, main)
    print(F, 'iter0 pos at q<=0.15:', int(((~S.dec[kept]) & (q <= 0.15)).sum()), 'at 0.5:', int(((~S.dec[kept]) & (q <= 0.5)).sum()), 'min q', q.min().round(3), 'decoy winners', int(S.dec[w].sum()), 'of', len(w))
    # per-feature AUC target-winner vs decoy-winner (main ranking), largest separators
    from scipy.stats import rankdata
    tw = w[~S.dec[w]]; dw = w[S.dec[w]]; au = []
    for j, n in enumerate(S.NAMES):
        r = rankdata(np.r_[Z[tw, j], Z[dw, j]]); au.append(((r[:len(tw)].sum() - len(tw) * (len(tw) + 1) / 2) / (len(tw) * len(dw)), n))
    print('   AUC T vs D winners:', ', '.join(f'{n} {a:.2f}' for a, n in sorted(au, reverse=True)[:6]))
