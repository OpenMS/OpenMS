"""Single-subscore rankings (no learning) on the XIC-maxima candidates, to sanity-check the 'minus xcorr' collapse and
see which subscores alone separate. Sign oriented a priori (higher = better unless listed in NEG)."""
import sys, numpy as np
sys.path.insert(0, '.')
import oswlike_score as S
NEG = {'xcorr_coel', 'xcorr_coel_w', 'lib_manhattan', 'lib_rmsd', 'log_width', 'norm_rt', 'ms1_xcorr_coel'}
for j, n in enumerate(S.NAMES):
    sc = np.nan_to_num(S.X0[:, j], nan=-1e9 if n not in NEG else 1e9)
    S.report(f'single {n}', -sc if n in NEG else sc)
# also check the minus-xcorr LDA: how many positives at the first selection step
cols = [c for c in S.ALL if c not in S.GROUPS['xcorr']]
X = S.X0[:, [S.NAMES.index(c) for c in cols]]
Z = np.nan_to_num((X - np.nanmean(X, 0)) / (np.nanstd(X, 0) + 1e-9))
tr = np.where(S.fold != 0)[0]
sc = Z[:, cols.index('lib_dotprod')]
for fdr in [0.15, 0.02, 0.02]:
    pos, neg = S.select(tr, sc, fdr); print('minus-xcorr SS step: fdr', fdr, 'positives', len(pos), 'negatives', len(neg))
    if len(pos) < 3: break
    sc = Z @ S.lda_fit(Z, pos, neg)
