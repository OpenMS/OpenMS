"""Nested label selection (predeclared single variant): in the semi-supervised loop, positives/negatives for the
training folds are chosen from INNER out-of-fold scores (each training fold scored by a model fitted on the other
training fold), never from the model's own in-sample scores. Same MLP as lag.py, 5-seed ensemble, entrapment FDP."""
import numpy as np, torch
src = open('overfit.py').read().split("SB = g['SEARCH']")[0]
exec(src)
fold_of, folds, select, lda, mlp_fit = g['fold_of'], g['folds'], g['select'], g['lda'], g['mlp_fit']

def fit(kind, X, pos, neg):
    return lda(X, pos, neg) if kind == 'lda' else mlp_fit(X, pos, neg)

def crossfit_nested(cols):
    X = df[cols].to_numpy(np.float64)
    out = np.full(len(df), np.nan)
    for k in folds:
        inner = [f for f in folds if f != k]
        tr = np.where(fold_of != k)[0]
        score = X[:, 0].copy()
        for kind, fdr in [('lda', 0.15), ('lda', 0.05), ('lda', 0.05), ('mlp', 0.05), ('mlp', 0.05)]:
            oof = np.full(len(df), np.nan)
            for j in inner:            # score inner fold j with a model trained on the other inner fold(s)
                trj = np.where((fold_of != k) & (fold_of != j))[0]
                pos, neg = select(trj, score, fdr)
                oof[fold_of == j] = fit(kind, X, pos, neg)[fold_of == j]
                for jj in inner:       # put inner folds on one scale (decoy-winner median/MAD), as pooled TDC does
                    pass
            for j in inner:
                m = fold_of == j
                d = m & g['is_dec']
                med = np.median(oof[d]); mad = 1.4826 * np.median(np.abs(oof[d] - med)) + 1e-9
                oof[m] = (oof[m] - med) / mad
            score = np.where(fold_of != k, oof, score)
            last = (kind, fdr)
        pos, neg = select(tr, score, 0.05)   # final model: all training folds, labels from out-of-fold scores
        out[fold_of == k] = mlp_fit(X, pos, neg)[fold_of == k]
    return out

SB = g['SEARCH'] + g['BOUND']
acc = 0
for s in range(10, 15):
    torch.manual_seed(s); np.random.seed(s)
    sc = crossfit_nested(SB)
    for k in folds:
        m = fold_of == k
        sc[m] = (sc[m] - sc[m].mean()) / sc[m].std()
    acc = acc + sc
tmp = df[['sample_id', 'precursor_id', 'is_decoy', 'fold', 'in_rt_window', 'pair_key', 'all_ent']].assign(score=acc)
best, w = g['paired_pooled'](tmp, 'score')
nt, nd = (w.is_decoy == 0).sum(), (w.is_decoy == 1).sum()
low = w[w.z < w[w.is_decoy == 1].z.median()]
R = ((low.is_decoy == 0).sum() / (low.is_decoy == 1).sum()) / (nt / nd)
print(f"{'MLP nested-select ens5 s+b':34s} {fdp(best)} | R {R:.3f}")
