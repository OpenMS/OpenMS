"""V18 charge-sibling features, pure helpers (label-free). See PLAN_V18.md."""
import numpy as np

TOL, CLIP = 3.0, (-5.0, 20.0)

def fold_z(full, w_scores, precursor_id, sample_id, is_decoy, fold):
    """Per-fold decoy-winner median/MAD normalisation of the full score of every window.
    Winners: each precursor's best window under w_scores (-1e30 rows never win), ties by sample_id."""
    full = np.asarray(full, np.float64); z = np.full(len(full), np.nan)
    o = np.lexsort((sample_id, -np.asarray(w_scores, np.float64)))
    _, first = np.unique(np.asarray(precursor_id)[o], return_index=True); win = o[first]
    for k in np.unique(fold):
        d = win[(fold[win] == k) & (is_decoy[win] == 1) & (w_scores[win] > -1e29)]
        med = np.median(full[d]); mad = 1.4826 * np.median(np.abs(full[d] - med))
        if not np.isfinite(mad) or mad <= 0: raise ValueError(f'degenerate MAD in fold {k}')
        m = fold == k; z[m] = (full[m] - med) / mad
    assert np.isfinite(z).all()
    return z

def sib_features(precursor_id, peptide_id, rt, z, tol=TOL, clip=CLIP):
    """For each window: has_sib (another precursor of the same peptide is in the population), has_near (a sibling
    window has apex within tol s), sib_z_near (max sibling z within tol, clipped; 0 when has_near is 0)."""
    n = len(z); out = np.zeros((n, 3), np.float64)
    pep = np.asarray(peptide_id); o = np.argsort(pep, kind='stable'); ps = pep[o]
    starts = np.r_[0, np.flatnonzero(np.diff(ps)) + 1]; ends = np.r_[starts[1:], n]
    for a, b in zip(starts, ends):
        g = o[a:b]; pr = np.asarray(precursor_id)[g]
        if len(np.unique(pr)) < 2: continue
        other = pr[:, None] != pr[None, :]
        near = other & (np.abs(rt[g][:, None] - rt[g][None, :]) <= tol)
        out[g, 0] = 1.0
        hn = near.any(1); out[g, 1] = hn
        zz = np.where(near, z[g][None, :], -np.inf).max(1)
        out[g, 2] = np.where(hn, np.clip(zz, *clip), 0.0)
    return out
