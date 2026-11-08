"""Predeclared PLAN.md metrics for one crossfit.tsv (target/decoy labels only)."""
import sys
import numpy as np
import pandas as pd


def paired_pooled(df, score):
    """Pooled paired TDC exactly as the tool: per-fold decoy-winner median/MAD standardisation, then (D+1)/T."""
    w = df[df.in_rt_window == 1]
    w = w.sort_values([score, 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id')
    z = w[score].copy()
    for f, g in w[w.is_decoy == 1].groupby('fold'):
        med = g[score].median()
        mad = 1.4826 * (g[score] - med).abs().median()
        z[w.fold == f] = (w.loc[w.fold == f, score] - med) / mad
    w = w.assign(z=z.values)
    w['key'] = w.pair_key
    best = w.sort_values(['z', 'is_decoy'], ascending=[False, False]).drop_duplicates('key')
    best = best.sort_values('z', ascending=False)
    # tie groups
    d = best.is_decoy.values
    zz = best.z.values
    T = np.cumsum(1 - d)
    D = np.cumsum(d)
    # last index per tie value
    last = np.r_[zz[1:] != zz[:-1], True]
    q = np.where(T > 0, np.minimum(1, (D + 1) / np.maximum(T, 1)), 1.0)
    qt = pd.Series(q).where(last).bfill().values
    qt = np.minimum.accumulate(qt[::-1])[::-1]
    best = best.assign(q=qt)
    return best, w


def main(path, lib_pairs):
    df = pd.read_csv(path, sep='\t')
    df = df.merge(lib_pairs, on='precursor_id', how='left')
    out = {}
    for score in ['logit', 'stage2_score']:
        best, w = paired_pooled(df, score)
        t = best[best.is_decoy == 0]
        for cut in (0.005, 0.01, 0.02, 0.05):
            out[f'{score}_q{cut}'] = int((t.q <= cut).sum())
        out[f'{score}_q0.01_per_fold'] = t[t.q <= 0.01].groupby('fold').size().to_dict()
    # Null-symmetry R on CNN precursor winners (in window), per fold below decoy-winner median.
    w = df[df.in_rt_window == 1].sort_values(['logit', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id')
    tl = dl = 0
    for f, g in w.groupby('fold'):
        med = g[g.is_decoy == 1].logit.median()
        low = g[g.logit < med]
        tl += int((low.is_decoy == 0).sum())
        dl += int((low.is_decoy == 1).sum())
        out[f'decoy_z999_fold{f}'] = None
    nt, nd = int((w.is_decoy == 0).sum()), int((w.is_decoy == 1).sum())
    out['T_low'], out['D_low'], out['N_T'], out['N_D'] = tl, dl, nt, nd
    out['R'] = round((tl / dl) / (nt / nd), 4)
    tl5 = dl5 = 0
    for f, g in w.groupby('fold'):
        c = g[g.is_decoy == 1].logit.quantile(0.05)
        tl5 += int(((g.logit < c) & (g.is_decoy == 0)).sum())
        dl5 += int(((g.logit < c) & (g.is_decoy == 1)).sum())
    out['R_q05'] = round((tl5 / dl5) / (nt / nd), 4)
    out['guard_R_median_ge_0.83'] = out['R'] >= 0.83
    for f, g in w[w.is_decoy == 1].groupby('fold'):
        med = g.logit.median()
        mad = 1.4826 * (g.logit - med).abs().median()
        out[f'decoy_z999_fold{f}'] = round(float(((g.logit - med) / mad).quantile(0.999)), 2)
    for k, v in out.items():
        print(f'{k}\t{v}')


if __name__ == '__main__':
    import sqlite3
    c = sqlite3.connect('file:' + sys.argv[2] + '?mode=ro', uri=True)
    p = pd.read_sql('select ID precursor_id, TRAML_ID, DECOY from PRECURSOR', c)
    p['pair_key'] = np.where(p.DECOY == 1, p.TRAML_ID.str.replace('^DECOY_', '', regex=True), p.TRAML_ID)
    main(sys.argv[1], p[['precursor_id', 'pair_key']])
