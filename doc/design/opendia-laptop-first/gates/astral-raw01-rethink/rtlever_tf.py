"""Validation of the RT lever (rtlever_full.py) with a TEACHER-FREE base score and robustness variants.

Base score: the tool's own stage-2 score (search+boundary LDA, cross-fitted, decoy-trained; no DIA-NN information).
RT term: log N(rt_apex - pred; sigma) added with weight 1 in units of the per-fold decoy-winner MAD of the base score
is NOT appropriate (different scales), so the base score is first converted to a log-odds proxy by standardising
with the decoy-winner median/MAD (as the tool's pooled TDC does); the RT term is in nats. Candidates limited to +-3 sigma.
Predictor noise: Gaussian or Student-t (3 dof, same scale), 3 seeds. Predicted RT = DIA-NN RT + noise for confirmed
targets, calibrated library RT for all others (same optimistic asymmetry as before; reported as such).
Also reports the 'sharpen everyone with the library RT' control: sigma shrunk but centres = library RT for all,
which shows the effect of a narrower window without better centres (should not help).
"""
import sqlite3, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; D = B + '/astral-design-raw01-20261107'
f = pd.read_csv(D + '/stage2-boundary-design/boundary/features.tsv', sep='\t', usecols=['sample_id', 'precursor_id', 'is_decoy', 'fold', 'stage2_score'])
cr = pd.read_csv(D + '/stage2-boundary-design/boundary/crossfit.tsv', sep='\t', usecols=['sample_id', 'in_rt_window'])
s = pd.read_csv(D + '/pipeline/search.features.tsv', sep='\t', usecols=['candidate_row', 'rt_apex', 'apex_rt_error_seconds']).rename(columns={'candidate_row': 'sample_id'})
df = f.merge(cr, on='sample_id').merge(s, on='sample_id'); df['cal_rt'] = df.rt_apex - df.apex_rt_error_seconds
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
p = pd.read_sql('select ID precursor_id, TRAML_ID, DECOY from PRECURSOR', c)
p['pair_key'] = np.where(p.DECOY == 1, p.TRAML_ID.str.replace('^DECOY_', '', regex=True), p.TRAML_ID)
df = df.merge(p[['precursor_id', 'pair_key']], on='precursor_id')
ref = pd.read_parquet(B + '/astral-diann-gap-20261108/ref_exact.parquet'); ref = ref[ref.Run.str.endswith('_01')]
df = df.merge(ref[['precursor_id']].assign(diann_rt=ref.RT.values * 60), on='precursor_id', how='left')
# standardise base score per fold by decoy-winner median/MAD
w0 = df.sort_values('stage2_score', ascending=False).drop_duplicates('precursor_id')
df['z'] = 0.0
for k, g in w0[w0.is_decoy == 1].groupby('fold'):
    med = g.stage2_score.median(); mad = 1.4826 * (g.stage2_score - med).abs().median()
    df.loc[df.fold == k, 'z'] = (df.loc[df.fold == k, 'stage2_score'] - med) / mad

def tdc(t, col):
    w = t.sort_values(col, ascending=False).drop_duplicates('precursor_id')
    w = w.sort_values([col, 'is_decoy'], ascending=[False, False]).drop_duplicates('pair_key').sort_values(col, ascending=False)
    d = w.is_decoy.values; T = np.cumsum(1 - d); Dc = np.cumsum(d)
    q = np.minimum.accumulate((np.minimum(1, (Dc + 1) / np.maximum(T, 1)))[::-1])[::-1]
    a = w[(d == 0) & (q <= 0.01)]
    return len(a), int(((a.rt_apex - a.diann_rt).abs() <= 10).sum())

print('base teacher-free (in_rt_window only, no RT term):', tdc(df[df.in_rt_window == 1], 'z'))
pr0 = df[['precursor_id', 'cal_rt', 'diann_rt']].drop_duplicates('precursor_id')
for noise in ('gauss', 't3'):
    for sigma in (24.0, 12.0, 6.0, 3.0):
        res = []
        for seed in range(3):
            rng = np.random.default_rng(seed)
            e = rng.normal(0, sigma, len(pr0)) if noise == 'gauss' else sigma * rng.standard_t(3, len(pr0))
            pred = np.where(pr0.diann_rt.notna(), pr0.diann_rt + e, pr0.cal_rt) if sigma < 24 else pr0.cal_rt.values
            t = df.merge(pr0[['precursor_id']].assign(pred=pred), on='precursor_id')
            zz = (t.rt_apex - t.pred) / sigma
            t = t[zz.abs() <= 3].assign(fs=lambda d_: d_.z - 0.5 * ((d_.rt_apex - d_.pred) / sigma) ** 2)
            res.append(tdc(t, 'fs'))
            if sigma == 24: break
        print(f'{noise:5s} sigma {sigma:4.1f}: (accepted, confirmed-right) per seed {res}', flush=True)
# control: narrower window, library centres for everyone
for sigma in (12.0, 6.0):
    t = df.assign(pred=df.cal_rt); zz = (t.rt_apex - t.pred) / sigma
    t = t[zz.abs() <= 3].assign(fs=lambda d_: d_.z - 0.5 * ((d_.rt_apex - d_.pred) / sigma) ** 2)
    print(f'control library centres, sigma {sigma}: {tdc(t, "fs")}')
