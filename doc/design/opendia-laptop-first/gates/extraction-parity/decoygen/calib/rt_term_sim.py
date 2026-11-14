"""Calibration probe: stage-3 anchors at q <= 0.01 when the paired competition ranks by the co-elution score plus an RT
term, from a diag.tsv (label-free; only target/decoy). RT term: log-likelihood ratio of the apex residual to the window
centre under a Gaussian of width sigma = half_window / k against a uniform over the window, weighted by w.
Usage: rt_term_sim.py <library.pqp> <diag.tsv> <half_window_s>"""
import sys, sqlite3, numpy as np, pandas as pd
lib, diag, hw = sys.argv[1], sys.argv[2], float(sys.argv[3])
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
pre = pd.read_sql('select ID precursor_id, TRAML_ID traml, DECOY decoy from PRECURSOR', c)
tid = pre[pre.decoy == 0].set_index('traml').precursor_id
pre['pair'] = np.where(pre.decoy == 1, pre.traml.str.replace('^DECOY_', '', regex=True).map(tid), pre.precursor_id)
d = pd.read_csv(diag, sep='\t', comment='#').merge(pre[['precursor_id', 'pair']], on='precursor_id')
d = d[d.stage3_score.notna()].copy()
d['res'] = d.stage3_rt - d.stage3_centre_rt


def anchors(s):
    d['s'] = s
    w = d.sort_values('s', ascending=False).drop_duplicates('pair')  # pair winner (ties: first, as sort is stable-ish)
    w = w.sort_values('s', ascending=False)
    D = np.cumsum(w.decoy.values); T = np.cumsum(1 - w.decoy.values)
    q = np.minimum.accumulate(np.minimum(1, (D + 1) / np.maximum(T, 1))[::-1])[::-1]
    return int(((q <= 0.01) & (w.decoy.values == 0)).sum())


print('score only:', anchors(d.stage3_score.values))
sigma1 = float(sys.argv[4]) if len(sys.argv) > 4 else None
if sigma1:
    # Parameter-free: score as evidence against the decoy score distribution (-log empirical p-value among decoys),
    # plus the RT log-likelihood ratio of a Gaussian with the stage-1 RT sigma against a uniform over the window.
    dec = np.sort(d[d.decoy == 1].stage3_score.values)
    p = (len(dec) - np.searchsorted(dec, d.stage3_score.values, side='left') + 1) / (len(dec) + 1)
    llr = -0.5 * (d.res / sigma1) ** 2 - np.log(sigma1 * np.sqrt(2 * np.pi)) + np.log(2 * hw)
    print(f'-log p(score | decoys) only: {anchors(-np.log(p))}')
    print(f'-log p + RT LLR (stage-1 sigma {sigma1:.1f} s): {anchors(-np.log(p) + llr.values)}')
    sys.exit()
for k in (2, 3, 4):
    sig = hw / k
    llr = -0.5 * (d.res / sig) ** 2 - np.log(sig * np.sqrt(2 * np.pi)) + np.log(2 * hw)
    for wt in (0.05, 0.1, 0.2):
        print(f'sigma = hw/{k}, weight {wt}: anchors {anchors(d.stage3_score.values + wt * llr.values)}')
