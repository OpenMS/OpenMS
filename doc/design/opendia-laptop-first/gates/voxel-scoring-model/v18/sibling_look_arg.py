"""V18 pre-plan look, raw01 only, descriptive. Charge-sibling co-elution.
For each precursor winner (D11FWT window), find the other charge states of the same peptide in the
population and ask whether they have a candidate window near the same apex, and how well it scores.
Label-free features; group labels joined only for summary. raw02 is not read."""
import sqlite3, sys, numpy as np, pandas as pd
V = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2'
L = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-entrapment-20261107/lib/library.pqp'
raw = sys.argv[1]
c = sqlite3.connect(f'file:{L}?mode=ro', uri=True)
m = pd.read_sql("select PRECURSOR_ID precursor_id, PEPTIDE_ID peptide_id from PRECURSOR_PEPTIDE_MAPPING", c)
P = pd.read_csv(f'{V}/{raw}/population.tsv', sep='\t')
full = np.load(f'{V}/{raw}/scores_D11F.npy'); lin = np.load(f'{V}/{raw}/lin_D11F.npy')
P['full'] = full; P['lin'] = lin
P = P.merge(m, on='precursor_id', how='left', validate='many_to_one')
assert P.peptide_id.notna().all()
# per-fold decoy-winner normalisation of full score, label-free (decoy flag is a search label, not an evaluation label)
W = pd.read_csv(f'{V}/v4/stageD_v16/score/reference_winners_{raw}.tsv.gz', sep='\t')
W = W.merge(m, on='precursor_id', how='left')
grp = np.where(W.stage.eq('D right window, below cut'), 'D', np.where(W.stage.eq('E accepted'), 'E', W.group))
W['grp'] = grp
# Normalise every candidate window's full score to the decoy-winner distribution of its fold.
for f, g in W[W.is_decoy == 1].groupby('fold'):
    med = g.full.median(); mad = 1.4826 * (g.full - med).abs().median()
    P.loc[P.fold == f, 'z'] = (P.loc[P.fold == f, 'full'] - med) / mad
by_pep = {k: v for k, v in P.groupby('peptide_id')}
rows = []
for r in W.itertuples():
    g = by_pep[r.peptide_id]; s = g[g.precursor_id != r.precursor_id]
    if len(s) == 0:
        rows.append((np.nan, np.nan, np.nan, 0)); continue
    d = np.abs(s.rt_apex.values - r.rt_apex)
    near = d <= 3.0
    best_near = s.z.values[near].max() if near.any() else -np.inf
    # where does the sibling's own best window sit relative to this apex
    sb = s.sort_values(['z', 'sample_id'], ascending=[False, True]).iloc[0]
    rows.append((best_near, abs(sb.rt_apex - r.rt_apex), near.sum(), len(s)))
S = pd.DataFrame(rows, columns=['sib_z_near', 'sib_best_dist', 'sib_n_near', 'sib_n'])
W = pd.concat([W.reset_index(drop=True), S], axis=1)
W.to_csv(f'sibling_{raw}.tsv.gz', sep='\t', index=False)
has = W.sib_n > 0
def auc(a, b):
    a = a[np.isfinite(a)] if False else a; x = np.r_[a, b]; rk = pd.Series(x).rank().values
    return (rk[:len(a)].sum() - len(a) * (len(a) + 1) / 2) / (len(a) * len(b))
out = []
for gname in ['D', 'E', 'H', 'ENT', 'DEC']:
    k = has & (W.grp == gname)
    x = W.loc[k]
    out.append(dict(group=gname, n=int((W.grp == gname).sum()), with_sibling=int(k.sum()),
                    sib_best_within3s=float((x.sib_best_dist <= 3).mean()),
                    sib_any_window_within3s=float((x.sib_n_near > 0).mean()),
                    sib_z_near_median=float(np.median(np.where(np.isfinite(x.sib_z_near), x.sib_z_near, -50)))))
print(pd.DataFrame(out).to_string(index=False))
# AUCs vs decoy winners with siblings, missing near window -> lowest value (tie handled by rank)
dec = W.loc[has & (W.grp == 'DEC'), 'sib_z_near'].replace(-np.inf, -50).values
for gname in ['D', 'E', 'H', 'ENT']:
    a = W.loc[has & (W.grp == gname), 'sib_z_near'].replace(-np.inf, -50).values
    print(gname, 'vs DEC AUC sib_z_near', round(auc(a, dec), 4))
# score-stratified: within deciles of own full_pct among decoys
W['stratum'] = pd.cut(W.full_pct, [0, .5, .9, .99, .999, 1.0001], right=False)
for st, g in W[has].groupby('stratum', observed=True):
    dd = g.loc[g.grp == 'DEC', 'sib_z_near'].replace(-np.inf, -50).values
    line = [str(st), f'DEC n={len(dd)}']
    for gname in ['D', 'H', 'ENT']:
        a = g.loc[g.grp == gname, 'sib_z_near'].replace(-np.inf, -50).values
        if len(a) >= 20 and len(dd) >= 20: line.append(f'{gname} n={len(a)} AUC={auc(a, dd):.3f}')
    print(' | '.join(line))
