"""Bulk-null density-ratio TDC on the NATIVE raw 01 stage-2 output (declared before running).
Input: byion/crossfit.tsv + byion/features.tsv (search+boundary+byion, commit 8c49a5bbb2). Score = tool stage2_score,
standardised per fold by decoy precursor winners (median, 1.4826 MAD) as in the tool's pooled TDC.
Method (as CALIB section 9): winners = per-precursor best in RT window; bulk = winners below the median decoy-winner z.
HistGradientBoosting (200 iter, lr 0.05, 31 leaves) on the 16 stage-2 features, bulk target winners (1) vs bulk decoy
winners (0), cross-fitted by the tool's fold. Decoy weight p/(1-p)*n_d/n_t clipped to [0.2, 5], renormalised to mean 1 over
bulk decoys. Weighted paired TDC: q = (sum kept decoy weights + 1)/kept targets, monotone.
Rule: PASS if weighted q<=0.05 FDP_c <= 0.06 AND weighted q<=0.01 FDP_c <= 0.015 for >= 2 of 3 classifier seeds.
Entrapment labels are used only to evaluate."""
import sqlite3, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; D = B + '/astral-design-raw01-20261107/stage2-byion-design/byion'
R = 148983 / 151239
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
p = pd.read_sql('select ID precursor_id, TRAML_ID, DECOY from PRECURSOR', c)
p['pair_key'] = np.where(p.DECOY == 1, p.TRAML_ID.str.replace('^DECOY_', '', regex=True), p.TRAML_ID)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c)
f = pd.read_csv(D + '/features.tsv', sep='\t')
FEAT = [x for x in f.columns if x not in ('sample_id', 'precursor_id', 'is_decoy', 'fold', 'stage2_score')]
w = f.sort_values(['stage2_score', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').copy()
w['z'] = 0.0
for k, g in w[w.is_decoy == 1].groupby('fold'):
    med = g.stage2_score.median(); mad = 1.4826 * (g.stage2_score - med).abs().median()
    w.loc[w.fold == k, 'z'] = (w.loc[w.fold == k, 'stage2_score'] - med) / mad
w = w.merge(p[['precursor_id', 'pair_key']], on='precursor_id').merge(ent, on='precursor_id', how='left')
w['e'] = w.e.fillna(0).astype(int)
dec = w.is_decoy.values == 1; X = w[FEAT].values; z = w.z.values; fold = w.fold.values
thr = np.median(z[dec]); bulk = z < thr; nt, nd = int((bulk & ~dec).sum()), int((bulk & dec).sum())
def tdc(wdec):
    d = w.assign(wd=wdec).sort_values(['z', 'is_decoy'], ascending=[False, False]).drop_duplicates('pair_key')
    dd = d.is_decoy.values == 1; T = np.cumsum(~dd); Dn = np.cumsum(np.where(dd, d.wd.values, 0.0))
    q = np.minimum.accumulate(np.minimum(1, (Dn + 1) / np.maximum(T, 1))[::-1])[::-1]
    out = []
    for qq in (0.01, 0.05):
        sel = (~dd) & (q <= qq); n = int(sel.sum()); ne = int(d.e.values[sel].sum())
        out.append(f'q<={qq}: {n} FDP_c {ne * (1 + 1 / R) / max(n, 1):.3f}')
    return ' | '.join(out)
print(f'winners {len(w)}, bulk targets {nt}, bulk decoys {nd}, features {FEAT}', flush=True)
print('plain    ', tdc(np.ones(len(w))), flush=True)
for seed in (0, 1, 2):
    pr = np.zeros(len(w)); idx = np.where(bulk)[0]; y = (~dec[idx]).astype(int)
    for k in np.unique(fold):
        tr = idx[fold[idx] != k]
        m = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.05, max_leaf_nodes=31, random_state=seed).fit(X[tr], (~dec[tr]).astype(int))
        te = np.where(dec & (fold == k))[0]; pr[te] = m.predict_proba(X[te])[:, 1]
    r = np.ones(len(w)); r[dec] = np.clip(pr[dec] / np.maximum(1 - pr[dec], 1e-9) * nd / nt, 0.2, 5)
    r[dec] = r[dec] / r[dec & bulk].mean()
    print(f'weighted seed {seed}', tdc(r), flush=True)
