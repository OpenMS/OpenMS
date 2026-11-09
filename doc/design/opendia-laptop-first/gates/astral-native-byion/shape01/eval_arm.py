"""Evaluate one native stage-2 arm (tool stage2_score): entrapment-matched N at FDP_c, plain pooled paired q and bulk-null
weighted q (HGB on the arm's own stage-2 features, bulknull_native method), FDP_c combined. Also prints the fold weights.
Usage: python eval_arm.py <arm dir>"""
import sys, re, sqlite3, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; R = 148983 / 151239; ARM = sys.argv[1]
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
p = pd.read_sql('select ID precursor_id, TRAML_ID, DECOY from PRECURSOR', c)
p['pair_key'] = np.where(p.DECOY == 1, p.TRAML_ID.str.replace('^DECOY_', '', regex=True), p.TRAML_ID)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c)
f = pd.read_csv(ARM + '/features.tsv', sep='\t').merge(p[['precursor_id', 'pair_key']], on='precursor_id').merge(ent, on='precursor_id', how='left')
FEAT = [x for x in f.columns if x not in ('sample_id', 'precursor_id', 'is_decoy', 'fold', 'stage2_score', 'pair_key', 'e')]
prec = pd.factorize(f.precursor_id)[0]; pair = pd.factorize(f.pair_key)[0]; dec = f.is_decoy.values == 1; fold = f.fold.values
entr = f.e.fillna(0).values.astype(bool); X = f[FEAT].values.astype(np.float64); s = f.stage2_score.values.astype(np.float64)
def winners(idx, s):
    o = idx[np.lexsort((idx, -s[idx]))]; _, first = np.unique(prec[o], return_index=True); return o[first]
w = winners(np.arange(len(s)), s); z = s.copy()
for k in np.unique(fold):
    dk = w[dec[w] & (fold[w] == k)]; med = np.median(s[dk]); mad = 1.4826 * np.median(np.abs(s[dk] - med))
    z[fold == k] = (s[fold == k] - med) / mad
w = winners(np.arange(len(z)), z); tg = w[~dec[w]]; o = tg[np.argsort(-z[tg], kind='stable')]
fd = np.cumsum(entr[o]) * (1 + 1 / R) / np.arange(1, len(o) + 1)
mt = [int(np.where(fd <= t)[0].max() + 1) if (fd <= t).any() else 0 for t in (0.005, 0.01, 0.015, 0.02)]
wdv = w[dec[w]]; twv = w[~dec[w]]; thr = np.median(z[wdv]); bt, bd = twv[z[twv] < thr], wdv[z[wdv] < thr]
out = [f'N at FDP_c 0.5/1/1.5/2%: {mt}']
for seed in (0, 1, 2):
    r = np.ones(len(z)); pr = np.zeros(len(z))
    for k in np.unique(fold):
        trn = np.r_[bt[fold[bt] != k], bd[fold[bd] != k]]
        m = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.05, max_leaf_nodes=31, random_state=seed).fit(X[trn], (~dec[trn]).astype(int))
        te = wdv[fold[wdv] == k]; pr[te] = m.predict_proba(X[te])[:, 1]
    r[wdv] = np.clip(pr[wdv] / np.maximum(1 - pr[wdv], 1e-9) * len(bd) / len(bt), 0.2, 5); r[wdv] /= r[bd].mean()
    for nm, wv in ((('plain', np.ones(len(z))),) if seed == 0 else ()) + ((f'bulk-null s{seed}', r),):
        d = pd.DataFrame({'i': w, 'p': pair[w], 's': z[w], 'd': dec[w]}).sort_values(['s', 'd'], ascending=[False, False]).drop_duplicates('p').sort_values('s', ascending=False)
        dd = d.d.values; T = np.cumsum(~dd); Dn = np.cumsum(np.where(dd, wv[d.i.values], 0.0))
        q = np.minimum.accumulate(np.minimum(1, (Dn + 1) / np.maximum(T, 1))[::-1])[::-1]; kept = d.i.values
        qs = []
        for qq in (0.01, 0.05):
            a = kept[(~dec[kept]) & (q <= qq)]; qs.append(f'q<={qq}: {len(a)} FDP_c {entr[a].sum() * (1 + 1 / R) / max(1, len(a)):.3f}')
        out.append(f'{nm} ' + ' '.join(qs))
print(ARM); print('\n'.join(out))
log = open(ARM + '/stdout.log').read()
for line in re.findall(r'second stage fold=.*', log): print(' ', line[:900])
