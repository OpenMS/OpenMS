"""h15a: which density-ratio classifier should native bulk-null use? (declared before running)
Same bulk-null method as bulknull_native.py, on the native search+boundary+byion features of raws 01/02, tool stage2_score
standardised per fold by decoy winners. Only the classifier for p(target | bulk winner, x) changes:
  hgb     : sklearn HistGradientBoosting 200 iter, 31 leaves (reference; no native equivalent)
  stumps  : gradient-boosted depth-1 trees, 100 rounds, lr 0.1, 32 quantile bins (native port is about 100 lines)
  depth2  : the same with depth-2 trees (HGB max_depth=2, 100 iter)
  logit   : logistic regression on standardised features (failed in CALIB section 9; control)
  none    : plain TDC
Rule: choose the SIMPLEST classifier (order logit < stumps < depth2 < hgb) whose bulk-null q<=0.05 FDP_c is within 0.005
of hgb on BOTH raws AND whose q<=0.01 FDP_c is <= 0.015 on both. Usage: python ratio_classifier.py"""
import sqlite3, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; R = 148983 / 151239
ARMS = {'raw01': B + '/astral-design-raw01-20261107/stage2-byion-design/byion', 'raw02': B + '/astral-raw02-byion-design-20261109/byion'}
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
p = pd.read_sql('select ID precursor_id, TRAML_ID, DECOY from PRECURSOR', c)
p['pair_key'] = np.where(p.DECOY == 1, p.TRAML_ID.str.replace('^DECOY_', '', regex=True), p.TRAML_ID)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c)
def make(kind):
    if kind == 'hgb': return HistGradientBoostingClassifier(max_iter=200, learning_rate=0.05, max_leaf_nodes=31, random_state=0)
    if kind == 'stumps': return HistGradientBoostingClassifier(max_iter=100, learning_rate=0.1, max_depth=1, max_bins=32, random_state=0)
    if kind == 'depth2': return HistGradientBoostingClassifier(max_iter=100, learning_rate=0.1, max_depth=2, max_bins=32, random_state=0)
    return LogisticRegression(max_iter=1000)
for run, arm in ARMS.items():
    f = pd.read_csv(arm + '/features.tsv', sep='\t').merge(p[['precursor_id', 'pair_key']], on='precursor_id').merge(ent, on='precursor_id', how='left')
    FEAT = [x for x in f.columns if x not in ('sample_id', 'precursor_id', 'is_decoy', 'fold', 'stage2_score', 'pair_key', 'e')]
    prec = pd.factorize(f.precursor_id)[0]; pair = pd.factorize(f.pair_key)[0]; dec = f.is_decoy.values == 1; fold = f.fold.values
    entr = f.e.fillna(0).values.astype(bool); X = f[FEAT].values.astype(np.float64); s = f.stage2_score.values.astype(np.float64)
    def winners(idx, s):
        o = idx[np.lexsort((idx, -s[idx]))]; _, first = np.unique(prec[o], return_index=True); return o[first]
    w = winners(np.arange(len(s)), s); z = s.copy()
    for k in np.unique(fold):
        dk = w[dec[w] & (fold[w] == k)]; med = np.median(s[dk]); z[fold == k] = (s[fold == k] - med) / (1.4826 * np.median(np.abs(s[dk] - med)))
    w = winners(np.arange(len(z)), z); wdv = w[dec[w]]; twv = w[~dec[w]]; thr = np.median(z[wdv]); bt, bd = twv[z[twv] < thr], wdv[z[wdv] < thr]
    Xs = (X - X[np.r_[bt, bd]].mean(0)) / (X[np.r_[bt, bd]].std(0) + 1e-9)
    d0 = pd.DataFrame({'i': w, 'p': pair[w], 's': z[w], 'd': dec[w]}).sort_values(['s', 'd'], ascending=[False, False]).drop_duplicates('p').sort_values('s', ascending=False)
    for kind in ('none', 'logit', 'stumps', 'depth2', 'hgb'):
        r = np.ones(len(z))
        if kind != 'none':
            pr = np.zeros(len(z))
            for k in np.unique(fold):
                trn = np.r_[bt[fold[bt] != k], bd[fold[bd] != k]]
                m = make(kind).fit(Xs[trn], (~dec[trn]).astype(int)); te = wdv[fold[wdv] == k]; pr[te] = m.predict_proba(Xs[te])[:, 1]
            r[wdv] = np.clip(pr[wdv] / np.maximum(1 - pr[wdv], 1e-9) * len(bd) / len(bt), 0.2, 5); r[wdv] /= r[bd].mean()
        dd = d0.d.values; Dn = np.cumsum(np.where(dd, r[d0.i.values], 0.0))
        q = np.minimum.accumulate(np.minimum(1, (Dn + 1) / np.maximum(np.cumsum(~dd), 1))[::-1])[::-1]; kept = d0.i.values
        qs = []
        for qq in (0.01, 0.05):
            a = kept[(~dec[kept]) & (q <= qq)]; qs.append(f'q<={qq}: {len(a)} FDP_c {entr[a].sum() * (1 + 1 / R) / max(1, len(a)):.3f}')
        print(f'{run} {kind:7s} ' + ' | '.join(qs), flush=True)
