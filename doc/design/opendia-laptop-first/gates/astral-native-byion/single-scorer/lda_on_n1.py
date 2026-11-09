"""Diagnostic (reported only, no rule): does a fold-crossed linear refit (LDA, PyProphet-style 3 rounds) on top of the
N1 TOPAZ logit plus the 11 native stage-2 heuristics restore TDC validity (entrapment-vs-decoy winner AUC ~0.5) and
keep N1's ranking? Uses the native byion features.tsv (same windows/sample_id order) with cnn_z replaced by the N1 logit."""
import sys, sqlite3, numpy as np, pandas as pd
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from scipy.stats import mannwhitneyu
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; R = 148983 / 151239
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
p = pd.read_sql('select ID precursor_id, TRAML_ID, DECOY from PRECURSOR', c)
p['pair_key'] = np.where(p.DECOY == 1, p.TRAML_ID.str.replace('^DECOY_', '', regex=True), p.TRAML_ID)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c)
FEATS = {'raw01': B + '/astral-design-raw01-20261107/stage2-byion-design/byion', 'raw02': B + '/astral-raw02-byion-design-20261109/byion'}
for raw, base in FEATS.items():
    f = pd.read_csv(base + '/features.tsv', sep='\t'); cf = pd.read_csv(base + '/crossfit.tsv', sep='\t', usecols=['sample_id', 'in_rt_window'])
    n1 = pd.read_csv(f'{raw}/N1/crossfit.tsv', sep='\t', usecols=['sample_id', 'logit', 'fold'])
    f = f.merge(cf, on='sample_id').merge(n1, on='sample_id', suffixes=('', '_n1'))
    assert (f.fold == f.fold_n1).all(), 'fold mismatch'
    f = f[f.in_rt_window == 1].merge(p[['precursor_id', 'pair_key']], on='precursor_id').merge(ent, on='precursor_id', how='left').reset_index(drop=True)
    cols = ['logit'] + [x for x in f.columns if x not in ('sample_id', 'precursor_id', 'is_decoy', 'fold', 'fold_n1', 'cnn_z', 'stage2_score', 'in_rt_window', 'pair_key', 'e', 'logit')]
    X = f[cols].values.astype(float); X = (X - X.mean(0)) / (X.std(0) + 1e-9)
    prec = pd.factorize(f.precursor_id)[0]; dec = f.is_decoy.values == 1; fold = f.fold.values; sid = f.sample_id.values
    def winners(idx, s):
        o = idx[np.lexsort((sid[idx], -s[idx]))]; _, first = np.unique(prec[o], return_index=True); return o[first]
    def tdc(w, s):
        d = pd.DataFrame({'i': w, 'p': f.pair_key.values[w], 's': s[w], 'd': dec[w]}).sort_values(['s', 'd'], ascending=[False, False]).drop_duplicates('p').sort_values('s', ascending=False)
        dd = d.d.values; q = np.minimum(1, (np.cumsum(dd) + 1) / np.maximum(np.cumsum(~dd), 1)); return d.i.values, np.minimum.accumulate(q[::-1])[::-1]
    out = np.zeros(len(f))
    for k in np.unique(fold):
        tr = np.where(fold != k)[0]; s = X[:, 0].copy()
        for it in range(3):
            w = winners(tr, s); kept, q = tdc(w, s); pos = kept[(~dec[kept]) & (q <= 0.01)]
            if len(pos) < 50: tw = w[~dec[w]]; pos = tw[np.argsort(-s[tw])[:200]]
            neg = w[dec[w]]; m = LinearDiscriminantAnalysis().fit(X[np.r_[pos, neg]], np.r_[np.ones(len(pos)), np.zeros(len(neg))])
            s = m.decision_function(X)
        te = np.where(fold == k)[0]; wd = winners(te, s); wd = wd[dec[wd]]; med = np.median(s[wd])
        out[te] = (s[te] - med) / (1.4826 * np.median(np.abs(s[wd] - med)) + 1e-12)
    w = winners(np.arange(len(f)), out); entr = f.e.fillna(0).values.astype(bool)
    tg = w[~dec[w]]; o = tg[np.argsort(-out[tg], kind='stable')]; fd = np.cumsum(entr[o]) * (1 + 1 / R) / np.arange(1, len(o) + 1)
    n1p = int(np.where(fd <= 0.01)[0].max() + 1)
    kept, q = tdc(w, out); res = []
    for qq in (0.01, 0.05):
        a = kept[(~dec[kept]) & (q <= qq)]; res.append(f'q<={qq} {len(a)} FDP_c {entr[a].sum() * (1 + 1 / R) / max(1, len(a)):.3f}')
    aucs = []
    for k in np.unique(fold):
        ww = w[fold[w] == k]; d = out[ww[dec[ww]]]; e = out[ww[(~dec[ww]) & entr[ww]]]; aucs.append(mannwhitneyu(e, d).statistic / (len(e) * len(d)))
    print(f'{raw} LDA(N1 logit + 11 heuristics): N@1% {n1p} | ' + ' | '.join(res) + ' | entrapment-vs-decoy AUC ' + ' '.join(f'{a:.3f}' for a in aucs), flush=True)
