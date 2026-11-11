# Sequence "decoy tell" check: can sequence-only features separate target from decoy library entries
# (and entrapment from decoy)? Library only, no run data. 3-fold by pair, logistic regression.
import sqlite3, numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
c = sqlite3.connect('file:astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
p = pd.read_sql("""select p.ID pid, p.TRAML_ID t, p.DECOY d, p.CHARGE z, pe.UNMODIFIED_SEQUENCE s,
 min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
 join PEPTIDE pe on pe.ID=m.PEPTIDE_ID join PEPTIDE_PROTEIN_MAPPING pm on pm.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pm.PROTEIN_ID group by p.ID""", c)
AA = 'ACDEFGHIKLMNPQRSTVWY'
def f(s):
    v = [len(s)] + [s.count(a) / len(s) for a in AA]
    for pos in (0, 1, -2, -1): v += [float(s[pos] == a) for a in AA]
    v += [s.count(a + b) / len(s) for a, b in ('KP', 'RP', 'DP', 'PP', 'GG', 'EE', 'KK', 'LL')]
    return v
p['key'] = p.t.str.replace('^DECOY_', '', regex=True); p['fold'] = p.key.map(hash) % 3
X = np.array([f(s) for s in p.s]); y = (p.d == 0).astype(int).values; sc = np.zeros(len(p))
for k in range(3):
    m = LogisticRegression(max_iter=2000).fit(X[p.fold != k], y[p.fold != k]); sc[p.fold == k] = m.predict_proba(X[p.fold == k])[:, 1]
d = p.d.values == 1
print('target(all) vs decoy AUC', round(roc_auc_score(y, sc), 3))
print('human vs decoy', round(roc_auc_score(np.r_[np.ones(((~d) & (p.e == 0)).sum()), np.zeros(d.sum())], np.r_[sc[(~d) & (p.e == 0)], sc[d]]), 3),
      ' entrapment vs decoy', round(roc_auc_score(np.r_[np.ones(((~d) & (p.e == 1)).sum()), np.zeros(d.sum())], np.r_[sc[(~d) & (p.e == 1)], sc[d]]), 3))
print('C-term residue share  target:', p[~d].s.str[-1].value_counts(normalize=True).head(3).round(3).to_dict(), ' decoy:', p[d].s.str[-1].value_counts(normalize=True).head(3).round(3).to_dict())
print('N-term residue P share target', round((p[~d].s.str[0] == 'P').mean(), 4), 'decoy', round((p[d].s.str[0] == 'P').mean(), 4))
