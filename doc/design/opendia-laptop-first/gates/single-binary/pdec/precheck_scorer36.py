"""PRECHECK_SCORER36 (diagnostic): is the scorer the limit now that the default has 36 features?
Same method as PRECHECK_SCORER (24 features), on r~7 D seed 0 with the current default (evidence block on):
  A  tool score as written
  B  label-free cross-fitted GBM on the same 36 features, 3 Percolator rounds from the tool score
  B2 label-free GBM with 2x capacity (600 rounds, 63 leaves), to check learner capacity
  C  GBM trained on DIA-NN labels (positives = windows within 6 s of DIA-NN RT of DIA-NN non-entrapment q<=0.01
     IDs, negatives = decoy windows); a ceiling, not a method
Plus per-quartile AUCs of the tool score and of B. DIA-NN and entrapment labels for evaluation only.
Usage: precheck_scorer36.py <features.tsv.gz> <tag>"""
import os, sys, sqlite3, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score
SRC, TAG = sys.argv[1:3]
G = '/home/singjc/.jcode/scratch/dia2/gap7'
L = os.environ.get('PQP', '/home/singjc/.jcode/scratch/dia2/entlib10/cut/hela_r5.predicted_decoys.pqp')
DN = '/home/singjc/.jcode/scratch/dia2/entlib12/D/report.parquet'
PQ = f'{G}/{TAG}_feat.parquet'
if not os.path.exists(PQ):
    F = pd.read_csv(SRC, sep='\t', engine='c', dtype={'candidate_row': np.int32, 'precursor_id': np.int32, 'is_decoy': np.int8,
                    'fold': np.int8, 'rt_apex': np.float64, 'q_value': np.float64}, low_memory=False)
    for c in F.columns:
        if F[c].dtype == np.float64 and c not in ('rt_apex', 'q_value'): F[c] = F[c].astype(np.float32)
    F.to_parquet(PQ)
F = pd.read_parquet(PQ)
meta = ['candidate_row', 'precursor_id', 'is_decoy', 'fold', 'rt_apex', 'score', 'q_value']
X = [c for c in F.columns if c not in meta]
print(f'{TAG}: windows {len(F)}, features {len(X)}', flush=True)
c = sqlite3.connect(f'file:{L}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID pid, p.TRAML_ID t, p.DECOY d, p.CHARGE z, pe.MODIFIED_SEQUENCE ms,
  min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
pq['key'] = np.where(pq.d == 1, pq.t.str.replace('^DECOY_', '', regex=True), pq.t)
tg = pq[pq.d == 0]; r = (tg.e == 1).sum() / (tg.e == 0).sum(); k = 1 + 1 / r
pair_id = pd.Series(pd.factorize(pq.key)[0], index=pq.pid); ent = pq.set_index('pid').e
d = pd.read_parquet(DN, columns=['Modified.Sequence', 'Precursor.Charge', 'Q.Value', 'RT', 'Precursor.Quantity'])
d = d[d['Q.Value'] <= 0.01]
tgs = tg.assign(sk=tg.ms + '/' + tg.z.astype(str)).drop_duplicates('sk').set_index('sk')
d['pid'] = tgs.pid.reindex(d['Modified.Sequence'] + '/' + d['Precursor.Charge'].astype(str)).values
d = d.dropna(subset=['pid']); d['pid'] = d.pid.astype(int); d = d[ent.reindex(d.pid).values == 0].drop_duplicates('pid').set_index('pid')
pid = F.precursor_id.values; dec = F.is_decoy.values.astype(bool); fold = F.fold.values; pair = pair_id.reindex(pid).values

def best_rows(score):
    o = np.lexsort((-score, pid)); first = np.r_[True, pid[o][1:] != pid[o][:-1]]; rows = o[first]
    b = pd.DataFrame({'row': rows, 'pid': pid[rows], 's': score[rows], 'dec': dec[rows], 'pair': pair[rows]})
    b = b.sort_values(['pair', 's', 'dec'], ascending=[True, False, False]).drop_duplicates('pair')
    b = b.sort_values(['s', 'dec'], ascending=[False, False])
    nd = np.cumsum(b.dec.values); nt = np.cumsum(~b.dec.values)
    b['q'] = np.minimum.accumulate((nd / np.maximum(nt, 1))[::-1])[::-1]
    return b

def evaluate(score, name):
    b = best_rows(score); t = b[~b.dec.values].copy(); t['e'] = ent.reindex(t.pid).values
    acc = t[t.q <= 0.01]; N, E = len(acc), int(acc.e.sum())
    fdp = np.cumsum(t.e.values) * k / np.arange(1, len(t) + 1); ok = np.where(fdp <= 0.01)[0]
    n1 = int(ok.max() + 1) if len(ok) else 0
    print(f'{name:46s} q<=0.01 N {N:6d} lower {E / max(N, 1):.4f} combined {E * k / max(N, 1):.4f} | N@FDP_c 1 % {n1:6d} | in DIA-NN {acc.pid.isin(d.index).mean():.3f}', flush=True)
    return n1

Xv = F[X].values.astype(np.float32)
def gbm(pos_mask, neg_mask, seed=0, rounds=300, leaves=31, max_neg=1_500_000):
    s = np.zeros(len(F)); rng = np.random.default_rng(seed)
    for f in range(3):
        tp = np.where(pos_mask & (fold != f))[0]; tn = np.where(neg_mask & (fold != f))[0]
        if len(tn) > max_neg: tn = rng.choice(tn, max_neg, replace=False)
        idx = np.r_[tp, tn]; y = np.r_[np.ones(len(tp)), np.zeros(len(tn))]
        w = np.r_[np.full(len(tp), len(tn) / max(len(tp), 1)), np.ones(len(tn))]
        m = HistGradientBoostingClassifier(max_iter=rounds, learning_rate=0.05, max_leaf_nodes=leaves, random_state=seed).fit(Xv[idx], y, sample_weight=w)
        te = np.where(fold == f)[0]; s[te] = m.decision_function(Xv[te])
    return s

def winners(score):
    b = best_rows(score); m = np.zeros(len(F), bool); m[b.row.values[(b.q <= 0.01) & ~b.dec.values]] = True; return m

kept = F.score.values > -1e29
s_tool = np.where(kept, F.score.values, -1e30).astype(np.float64)
a = evaluate(s_tool, 'A tool score')
s = s_tool
for it in range(3):
    pos = winners(s); s = gbm(pos, dec, seed=it)
    bn = evaluate(s, f'B label-free GBM round {it + 1} ({pos.sum()} pos)')
sB = s
s = s_tool
for it in range(3):
    pos = winners(s); s = gbm(pos, dec, seed=10 + it, rounds=600, leaves=63)
b2 = evaluate(s, 'B2 label-free GBM 2x capacity (3 rounds)')
drt = np.abs(F.rt_apex.values - (d.RT * 60).reindex(pid).values)
pos_c = (~dec) & (drt <= 6)
sC = gbm(pos_c, dec, seed=9)
cn = evaluate(sC, f'C DIA-NN-labelled GBM ceiling ({pos_c.sum()} pos)')

# Per-quartile AUC of DIA-NN-confirmed best windows against decoy precursor best windows.
P = pd.DataFrame({'row': np.where((~dec) & (drt <= 6))[0]}); P['pid'] = pid[P.row]; P['s'] = s_tool[P.row]
P = P.sort_values('s', ascending=False).drop_duplicates('pid')
P['qq'] = pd.qcut(np.log10(d['Precursor.Quantity'].reindex(P.pid).clip(lower=1).values), 4, labels=['Q1', 'Q2', 'Q3', 'Q4'])
Dd = best_rows(s_tool); Dd = Dd[Dd.dec.values]
rng = np.random.default_rng(0); Dd = Dd.iloc[rng.choice(len(Dd), min(200000, len(Dd)), replace=False)]
for name, sc in (('tool', s_tool), ('B', sB), ('C', sC)):
    print(f'AUC {name:5s}', ' '.join(f'{q} {roc_auc_score(np.r_[np.ones((P.qq == q).sum()), np.zeros(len(Dd))], np.r_[sc[P.row[P.qq == q]], sc[Dd.row]]):.3f}' for q in ['Q1', 'Q2', 'Q3', 'Q4']))
print(f'\nSUMMARY {TAG}: tool {a}, label-free GBM {bn} ({bn / a - 1:+.1%}), 2x GBM {b2} ({b2 / a - 1:+.1%}), DIA-NN-labelled ceiling {cn} ({cn / a - 1:+.1%})')
