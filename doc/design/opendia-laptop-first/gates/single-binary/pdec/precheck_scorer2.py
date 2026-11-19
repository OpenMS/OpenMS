"""PRECHECK_SCORER part 2: full table from saved scores, plus per-feature separation of DIA-NN-confirmed windows vs
decoy windows by DIA-NN quantity quartile (diagnostic; DIA-NN and entrapment labels for evaluation only)."""
import sqlite3, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
G = '/home/singjc/.jcode/scratch/dia2/gap7'
L = '/home/singjc/.jcode/scratch/dia2/entlib10/cut/hela_r5.predicted_decoys.pqp'
DN = '/home/singjc/.jcode/scratch/dia2/entlib12/D/report.parquet'
F = pd.read_parquet(f'{G}/D_feat.parquet')
meta = ['candidate_row', 'precursor_id', 'is_decoy', 'fold', 'rt_apex', 'score', 'q_value']
X = [c for c in F.columns if c not in meta]
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
tg = tg.assign(sk=tg.ms + '/' + tg.z.astype(str)).drop_duplicates('sk').set_index('sk')
d['pid'] = tg.pid.reindex(d['Modified.Sequence'] + '/' + d['Precursor.Charge'].astype(str)).values
d = d.dropna(subset=['pid']); d['pid'] = d.pid.astype(int); d = d[ent.reindex(d.pid).values == 0].drop_duplicates('pid').set_index('pid')
pid = F.precursor_id.values; dec = F.is_decoy.values.astype(bool); pair = pair_id.reindex(pid).values

def evaluate(score, name):
    o = np.lexsort((-score, pid)); first = np.r_[True, pid[o][1:] != pid[o][:-1]]
    best = pd.DataFrame({'pid': pid[o][first], 's': score[o][first], 'dec': dec[o][first], 'pair': pair[o][first]})
    best = best.sort_values(['pair', 's', 'dec'], ascending=[True, False, False]).drop_duplicates('pair')
    best = best.sort_values(['s', 'dec'], ascending=[False, False])
    nd = np.cumsum(best.dec.values); nt = np.cumsum(~best.dec.values)
    q = np.minimum.accumulate((nd / np.maximum(nt, 1))[::-1])[::-1]
    t = best[~best.dec.values].assign(q=q[~best.dec.values]); t['e'] = ent.reindex(t.pid).values
    acc = t[t.q <= 0.01]; N, E = len(acc), int(acc.e.sum())
    fdp = np.cumsum(t.e.values) * k / np.arange(1, len(t) + 1); ok = np.where(fdp <= 0.01)[0]
    print(f'{name:36s} q<=0.01 N {N:6d} lower {E / max(N, 1):.4f} combined {E * k / max(N, 1):.4f} | '
          f'N@FDP_c 1 % {int(ok.max() + 1) if len(ok) else 0:6d} | in DIA-NN {acc.pid.isin(d.index).mean():.3f}')

kept = F.score.values > -1e29
evaluate(np.where(kept, F.score.values, -1e30).astype(np.float64), 'A tool score')
evaluate(np.load(f'{G}/score_B.npy'), 'B label-free GBM (3 rounds)')
sc = np.load(f'{G}/score_C.npy'); evaluate(sc, 'C DIA-NN-labelled GBM ceiling')

# Separation: DIA-NN-confirmed target windows (within 6 s of DIA-NN RT) vs decoy windows, by quantity quartile.
drt = np.abs(F.rt_apex.values - (d.RT * 60).reindex(pid).values)
posw = (~dec) & (drt <= 6)
P = F[posw].assign(qty=d['Precursor.Quantity'].reindex(F.precursor_id[posw]).values)
P = P.sort_values('score', ascending=False).drop_duplicates('precursor_id')
P['qq'] = pd.qcut(np.log10(P.qty.clip(lower=1)), 4, labels=['Q1', 'Q2', 'Q3', 'Q4'])
# Decoy reference: best decoy window per decoy precursor (what a target competes against).
Dd = F[dec & kept].sort_values('score', ascending=False).drop_duplicates('precursor_id')
rng = np.random.default_rng(0); Dd = Dd.iloc[rng.choice(len(Dd), min(200000, len(Dd)), replace=False)]
print(f'\nDIA-NN-confirmed best windows {len(P)}; decoy precursor best windows (sample) {len(Dd)}')
print('AUC vs decoys by DIA-NN quantity quartile (0.5 = no separation)')
cols = ['score'] + X
rows = []
for col in cols:
    row = {'feature': col}
    for q in ['Q1', 'Q2', 'Q3', 'Q4']:
        a = P.loc[P.qq == q, col].fillna(0).values; b = Dd[col].fillna(0).values
        row[q] = roc_auc_score(np.r_[np.ones(len(a)), np.zeros(len(b))], np.r_[a, b])
    rows.append(row)
print(pd.DataFrame(rows).set_index('feature').round(3).to_string())
print('\nDIA-NN-confirmed best window present per quartile (of DIA-NN q<=0.01 non-entrapment):')
d['qq'] = pd.qcut(np.log10(d['Precursor.Quantity'].clip(lower=1)), 4, labels=['Q1', 'Q2', 'Q3', 'Q4'])
print(d.assign(has=d.index.isin(P.precursor_id)).groupby('qq', observed=True).has.mean().round(3).to_dict())
