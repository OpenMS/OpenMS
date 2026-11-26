#!/usr/bin/env python3
"""SEP (diagnosis): in which features do accepted-region entrapments differ from decoys that score the same?
Rows: kept window per precursor, top region (q <= 0.3). Within score deciles of decoys, compare the feature distributions of
entrapment targets vs decoys (both are false). A feature with a strong, consistent shift (standardised mean difference
within score bins) is a property that lets false targets beat their decoys at equal score, i.e. a non-exchangeability
source. Human targets shown for reference. Usage: sep.py <lib.pqp> <reduced features.tsv.gz> ...   env RELABEL"""
import sqlite3, sys, os, numpy as np, pandas as pd
lib = sys.argv[1]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
P = pd.read_sql("""select p.ID id, pe.UNMODIFIED_SEQUENCE useq,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID where p.DECOY=0 group by p.ID""", c)
if os.environ.get('RELABEL'):
    rl = set(pd.read_csv(os.environ['RELABEL'], sep='\t').seq); P.loc[P.useq.isin(rl), 'e'] = 0
ent = P.set_index('id').e
X = pd.concat([pd.read_csv(f, sep='\t') for f in sys.argv[2:]])
X = X.sort_values('score', ascending=False).drop_duplicates(['precursor_id']) if X.precursor_id.duplicated().any() and len(sys.argv) == 3 else X
X['cls'] = np.where(X.is_decoy == 1, 'D', np.where(X.precursor_id.map(ent).fillna(0) == 1, 'E', 'H'))
skip = {'candidate_row', 'precursor_id', 'is_decoy', 'fold', 'rt_apex', 'score', 'q_value', 'cls', 'network_score', 'network_q_value'}
feats = [f for f in X.columns if f not in skip]
for band, lo, hi in (('q<=0.01', 0, 0.01), ('0.01-0.05', 0.01, 0.05), ('0.05-0.3', 0.05, 0.3)):
    Z = X[(X.q_value > lo) & (X.q_value <= hi)]
    if (Z.cls == 'E').sum() < 30: continue
    dec = Z[Z.cls == 'D'].score
    edges = np.unique(np.quantile(dec, np.linspace(0, 1, 6))) if len(dec) > 30 else np.array([Z.score.min(), Z.score.max()])
    Z = Z.assign(b=np.clip(np.searchsorted(edges, Z.score, side='right') - 1, 0, len(edges) - 2))
    rows = []
    for f in feats:
        v = Z[f].replace([np.inf, -np.inf], np.nan)
        sd = v[Z.cls == 'D'].std()
        if not (sd > 0): continue
        d_ed, d_hd, w = 0.0, 0.0, 0
        for b, g in Z.assign(v=v).groupby('b'):
            ne, nd = (g.cls == 'E').sum(), (g.cls == 'D').sum()
            if ne < 5 or nd < 5: continue
            d_ed += ne * (g.v[g.cls == 'E'].mean() - g.v[g.cls == 'D'].mean()) / sd
            d_hd += ne * (g.v[g.cls == 'H'].mean() - g.v[g.cls == 'D'].mean()) / sd if (g.cls == 'H').sum() >= 5 else 0
            w += ne
        if w: rows.append((f, d_ed / w, d_hd / w))
    R = pd.DataFrame(rows, columns=['feature', 'E_minus_D', 'H_minus_D']).assign(a=lambda r: r.E_minus_D.abs()).sort_values('a', ascending=False)
    print(f'== {band}: E {(Z.cls == "E").sum()} D {(Z.cls == "D").sum()} H {(Z.cls == "H").sum()} (score-binned SMD, in decoy SD units)')
    print(R.head(14).drop(columns='a').round(3).to_string(index=False))
