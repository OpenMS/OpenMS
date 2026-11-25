#!/usr/bin/env python3
"""Are the extra GBM entrapments riding on their parent peptide's signal? Diagnosis only.
Entrapment accession sp|<pep>_p_target|..., parent sp|<pep>_target|... For each accepted entrapment (own q<=0.01) under
network and GBM: is its parent a confident ID, in another fold (so a training positive of the model that scores the entrapment),
co-eluting (|apex diff| <= 3 s), and how similar is the sequence (positional identity)."""
import sqlite3, sys, numpy as np, pandas as pd
lib = sys.argv[1]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID precursor_id, p.DECOY d, p.CHARGE z, pe.UNMODIFIED_SEQUENCE seq, min(pr.PROTEIN_ACCESSION) acc
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
pq = pq[pq.d == 0].copy()
pq['ent'] = pq.acc.str.contains('_p_target')
pq['root'] = pq.acc.str.extract(r'\|([A-Z]+)_(?:p_)?target')[0]
P = pq.set_index('precursor_id')
def ident(a, b):
    return np.mean([x == y for x, y in zip(a, b)]) if len(a) == len(b) else np.nan
for path in sys.argv[2:]:
    f = pd.read_csv(path, sep='\t', usecols=['precursor_id', 'is_decoy', 'fold', 'rt_apex', 'score', 'q_value', 'network_score', 'network_q_value'])
    f = f[f.is_decoy == 0]
    for s, q in (('network_score', 'network_q_value'), ('score', 'q_value')):
        k = f[f[s] > -1e29].set_index('precursor_id')
        k = k.join(P[['ent', 'root', 'z', 'seq']])
        par = k[~k.ent].reset_index().set_index(['root', 'z'])
        acc = k[k.ent & (k[q] <= 0.01)]
        rows = []
        for pid, r in acc.iterrows():
            key = (r.root, r.z)
            if key in par.index:
                pr = par.loc[key]
                pr = pr.iloc[0] if isinstance(pr, pd.DataFrame) else pr
                rows.append(dict(par_id=pr[q] <= 0.01, other_fold=pr.fold != r.fold, coelute=abs(pr.rt_apex - r.rt_apex) <= 3, ident=ident(pr.seq, r.seq)))
            else:
                rows.append(dict(par_id=False, other_fold=np.nan, coelute=False, ident=np.nan))
        R = pd.DataFrame(rows)
        lk = R.par_id & R.coelute
        print(f'{path.split("/")[0]} {s:14s}: accepted entrapments {len(R)}; parent same charge is an ID {R.par_id.mean():.2f}; '
              f'ID and co-eluting {lk.mean():.2f} (of those other fold {R[lk].other_fold.mean():.2f}); median positional identity {R.ident.median():.2f}')
    # baseline: all entrapments whose parent is an ID, fraction co-eluting
    k = f[f.score > -1e29].set_index('precursor_id').join(P[['ent', 'root', 'z']])
    par = k[~k.ent & (k.q_value <= 0.01)].reset_index().drop_duplicates(['root', 'z']).set_index(['root', 'z'])
    e = k[k.ent].reset_index().join(par[['rt_apex']].rename(columns={'rt_apex': 'prt'}), on=['root', 'z'])
    e = e[e.prt.notna()]
    print(f'   all entrapments with an ID parent: {len(e)}, co-eluting with parent {np.mean(abs(e.prt - e.rt_apex) <= 3):.3f}')
