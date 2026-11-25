#!/usr/bin/env python3
"""Summary table per run and pooled: own q<=0.01 for network and GBM scores.
Full entrapment set: N, E, lower N_E/N, combined 2 N_E/N (r = 1), N@FDP_c1%.
Parent-sharing split: entrapments whose library fragments match >= 2 of their isobaric human parent's (ov2+, 7.4 % of
entrapments) vs the rest; restricted bound excludes ov2+ entrapments and uses r' = r * share(ov<=1). Post hoc, descriptive only."""
import sqlite3, sys, numpy as np, pandas as pd
lib = sys.argv[1]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID precursor_id, min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
E = pq.set_index('precursor_id').e
ov = pd.read_csv('overlap.tsv', sep='\t', index_col=0).ov
ent_ids = E[E == 1].index
share_lo = float((ov.reindex(ent_ids).fillna(0) <= 1).mean())
r_lo = 1.0 * share_lo
print(f'share of entrapments with ov<=1: {share_lo:.4f}, r\' = {r_lo:.4f}')
rows = []
for path in sys.argv[2:]:
    f = pd.read_csv(path, sep='\t', usecols=['precursor_id', 'is_decoy', 'score', 'q_value', 'network_score', 'network_q_value'])
    f['e'] = f.precursor_id.map(E).fillna(0).astype(int)
    f['hi'] = (f.e == 1) & (f.precursor_id.map(ov).fillna(0) >= 2)
    for s, q in (('net', 'network_q_value'), ('gbm', 'q_value')):
        sc = 'network_score' if s == 'net' else 'score'
        t = f[(f[sc] > -1e29) & (f.is_decoy == 0)].sort_values([q, sc], ascending=[True, False])
        acc = t[t[q] <= 0.01]
        fdp = np.cumsum(t.e.values) * 2 / np.arange(1, len(t) + 1)
        ok = np.where(fdp <= 0.01)[0]
        lo = t[~t.hi]
        fdp2 = np.cumsum(lo.e.values) * (1 + 1 / r_lo) / np.arange(1, len(lo) + 1)
        ok2 = np.where(fdp2 <= 0.01)[0]
        acc2 = acc[~acc.hi]
        rows.append(dict(run=path.split('/')[0], score=s, N=len(acc), E=int(acc.e.sum()), E_hi=int(acc.hi.sum()),
                         N1=int(ok.max() + 1) if len(ok) else 0, N_lo=len(acc2), E_lo=int(acc2.e.sum()), N1_lo=int(ok2.max() + 1) if len(ok2) else 0))
R = pd.DataFrame(rows)
R['lower'] = R.E / R.N; R['combined'] = 2 * R.E / R.N; R['combined_lo'] = R.E_lo * (1 + 1 / r_lo) / R.N_lo
pd.set_option('display.width', 200)
print(R.to_string(index=False, float_format=lambda v: f'{v:.4f}'))
print('\npooled:')
for s, g in R.groupby('score'):
    N, Ee, Eh, Nl, El = g.N.sum(), g.E.sum(), g.E_hi.sum(), g.N_lo.sum(), g.E_lo.sum()
    print(f'  {s}: runs {len(g)} N {N} E {Ee} (ov2+ {Eh}) lower {Ee / N:.4f} combined {2 * Ee / N:.4f} | without ov2+: combined {El * (1 + 1 / r_lo) / Nl:.4f} | '
          f'median N@FDP_c1% {g.N1.median():.0f}, runs over 1 % combined {int((g.combined > 0.01).sum())}, worst {g.combined.max():.4f}')
