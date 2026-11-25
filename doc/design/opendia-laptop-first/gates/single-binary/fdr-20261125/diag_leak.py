#!/usr/bin/env python3
"""Cross-fold peak-sharing test. Cross-fitting assumes no row of the scored fold resembles a training row. A null
precursor whose kept window sits on the peak of a confident target (same apex RT, precursor in the same isolation
window) inherits that peak's MS2/shape features; if that target is in another fold it was a GBM training positive.
For tail decoys and tail entrapments (matched by decoy rank under each score) and for accepted entrapments: the share
whose apex is within 1.5 s of a network-confident target with |dmz| < 2 Th, and of those the share in another fold."""
import sqlite3, sys, numpy as np, pandas as pd
lib = sys.argv[1]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID precursor_id, p.PRECURSOR_MZ mz,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
P = pq.set_index('precursor_id')
for path in sys.argv[2:]:
    f = pd.read_csv(path, sep='\t', usecols=['precursor_id', 'is_decoy', 'fold', 'rt_apex', 'score', 'q_value', 'network_score', 'network_q_value'])
    nk = f[f.network_score > -1e29].set_index('precursor_id').join(P)
    ids = nk[(nk.is_decoy == 0) & (nk.e == 0) & (nk.network_q_value <= 0.01)].sort_values('rt_apex')
    irt, imz, ifo = ids.rt_apex.values, ids.mz.values, ids.fold.values
    def share(rows):
        hit = oth = 0
        for rt, mz, fo in zip(rows.rt_apex.values, rows.mz.values, rows.fold.values):
            lo, hi = np.searchsorted(irt, rt - 1.5), np.searchsorted(irt, rt + 1.5)
            m = np.abs(imz[lo:hi] - mz) < 2
            if m.any():
                hit += 1
                oth += bool((ifo[lo:hi][m] != fo).any())
        return hit, oth
    print(f'== {path.split("/")[0]} (confident targets {len(ids)})')
    for s, q in (('network_score', 'network_q_value'), ('score', 'q_value')):
        k = f[f[s] > -1e29].set_index('precursor_id').join(P)
        d = k[k.is_decoy == 1].sort_values(s, ascending=False)
        t = d[s].values[199]
        topD = d.head(200)
        topE = k[(k.is_decoy == 0) & (k.e == 1) & (k[s] >= t)]
        accE = k[(k.is_decoy == 0) & (k.e == 1) & (k[q] <= 0.01)]
        bulkE = k[(k.is_decoy == 0) & (k.e == 1)].sample(3000, random_state=0)
        bulkD = d.sample(3000, random_state=0)
        out = []
        for lab, r in (('top200 D', topD), ('E above D200', topE), ('accepted E', accE), ('random E', bulkE), ('random D', bulkD)):
            h, o = share(r)
            out.append(f'{lab} {h}/{len(r)} ({h / max(1, len(r)):.3f}, other fold {o})')
        print(f'  {s:14s}: ' + '; '.join(out))
