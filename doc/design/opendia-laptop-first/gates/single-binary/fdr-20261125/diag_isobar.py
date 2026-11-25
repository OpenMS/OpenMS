#!/usr/bin/env python3
"""Isobaric-parent test. Entrapments are shuffles of human peptides keeping K/R and the C-terminus, so each has the
same precursor m/z as its parent. Decoys (reverse, terminal mutated) are mass-shifted and cannot share a present
target's MS1/fragment signal. For accepted entrapments (own q<=0.01): is an isobaric human target (same composition,
charge, mods) an ID co-eluting within 5 s? Compare with matched tail decoys (their isobaric human partner, if any). Diagnosis only."""
import sqlite3, sys, numpy as np, pandas as pd
lib = sys.argv[1]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID precursor_id, p.DECOY d, p.CHARGE z, p.PRECURSOR_MZ mz, pe.UNMODIFIED_SEQUENCE seq, pe.MODIFIED_SEQUENCE mseq,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
pq['nmod'] = pq.mseq.str.count('UniMod')
pq['comp'] = pq.seq.map(lambda s: ''.join(sorted(s))) + '|' + pq.z.astype(str) + '|' + pq.nmod.astype(str)
P = pq.set_index('precursor_id')
hum = pq[(pq.d == 0) & (pq.e == 0)]
for path in sys.argv[2:]:
    f = pd.read_csv(path, sep='\t', usecols=['precursor_id', 'is_decoy', 'rt_apex', 'score', 'q_value', 'network_score', 'network_q_value'])
    print(f'== {path.split("/")[0]}')
    for s, q in (('network_score', 'network_q_value'), ('score', 'q_value')):
        k = f[f[s] > -1e29].set_index('precursor_id').join(P[['comp', 'e', 'mz']])
        ids = k[(k.is_decoy == 0) & (k.e == 0) & (k[q] <= 0.01)]
        idx = ids.reset_index().groupby('comp').rt_apex.apply(list)
        def has_parent(r):
            rts = idx.get(r.comp)
            return rts is not None and min(abs(np.array(rts) - r.rt_apex)) <= 5
        acc_e = k[(k.e == 1) & (k.is_decoy == 0) & (k[q] <= 0.01)]
        d = k[k.is_decoy == 1].sort_values(s, ascending=False).head(200)
        ne = sum(has_parent(r) for _, r in acc_e.iterrows())
        isobar = acc_e.comp.isin(idx.index).mean() if len(acc_e) else np.nan
        nd = sum(has_parent(r) for _, r in d.iterrows())
        # background: all scored entrapments whose isobaric human parent is an ID; fraction accepted
        allE = k[(k.e == 1) & (k.is_decoy == 0)]
        withp = allE[allE.comp.isin(idx.index)]
        wp_co = withp[[has_parent(r) for _, r in withp.iterrows()]]
        print(f'  {s:14s}: accepted entrapments {len(acc_e)}, isobaric human parent is an ID {isobar:.2f}, parent ID co-eluting <=5 s {ne} ({ne / max(1, len(acc_e)):.2f}); '
              f'top-200 decoys with an isobaric co-eluting human ID {nd}')
        print(f'      entrapments with an ID parent {len(withp)}, co-eluting {len(wp_co)}; acceptance rate co-eluting {np.mean(wp_co[q] <= 0.01) if len(wp_co) else float("nan"):.3f}, '
              f'not co-eluting {np.mean(withp.drop(wp_co.index)[q] <= 0.01):.4f}, no ID parent {np.mean(allE.drop(withp.index)[q] <= 0.01):.4f}')
