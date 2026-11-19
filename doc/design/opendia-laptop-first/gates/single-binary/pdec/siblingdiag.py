"""Sibling check (diagnostic): are 'fully explained' accepted winners explained by a co-eluting sibling, meaning the same
unmodified sequence at another charge or modification state? A sibling shares b/y fragments by construction, so the
interference graph treats it as a claimant. Per accepted winner (q<=0.01): has_sib = another accepted precursor of the same
UNMODIFIED_SEQUENCE in the same run with |rt_apex diff| <= 10 s."""
import sqlite3, glob, numpy as np, pandas as pd
c = sqlite3.connect('file:hela_r5.predicted_decoys.pqp?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID precursor_id, p.TRAML_ID t, p.DECOY d, pe.UNMODIFIED_SEQUENCE s,
 min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
 from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
 join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c).set_index('precursor_id')
pq['key'] = np.where(pq.d == 1, pq.t.str.replace('^DECOY_', '', regex=True), pq.t)
pq['pe'] = pq.key.map(pq[pq.d == 0].set_index('key').e)
rows = []
for f in sorted(glob.glob('runs/[DH]_cur_s0/*.features.tsv.gz')):
    F = pd.read_csv(f, sep='\t', usecols=['precursor_id', 'rt_apex', 'score', 'q_value', 'i_explained'])
    F = F[F.score > -1e29].sort_values('score', ascending=False).drop_duplicates('precursor_id')
    F = F[F.q_value <= 0.01].join(pq[['d', 'e', 'pe', 's']], on='precursor_id').reset_index(drop=True)
    F['run'] = f.split('/')[1][0]
    has = np.zeros(len(F), bool)
    for _, idx in F.groupby('s').groups.items():
        if len(idx) < 2: continue
        rt = F.loc[idx, 'rt_apex'].values
        has[np.asarray(idx)] = (np.abs(rt[:, None] - rt[None, :]) <= 10).sum(1) > 1
    F['has_sib'] = has; rows.append(F)
W = pd.concat(rows, ignore_index=True); W['full'] = W.i_explained >= 0.999
cls = {'entrapment targets': W[(W.d == 0) & (W.e == 1)], 'decoys (entrapment pairs)': W[(W.d == 1) & (W.pe == 1)],
       'HeLa targets': W[(W.d == 0) & (W.e == 0)]}
for k, v in cls.items():
    full = v.full
    print(f"{k:26s} n {len(v):6d}  fully explained {full.mean():.3f} | with co-eluting sibling: all {v.has_sib.mean():.3f}, "
          f"fully explained {v.has_sib[full].mean():.3f}, not fully explained {v.has_sib[~full].mean():.3f}")
H = cls['HeLa targets']
print(f"HeLa fully explained without a co-eluting sibling: {int((H.full & ~H.has_sib).sum())} of {int(H.full.sum())}")
