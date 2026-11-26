#!/usr/bin/env python3
"""Cross-check reconstructed family folds against actual production feature exports, first 1,000 rows each."""
import sqlite3,pandas as pd
lib='/home/singjc/.jcode/scratch/blend-20261125/eval/hela_r5.predicted_decoys.pqp'
c=sqlite3.connect(f'file:{lib}?mode=ro',uri=True)
def fnv(s):
    h=14695981039346656037
    for b in s.encode(): h=((h^b)*1099511628211)&((1<<64)-1)
    return h
for s in ('D','H'):
    t=pd.read_csv(f'/home/singjc/.jcode/scratch/calsearch/top/{s}_G0f_s0.tsv.gz',sep='\t',usecols=['precursor_id','fold'],nrows=1000)
    ids=','.join(map(str,t.precursor_id.unique()))
    q=f'''select p.ID precursor_id,p.DECOY decoy,p.TRAML_ID name,pe.UNMODIFIED_SEQUENCE seq from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID where p.ID in ({ids})'''
    rows=pd.read_sql(q,c)
    checks=0
    keys=[r.name[6:] for r in rows.itertuples() if r.decoy]
    target_seq=dict(c.execute('select p.TRAML_ID,pe.UNMODIFIED_SEQUENCE from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID where p.DECOY=0 and p.TRAML_ID in ('+','.join('?' for _ in keys)+')',keys))
    for r in rows.itertuples():
        seq=target_seq.get(r.name[6:],r.seq) if r.decoy else r.seq
        assert (t.loc[t.precursor_id==r.precursor_id,'fold']==fnv(seq)%3).all()
        checks+=1
    print(f'{s}: actual production fold agrees for {checks} distinct precursors in first {len(t)} feature rows')
