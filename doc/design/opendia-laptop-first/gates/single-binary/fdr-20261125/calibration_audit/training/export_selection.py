#!/usr/bin/env python3
"""Export full network precursor winners for selection diagnostics. No PROTEIN table or entrapment label is read.
Usage: export_selection.py library.pqp precursor.tsv out.tsv
"""
import sqlite3,sys,pandas as pd,numpy as np
lib,table,out=sys.argv[1:]
c=sqlite3.connect(f'file:{lib}?mode=ro',uri=True)
p=pd.read_sql('select p.ID precursor_id,p.TRAML_ID name,p.DECOY decoy,pe.UNMODIFIED_SEQUENCE seq from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID',c)
assert not p.precursor_id.duplicated().any()
p['key']=np.where(p.decoy==1,p.name.str.replace('^DECOY_','',regex=True),p.name)
seq=p[p.decoy==0].set_index('key').seq
p['family']=p.key.map(seq).fillna(p.seq)
def fnv(s):
    h=14695981039346656037
    for b in s.encode():h=((h^b)*1099511628211)&((1<<64)-1)
    return h
p['fold']=p.family.map(lambda s:fnv(s)%3)
p['pair']=pd.factorize(p.key)[0]
p['dense']=np.arange(len(p))
t=pd.read_csv(table,sep='\t',usecols=['precursor_id','score','is_decoy'])
t=t[t.score>-1e29].sort_values('score',ascending=False).drop_duplicates('precursor_id').merge(p,on='precursor_id',validate='one_to_one')
assert (t.is_decoy==t.decoy).all()
t[['dense','pair','fold','decoy','score']].to_csv(out,sep='\t',header=False,index=False)
print(f'{table}: {len(t)} rows, fold counts {t.fold.value_counts().to_dict()}')
