"""Per-run DIA-NN 1% precursors (full report) and how many are in our library."""
import pandas as pd
D='/home/singjc/Documents/datasets/proteobench-astral-entrapment/proteobench_diann_v231_res_20260717_123332/'
r=pd.read_parquet(D+'input_file.parquet',columns=['Run','Stripped.Sequence','Precursor.Charge','Decoy','Q.Value','RT','Protein.Ids','Precursor.Mz'])
r=r[(r.Decoy==0)&(r['Q.Value']<=0.01)].copy()
r['key']=r['Stripped.Sequence']+'/'+r['Precursor.Charge'].astype(str)
r['entrap']=r['Protein.Ids'].str.contains('_p_target',regex=False) & ~r['Protein.Ids'].str.contains(';',regex=False) if False else r['Protein.Ids'].apply(lambda s: all('_p_target' in p for p in s.split(';')))
t=pd.read_parquet('lib_targets.parquet'); keys=set(t.key)
print('unique DIA-NN 1% precursors over 3 runs', r.key.nunique())
for run,g in r.groupby('Run'):
    inl=g[g.key.isin(keys)]
    print(run,'DIA-NN 1%',len(g),'entrap',int(g.entrap.sum()),'| in our library',len(inl),'entrap',int(inl.entrap.sum()), 'charges', inl['Precursor.Charge'].value_counts().sort_index().to_dict())
r.to_parquet('diann_runs_1pct.parquet')
hs=t[~t.key.isin(set())]
print('library human-or-entrap target precursors', len(t), 'with charge 2-3', int(t.charge.isin([2,3]).sum()))
