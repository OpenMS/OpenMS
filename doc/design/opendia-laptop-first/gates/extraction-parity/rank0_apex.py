"""Peak choice independent of boundary width: is our rank-0 apex within 0 / 1 / 2 scans of DIA-NN's apex?
raw01, DIA-NN q <= 1%. Also entrapments among the matches (evaluation only). Usage: rank0_apex.py run [run ...]"""
import sys, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; X = B + '/extraction-parity-20261112'
idm = pd.read_csv(X + '/idmap.tsv', sep='\t')
rep = pd.read_parquet(X + '/diann/report_raw01.parquet', columns=['Precursor.Id', 'Q.Value', 'RT']).merge(idm[['Precursor.Id', 'precursor_id']], on='Precursor.Id')
ent = pd.read_csv(B + '/diann-ourlib-20261112/diann_ourlib_raw01.tsv.gz', sep='\t', usecols=['precursor_id', 'e']).drop_duplicates('precursor_id')
rep = rep.merge(ent, on='precursor_id'); rep = rep[rep['Q.Value'] <= 0.01].set_index('precursor_id'); rep['RT'] *= 60
out = open(X + '/RANK0_APEX_raw01.txt', 'a')
for run in sys.argv[1:]:
    f = pd.read_csv(X + '/' + run + '/search.features.tsv', sep='\t', usecols=['precursor_id', 'is_decoy', 'rt_apex', 'peak_rank', 'rt_start', 'rt_end'])
    f = f[(f.is_decoy == 0) & f.precursor_id.isin(rep.index)]
    top = f[f.peak_rank == 0].set_index('precursor_id')
    d = ((top.rt_apex - rep.RT.reindex(top.index)).abs() / 1.646).round()
    anyp = f.assign(d=((f.rt_apex - rep.RT.reindex(f.precursor_id).values).abs() / 1.646).round()).groupby('precursor_id').d.min()
    e = rep.e.reindex(top.index)
    s = (f'{run}: rank-0 apex at DIA-NN apex 0 scans {int((d == 0).sum())} | <=1 {int((d <= 1).sum())} (entrapments {int(e[d <= 1].sum())}) | '
         f'<=2 {int((d <= 2).sum())} | any of 5 peaks <=1 {int((anyp <= 1).sum())} | width median s {float((top.rt_end - top.rt_start).median()):.2f}')
    print(s); out.write(s + '\n')
