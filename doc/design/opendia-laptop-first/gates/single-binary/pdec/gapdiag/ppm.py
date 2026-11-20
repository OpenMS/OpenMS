import pandas as pd, numpy as np
d = pd.read_parquet('/nfs/research/juan/singj/opendia/entlib-diann5-20261119/D/report.parquet')
d = d[(d['Q.Value'] <= 0.01) & (d.Decoy == 0)].set_index('Precursor.Id')
cls = pd.read_csv('classes_D.tsv', sep='\t').drop_duplicates('Precursor.Id').set_index('Precursor.Id').cls
d = d.join(cls, how='inner')
d['fr_ppm'] = d['Best.Fr.Mz.Delta'] / d['Best.Fr.Mz'] * 1e6
d['ms1_ppm'] = d['Ms1.Apex.Mz.Delta'] / d['Precursor.Mz'] * 1e6
off, hw = 0.647, 8.338
print('DIA-NN best-fragment and MS1-apex mass errors by class; ours: MS2 offset 0.65 ppm, half width 8.34 ppm')
for k in ['A', 'B1', 'B2', 'C']:
    g = d[d.cls == k]; e = (g.fr_ppm - off).abs(); m = g.ms1_ppm.abs()
    print(f'{k:3s} n {len(g):6d} | MS2 |err-offset| median {e.median():5.2f} p90 {e.quantile(0.9):5.2f} p99 {e.quantile(0.99):6.2f} outside +-{hw:.2f}: {(e > hw).mean():.3f}'
          f' | MS1 |err| median {m.median():5.2f} p90 {m.quantile(0.9):5.2f} outside 20 ppm: {(m > 20).mean():.3f}')
print('signed MS2 error median by class:', {k: round(float(d.loc[d.cls == k, "fr_ppm"].median()), 2) for k in ['A', 'B1', 'B2', 'C']})
