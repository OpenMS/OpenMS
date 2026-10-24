"""Prototype check: features FFID extracted at the target seeds vs the decoy-probe seeds of one run. Seeds come in
(target, decoy) pairs in the seed file. A seed counts as extracted if the run's checkpoint has an unidentified
feature of that charge within 10 ppm and inside the seed's apex tolerance. Usage: seed_features.py seeds.featureXML run.featureParquet"""
import sys, numpy as np, pandas as pd, pyarrow.parquet as pq, pyopenms as oms
fm = oms.FeatureMap(); oms.FeatureXMLFile().load(sys.argv[1], fm)
seeds = pd.DataFrame([(f.getMZ(), f.getRT(), f.getCharge(), float(f.getMetaValue('seed_apex_rt_tolerance'))) for f in fm], columns=['mz', 'rt', 'z', 'tol'])
seeds['kind'] = np.where(np.arange(len(seeds)) % 2 == 0, 'target', 'decoy')
d = sys.argv[2]
f = pq.read_table(d + '/features.parquet', columns=['unique_id', 'rt', 'mz', 'charge', 'intensity', 'depth', 'metavalues', 'quality']).to_pandas()
f = f[f.depth == 0]
p = pq.read_table(d + '/psms.parquet', columns=['feature_unique_id']).to_pandas()
# float-stored ids: compare on rounded values
ids = set(np.round(p.feature_unique_id.values / 1e6))
f['noid'] = ~np.isin(np.round(f.unique_id.values.astype(float) / 1e6), list(ids))
f = f[f.noid].sort_values('mz')
def mv(row, k):
    for m in row:
        if m['name'] == k: return m['value']
    return None
f['masserr'] = [mv(r, 'masserror_ppm') for r in f.metavalues]
M, R, Z = f.mz.values, f.rt.values, f.charge.values
res = []
for s in seeds.itertuples():
    lo, hi = np.searchsorted(M, s.mz * (1 - 10e-6)), np.searchsorted(M, s.mz * (1 + 10e-6), 'right')
    best = None
    for i in range(lo, hi):
        if Z[i] == s.z and abs(R[i] - s.rt) <= s.tol:
            if best is None or f.intensity.values[i] > f.intensity.values[best]: best = i
    res.append((s.kind, best is not None, f.intensity.values[best] if best is not None else np.nan,
                abs(R[best] - s.rt) if best is not None else np.nan, f.quality.values[best] if best is not None else np.nan))
x = pd.DataFrame(res, columns=['kind', 'extracted', 'intensity', 'apex_dt', 'quality'])
print(x.groupby('kind').agg(n=('extracted', 'size'), extracted_pct=('extracted', lambda s: round(100 * s.mean(), 1)),
      median_intensity=('intensity', 'median'), median_apex_dt=('apex_dt', 'median'), median_quality=('quality', 'median')).to_string())
e = x[x.extracted]
for k in ('target', 'decoy'):
    print(k, 'intensity quartiles', np.nanpercentile(e[e.kind == k].intensity, [25, 50, 75]).round(0))
