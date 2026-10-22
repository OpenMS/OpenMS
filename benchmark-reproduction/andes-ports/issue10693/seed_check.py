"""For precursor x run pairs of one run without a PIP-ECHO candidate: was there a seed (FeatureFinderMultiplex,
ProteomicsLFQ -debug 667: debug_seeds_*.featureXML) of the same charge within 10 ppm inside the local RT window,
and did FeatureFinderIdentification keep a feature there (debug_fraction_*.featureXML)? The window centre is on the
aligned RT scale; it is mapped back to the run's own RT with the identified features both scales share.
Usage: seed_check.py link.tsv cand_pairs.tsv.gz acceptor_map seeds.featureXML ffid.featureXML"""
import sys, numpy as np, pandas as pd, pyopenms as oms
link, pairs, amap, seedf, ffidf = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5]
def load(fn):
    fm = oms.FeatureMap(); oms.FeatureXMLFile().load(fn, fm)
    rows = []
    for f in fm:
        ids = f.getPeptideIdentifications()
        key = None
        if ids and ids[0].getHits():
            key = f"{ids[0].getHits()[0].getSequence().toString()}_{f.getCharge()}"
        rows.append((f.getMZ(), f.getRT(), f.getCharge(), key, f.getIntensity()))
    return pd.DataFrame(rows, columns=['mz', 'rt', 'charge', 'key', 'intensity'])
seeds, ffid = load(seedf), load(ffidf)
print('seeds', len(seeds), '| FFID features', len(ffid), 'with ID', int(ffid.key.notna().sum()), 'without', int(ffid.key.isna().sum()))
# raw RT <-> aligned RT from the identified features of this run
d = pd.read_csv(link, sep='\t', usecols=['donor_map', 'key', 'donor_rt'])
al = d[d.donor_map == amap].groupby('key').donor_rt.median()
raw = ffid.dropna(subset=['key']).groupby('key').rt.median()
j = pd.concat([raw.rename('raw'), al.rename('al')], axis=1, join='inner').sort_values('al')
bins = np.array_split(np.arange(len(j)), 200)
xa = np.array([j.al.iloc[b].median() for b in bins]); xr = np.array([j.raw.iloc[b].median() for b in bins])
to_raw = lambda t: np.interp(t, xa, xr)
print('RT map from', len(j), 'shared IDs; median |raw-aligned|', round(float((j.raw - j.al).abs().median()), 1), 's')
g = pd.read_csv(pairs, sep='\t'); g = g[g.acceptor_map == amap].copy()
def index(df):
    df = df.sort_values('mz'); return df.mz.values, df.rt.values, df.charge.values, df.key.values
S, F = index(seeds), index(ffid)
def hit(ix, mz, rt, z, gate, need_noid=False):
    m, r, c, k = ix
    lo, hi = np.searchsorted(m, mz * (1 - 10e-6)), np.searchsorted(m, mz * (1 + 10e-6), 'right')
    for i in range(lo, hi):
        if c[i] == z and abs(r[i] - rt) <= gate and (not need_noid or pd.isna(k[i])): return True
    return False
theo = lambda key: oms.AASequence.fromString(key.rsplit('_', 1)[0]).getMZ(int(key.rsplit('_', 1)[1]))
g['mz_t'] = [theo(k) for k in g.key]
g['rt_raw'] = to_raw(g.rt.values)
g['seed'] = [hit(S, m, t, z, w) for m, t, z, w in zip(g.mz_t, g.rt_raw, g.charge, g.gate)]
g['ffid_noid'] = [hit(F, m, t, z, w, True) for m, t, z, w in zip(g.mz_t, g.rt_raw, g.charge, g.gate)]
g['ffid_any'] = [hit(F, m, t, z, w) for m, t, z, w in zip(g.mz_t, g.rt_raw, g.charge, g.gate)]
t = g.groupby('status')[['seed', 'ffid_noid', 'ffid_any']].mean().mul(100).round(1)
t['n'] = g.groupby('status').size()
print(f'run map {amap}: share of pairs with a seed / an FFID feature without ID / any FFID feature at the spot (%)')
print(t.to_string())
