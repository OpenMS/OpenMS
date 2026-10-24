"""Prototype (#10693, lever 3): seeds at predicted positions for every reported precursor in every run where it is not
identified, plus one decoy probe per target at the same m/z and a fixed RT shift whose sign follows from the precursor
key (sum of its characters even: +, odd: -), as PIP-ECHO uses it with PIPECHO_PROBE_DECOYS=<shift> (prototype v2): the
shift is max(<shift>, 4 x the local half-window of the donor/acceptor pair), median over the donors.

Positions come from a PIPECHO_DIAG dump of the same data (local RT model centre per donor and acceptor run, aligned RT
scale; median over the donors). They are mapped to each run's own RT with the run's identified features: donor rows of
that run (aligned RT) joined on m/z (rounded to 4 decimals) and charge to the checkpoint features (own RT). Seed m/z is the precursor's
theoretical m/z. Each seed carries seed_apex_rt_tolerance (FFID drops seed features whose apex is further away).
Usage: make_probe_seeds.py link.tsv cand_pairs.tsv.gz feat_dir out_dir shift_s apex_tol_s"""
import sys, os, glob
import numpy as np, pandas as pd, pyarrow.parquet as pq, pyopenms as oms
link, pairs, feat_dir, out_dir, shift, tol = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], float(sys.argv[5]), float(sys.argv[6])
os.makedirs(out_dir, exist_ok=True)
d = pd.read_csv(link, sep='\t', usecols=['donor_map', 'acceptor_map', 'key', 'charge', 'mz', 'donor_rt', 'center', 'gate'])
runs = sorted(glob.glob(os.path.join(feat_dir, '*.featureParquet')))          # design order = map index
g = pd.read_csv(pairs, sep='\t', usecols=['key', 'acceptor_map'])               # reported precursors x runs without ID
sgn = d.key.map(lambda k: 1.0 if sum(map(ord, k)) % 2 == 0 else -1.0)
d['decoy_center'] = d.center + sgn * np.maximum(shift, 4.0 * d.gate)
cent = d.groupby(['key', 'acceptor_map']).center.median()
dcent = d.groupby(['key', 'acceptor_map']).decoy_center.median()
for m, path in enumerate(runs):
    name = os.path.basename(path).split('_fg')[0]
    f = pq.read_table(os.path.join(path, 'features.parquet'), columns=['rt', 'mz', 'charge', 'depth']).to_pandas()
    f = f[f.depth == 0]
    don = d[d.donor_map == m].drop_duplicates(['key'])[['mz', 'charge', 'donor_rt']]
    don = don.assign(mzr=don.mz.round(4)); f = f.assign(mzr=f.mz.round(4))
    j = don.merge(f[['mzr', 'charge', 'rt']], on=['mzr', 'charge']).drop_duplicates('donor_rt')
    j = j[(j.rt - j.donor_rt).abs() < 300].sort_values('donor_rt')       # drop coincidental m/z matches
    bins = np.array_split(np.arange(len(j)), 200)
    xa = np.array([j.donor_rt.iloc[b].median() for b in bins]); xr = np.array([j.rt.iloc[b].median() for b in bins])
    to_raw = lambda t: np.interp(t, xa, xr)
    sel = g[g.acceptor_map == m]
    fm = oms.FeatureMap()
    n = 0
    for key in sel.key:
        c = cent.get((key, m))
        if c is None or np.isnan(c): continue
        seq, z = key.rsplit('_', 1); z = int(z)
        mz = oms.AASequence.fromString(seq).getMZ(z)
        for rt in (c, dcent.get((key, m))):
            ft = oms.Feature(); ft.setMZ(mz); ft.setRT(float(to_raw(rt))); ft.setCharge(z)
            ft.setMetaValue('seed_apex_rt_tolerance', tol); fm.push_back(ft)
        n += 1
    fm.ensureUniqueId()
    oms.FeatureXMLFile().store(os.path.join(out_dir, name + '.featureXML'), fm)
    print(name, 'shared IDs for the RT map', len(j), 'targets', n, 'seeds', fm.size(), flush=True)
