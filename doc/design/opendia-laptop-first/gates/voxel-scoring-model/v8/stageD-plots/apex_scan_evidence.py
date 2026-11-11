"""Apex-scan evidence at native sampling (raw 01, old population, D11FWT). Descriptive; labels only define groups.
For each winner: the MS2 scans whose isolation window holds the precursor m/z, nearest our apex (k = 0) and the two
neighbours (k = -1, +1; revisit 1.65 s). Per scan: library fragments (6) found within 10 ppm of the calibrated m/z, and
their summed share of the scan's TIC. Also: points across the peak = scans in apex +- 3 s with >= 3 fragments present.
Groups: D (all stage-D misses), DEC and ENT winners with score in the D interquartile range, E (accepted, sample 600).
Question: at matched score, do D windows show more library fragments in the raw apex scans than decoys and entrapment?"""
import json, sqlite3, numpy as np, pandas as pd, pyopenms as oms
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; V = B + '/astral-voxel-pilot-20261110/v2'; OUT = B + '/dd-plots-20261111'
OFF = float(json.load(open(B + '/astral-design-raw01-20261107/window300k/cal.json'))['ms2_mz']['offset_ppm']); PPM = 10.0
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
pm = pd.read_sql('select ID precursor_id, PRECURSOR_MZ pmz from PRECURSOR', c).set_index('precursor_id').pmz
tr = pd.read_sql("""select m.PRECURSOR_ID precursor_id, t.PRODUCT_MZ mz from TRANSITION_PRECURSOR_MAPPING m join TRANSITION t on t.ID=m.TRANSITION_ID""", c)
fr = {p: g.mz.values * (1 + OFF * 1e-6) for p, g in tr.groupby('precursor_id')}
pop = pd.read_csv(f'{V}/raw01/population.tsv', sep='\t'); s = np.load(f'{V}/raw01/scores_D11FWT.npy').astype(np.float64)
w = pop.assign(s=s).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').copy()
for f_, g in w[w.is_decoy == 1].groupby('fold'):
    med = g.s.median(); mad = 1.4826 * (g.s - med).abs().median(); w.loc[w.fold == f_, 's'] = (w.loc[w.fold == f_, 's'] - med) / mad
w = w.set_index('precursor_id')
st = pd.read_csv(f'{V}/v4/dd_diann_raw01.tsv', sep='\t'); band = pd.read_csv(f'{V}/v4/dd_band_raw01.tsv', sep='\t')
D = st[st.grp == 'D'].precursor_id.values; lo, hi = w.loc[D, 's'].quantile([0.25, 0.75])
inb = lambda g: band[(band.grp == g) & band.precursor_id.map(w.s).between(lo, hi)].precursor_id.values
rng = np.random.default_rng(3)
G = {'D': D, 'DEC': inb('DEC'), 'ENT': inb('ENT'), 'E': rng.choice(st[st.grp == 'E'].precursor_id.values, 600, replace=False)}
idx = pd.read_parquet(f'{OUT}/scan_index.parquet'); m2 = idx[idx.ms == 2].sort_values('rt')
od = oms.OnDiscMSExperiment(); od.openFile(f'{OUT}/raw01.mzML', False)
rows = []
for grp, ids in G.items():
    for pid in ids:
        pmz, apex, f = pm[pid], w.loc[pid, 'rt_apex'], fr[pid]
        sc = m2[(m2.lo <= pmz) & (m2.hi >= pmz)]; sc = sc[(sc.rt > apex - 6) & (sc.rt < apex + 6)]
        if sc.empty: continue
        j0 = int((sc.rt - apex).abs().values.argmin()); rec = {'grp': grp, 'precursor_id': pid, 's': w.loc[pid, 's']}
        npts = 0
        for k, (_, r) in enumerate(sc.iterrows()):
            mz, it = od.getSpectrum(int(r.idx)).get_peaks()
            hit = [np.any(np.abs(mz - x) / x * 1e6 <= PPM) for x in f]
            share = sum(it[np.abs(mz - x) / x * 1e6 <= PPM].sum() for x in f) / max(it.sum(), 1)
            if abs(r.rt - apex) <= 3 and sum(hit) >= 3: npts += 1
            if k - j0 in (-1, 0, 1): rec[f'n_frag_{k - j0}'] = sum(hit); rec[f'tic_share_{k - j0}'] = share
        rec['pts_ge3'] = npts; rows.append(rec)
df = pd.DataFrame(rows); df['n_frag_best'] = df[['n_frag_-1', 'n_frag_0', 'n_frag_1']].max(1)
df.to_csv(f'{OUT}/apex_scan_evidence.tsv', sep='\t', index=False)
print(f'score band (D interquartile range) {lo:.1f}..{hi:.1f}; groups', df.groupby('grp').size().to_dict())
print('medians:'); print(df.groupby('grp')[['s', 'n_frag_0', 'n_frag_best', 'tic_share_0', 'pts_ge3']].median().round(3).to_string())
print('share with >= 4 of 6 library fragments in the apex scan (k = 0):'); print(df.groupby('grp').n_frag_0.apply(lambda v: round((v >= 4).mean(), 3)).to_string())
print('share with >= 2 scans (apex +- 3 s) holding >= 3 fragments:'); print(df.groupby('grp').pts_ge3.apply(lambda v: round((v >= 2).mean(), 3)).to_string())
def auc(a, b):
    r = pd.Series(np.r_[a, b]).rank().values; return (r[:len(a)].sum() - len(a) * (len(a) + 1) / 2) / (len(a) * len(b))
dd = df[df.s.between(lo, hi)]
for k in ('n_frag_0', 'n_frag_best', 'tic_share_0', 'pts_ge3'):
    a = {g: dd[dd.grp == g][k].dropna().values for g in ('D', 'DEC', 'ENT')}
    print(f'{k}: AUC D vs DEC {auc(a["D"], a["DEC"]):.3f}, ENT vs DEC {auc(a["ENT"], a["DEC"]):.3f} (scores within the D interquartile range)')
