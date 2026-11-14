"""Calibration decoy-sensitivity probe (label-free; library and calibration diagnostics only).
Recomputes the stage-3 paired q-values from a diag.tsv with the panel restricted to pairs whose decoy shares fewer than
K of its library fragments (within 20 ppm) with its own target, for K = inf (as run), 3, 2, 1. Reports the anchors at
q <= 0.01. Usage: overlap_sim.py <library.pqp> <diag.tsv>"""
import sys, sqlite3, numpy as np, pandas as pd
lib, diag = sys.argv[1:3]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
pre = pd.read_sql('select ID precursor_id, TRAML_ID traml, DECOY decoy from PRECURSOR', c)
tid = pre[pre.decoy == 0].set_index('traml').precursor_id
pre['partner'] = np.where(pre.decoy == 1, pre.traml.str.replace('^DECOY_', '', regex=True).map(tid), -1)
d = pd.read_csv(diag, sep='\t', comment='#').merge(pre[['precursor_id', 'partner']], on='precursor_id')
ids = set(d.precursor_id)
tr = pd.read_sql('select m.PRECURSOR_ID id, t.PRODUCT_MZ mz, t.LIBRARY_INTENSITY li from TRANSITION t join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID', c)
tr = tr[tr.id.isin(ids)]
frag = {i: g.sort_values('li', ascending=False).mz.values[:6] for i, g in tr.groupby('id')}
dec = d[d.decoy == 1]
shared = {}
for r in dec.itertuples():
    f1, f2 = frag.get(r.precursor_id), frag.get(int(r.partner))
    shared[int(r.partner)] = 0 if f1 is None or f2 is None else int(sum(np.min(np.abs(f2 - x)) / x * 1e6 < 20 for x in f1))
d['pair'] = np.where(d.decoy == 1, d.partner, d.precursor_id).astype(int)
d['shared'] = d.pair.map(shared).fillna(0).astype(int)
print('pairs by shared fragments:', d[d.decoy == 0].shared.value_counts().sort_index().to_dict())


def anchors(sub):
    s = sub.set_index('precursor_id')
    comp = []
    for r in sub[sub.stage3_score.notna()].itertuples():
        p = int(r.partner) if r.decoy == 1 else None
        if r.decoy == 1:
            if p in s.index and s.loc[p, 'stage3_score'] == s.loc[p, 'stage3_score']:
                continue
            comp.append((r.stage3_score, 1)); continue
        dp = sub[(sub.decoy == 1) & (sub.pair == r.precursor_id)]
        ds = dp.stage3_score.iloc[0] if len(dp) else np.nan
        comp.append((r.stage3_score, 0) if not (ds == ds) or r.stage3_score > ds else (ds, 1))
    comp.sort(key=lambda x: -x[0])
    t = dcount = 0; fdr = []
    for sc, isd in comp:
        dcount += isd; t += 1 - isd; fdr.append((dcount + 1) / max(t, 1))
    q = np.minimum.accumulate(np.array(fdr)[::-1])[::-1]
    return int(sum(1 for (sc, isd), qq in zip(comp, q) if not isd and qq <= 0.01)), int(sum(isd for sc, isd in comp[:100]))


for k in (99, 3, 2, 1):
    sub = d[d.shared < k]
    a, top = anchors(sub)
    print(f'exclude pairs sharing >= {k} fragments: pairs {int((sub.decoy == 0).sum())}, anchors at q<=0.01 {a}, decoys in top 100 {top}')
