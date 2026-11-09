"""k11 library build (PLAN.md). Label-free: uses library fragments only.
Writes library.tsv (targets unchanged + one chosen decoy per target), choice.tsv and build.txt."""
import sqlite3, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
HERE = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-commatched-decoys-20261109'
SRC = '/home/singjc/Documents/datasets/proteobench-astral-entrapment/idfree_calibration_research/library_exact_tsv_20261017_01/ProteoBenchFASTA_Entrapment_Human_with_contaminants_entrapment_pep_100k_randomly_sampled.decoys.tsv'
TPQP = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-design-raw01-20261107/rethink/decoyregen/targets.pqp'
COLS = ['PrecursorMz', 'ProductMz', 'LibraryIntensity', 'TransitionGroupId', 'Decoy']
log = open(HERE + '/build.txt', 'w')
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); log.write(s + '\n'); log.flush()
src = pd.read_csv(SRC, sep='\t', usecols=COLS)
tid2n = dict(pd.read_sql('select ID, TRAML_ID from PRECURSOR', sqlite3.connect('file:' + TPQP + '?mode=ro', uri=True)).values)
def top6(df):
    df = df.sort_values(['g', 'LibraryIntensity'], ascending=[True, False], kind='stable').groupby('g').head(6)
    return df
tg = src[src.Decoy == 0].rename(columns={'TransitionGroupId': 'g'}); tg = top6(tg)
tg['win'] = np.floor((tg.PrecursorMz - 380) / 2).astype(int)
ref = {w: np.sort(g.ProductMz.values) for w, g in tg.groupby('win')}
TF = tg.groupby('g').ProductMz.apply(np.array).to_dict(); TW = tg.groupby('g').win.first().to_dict()
def count(f, w):
    a = ref.get(w)
    if a is None: return 0
    return np.searchsorted(a, f * (1 + 17e-6)) - np.searchsorted(a, f * (1 - 17e-6))
def score(df, target_of):
    """commonness (excluding own target) and own-target overlap at 34 ppm, per candidate group."""
    out = {}
    for g, d in df.groupby('g'):
        t = target_of(g); f = d.ProductMz.values; w = int(np.floor((d.PrecursorMz.values[0] - 380) / 2))
        tf = TF.get(t, np.array([]))
        c = np.array([count(x, w) for x in f], dtype=float)
        if TW.get(t) == w and len(tf):
            c -= (np.abs(f[:, None] - tf[None, :]) <= 17e-6 * f[:, None]).sum(1)
        ov = int((np.abs(f[:, None] - tf[None, :]) <= 34e-6 * f[:, None]).any(1).sum()) if len(tf) else 0
        out[t] = (float(np.log1p(np.maximum(c, 0)).mean()), ov)
    return out
# target commonness (self excluded via the same rule)
ct = score(tg.assign(), lambda g: g)
say('targets', len(ct))
cand = {'PR': score(top6(src[src.Decoy == 1].rename(columns={'TransitionGroupId': 'g'})), lambda g: g[len('DECOY_'):])}
for s in (1, 2, 3, 4):
    d = pd.read_csv(f'{HERE}/shuf{s}.tsv', sep='\t', usecols=COLS).rename(columns={'TransitionGroupId': 'g'})
    cand[f'S{s}'] = score(top6(d), lambda g: tid2n[int(g[len('DECOY_'):])]); say('scored', f'S{s}', len(cand[f'S{s}']))
rows = []
for t, (c0, _) in ct.items():
    best = None
    for k in ('PR', 'S1', 'S2', 'S3', 'S4'):
        v = cand[k].get(t)
        if v is None or v[1] > 1: continue
        dist = abs(v[0] - c0)
        if best is None or dist < best[1] - 1e-12: best = (k, dist, v[0], v[1])
    if best is None and t in cand['PR']: best = ('PR', abs(cand['PR'][t][0] - c0), cand['PR'][t][0], cand['PR'][t][1])
    if best is not None: rows.append((t, best[0], c0, best[2], best[3]))
ch = pd.DataFrame(rows, columns=['target', 'choice', 'com_target', 'com_decoy', 'overlap']); ch.to_csv(HERE + '/choice.tsv', sep='\t', index=False)
pr = np.array([v[0] for v in cand['PR'].values()])
say('decoys chosen', len(ch), 'shares', ch.choice.value_counts(normalize=True).round(3).to_dict())
say(f'commonness mean: targets {np.mean([v[0] for v in ct.values()]):.3f} | old pseudo-reverse {pr.mean():.3f} | chosen {ch.com_decoy.mean():.3f}')
for nm, x in (('old', pr), ('chosen', ch.com_decoy.values)):
    tv = np.array([v[0] for v in ct.values()]); y = np.r_[np.ones(len(tv)), np.zeros(len(x))]
    say(f'AUC target vs {nm} decoy commonness: {roc_auc_score(y, np.r_[tv, x]):.3f}')
say('own-target overlap of chosen: share >=2', round((ch.overlap >= 2).mean(), 4))
# write the library TSV: targets unchanged, then the chosen decoys renamed DECOY_<target>
full = pd.read_csv(SRC, sep='\t', dtype=str, keep_default_na=False)
parts = [full[full.Decoy == '0']]
pick = ch.groupby('choice').target.apply(set).to_dict()
d = full[full.Decoy == '1']; d = d[d.TransitionGroupId.str[len('DECOY_'):].isin(pick.get('PR', set()))]; parts.append(d); del full
for s in (1, 2, 3, 4):
    want = pick.get(f'S{s}', set())
    d = pd.read_csv(f'{HERE}/shuf{s}.tsv', sep='\t', dtype=str, keep_default_na=False)
    t = d.TransitionGroupId.str[len('DECOY_'):].astype(int).map(tid2n); d = d[t.isin(want)].copy(); t = t[d.index]
    # Shuffles can repeat an annotation at the same m/z within a group; the original row id keeps transition IDs unique.
    d['TransitionGroupId'] = 'DECOY_' + t; d['TransitionId'] = d.TransitionGroupId + '_' + d.Annotation + f'_S{s}_' + d.TransitionId
    parts.append(d)
out = pd.concat(parts, ignore_index=True)
say('rows', len(out), 'decoy groups', out[out.Decoy == '1'].TransitionGroupId.nunique(), 'duplicate transition ids', int(out.TransitionId.duplicated().sum()))
out.to_csv(HERE + '/library.tsv', sep='\t', index=False)
say('done')
