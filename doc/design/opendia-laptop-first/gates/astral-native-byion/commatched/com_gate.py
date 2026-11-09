"""k11 gate: commonness-matched decoys, emulated by reweighting (declared before running).
A commonness-matched decoy library makes the decoy commonness distribution equal to the targets'. Re-weighting the existing
decoys by w(bin) = P_target(bin)/P_decoy(bin) (10 quantile bins, library only, no labels, no run data) gives that same
marginal. If commonness drives the excess, high-commonness decoys already score like false targets and the weighted q is
calibrated. This bounds what a regenerated library can deliver through commonness alone.
Commonness as com_tdc.py: mean over top-6 fragments of log1p(#top-6 fragments of OTHER target precursors in the same 2 Th
window within +-17 ppm), on the native library (astral-entrapment-20261107/lib/library.pqp).
Native search+boundary+byion score (tool stage2_score, pooled decoy-standardised), plain pooled paired TDC.
Rule: GO for k11 (library regeneration) only if weighted q<=0.05 FDP_c <= 0.06 on BOTH raws, or drops by >= 0.01 on both.
Otherwise k11 is NO-GO: commonness-matched decoys are not expected to fix the 5% excess."""
import os, sqlite3, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; R = 148983 / 151239
HERE = os.path.dirname(os.path.abspath(__file__)); CACHE = HERE + '/commonness_native.parquet'
ARMS = {'raw01': B + '/astral-design-raw01-20261107/stage2-byion-design/byion', 'raw02': B + '/astral-raw02-byion-design-20261109/byion'}
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
p = pd.read_sql('select ID precursor_id, TRAML_ID, DECOY, PRECURSOR_MZ mz from PRECURSOR', c)
p['pair_key'] = np.where(p.DECOY == 1, p.TRAML_ID.str.replace('^DECOY_', '', regex=True), p.TRAML_ID)
if not os.path.exists(CACHE):
    tr = pd.read_sql('select m.PRECURSOR_ID pid, t.PRODUCT_MZ f, t.LIBRARY_INTENSITY li, t.ID t from TRANSITION t join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID', c)
    tr = tr.sort_values(['pid', 'li', 't'], ascending=[True, False, True]).groupby('pid').head(6)
    pr = p.rename(columns={'precursor_id': 'pid'}); pr['win'] = np.floor((pr.mz - 380) / 2).astype(int)
    tr = tr.merge(pr[['pid', 'win', 'DECOY']], on='pid')
    tg = tr[tr.DECOY == 0]
    ref = {w: (np.sort(g.f.values), g.pid.values[np.argsort(g.f.values)]) for w, g in tg.groupby('win')}
    out = []
    for (pid, w), g in tr.groupby(['pid', 'win']):
        a, ids = ref.get(w, (np.array([]), np.array([])))
        cnt = [int((ids[np.searchsorted(a, f * (1 - 17e-6)):np.searchsorted(a, f * (1 + 17e-6))] != pid).sum()) for f in g.f.values]
        out.append((pid, float(np.log1p(np.array(cnt)).mean()), int(g.DECOY.values[0])))
    pd.DataFrame(out, columns=['pid', 'com', 'd']).to_parquet(CACHE)
cm = pd.read_parquet(CACHE); COM = dict(zip(cm.pid, cm.com))
ct, cd = cm.com[cm.d == 0].values, cm.com[cm.d == 1].values
edges = np.unique(np.quantile(np.r_[ct, cd], np.linspace(0, 1, 11)))
binof = lambda x: np.clip(np.searchsorted(edges, x, 'right') - 1, 0, len(edges) - 2)
WB = (np.bincount(binof(ct), minlength=len(edges) - 1) / len(ct)) / np.maximum(np.bincount(binof(cd), minlength=len(edges) - 1) / len(cd), 1e-12)
print(f'library commonness mean target {ct.mean():.3f} decoy {cd.mean():.3f} | bin weights {np.round(WB, 2).tolist()}', flush=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c)
for run, arm in ARMS.items():
    f = pd.read_csv(arm + '/features.tsv', sep='\t', usecols=['sample_id', 'precursor_id', 'is_decoy', 'fold', 'stage2_score'])
    f = f.merge(p[['precursor_id', 'pair_key']], on='precursor_id').merge(ent, on='precursor_id', how='left'); f['e'] = f.e.fillna(0).astype(int)
    s = f.stage2_score.values; fold = f.fold.values
    w0 = f.sort_values(['stage2_score', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id'); z = s.copy()
    for k in np.unique(fold):
        g = w0[(w0.is_decoy == 1) & (w0.fold == k)].stage2_score; med = g.median(); z[fold == k] = (s[fold == k] - med) / (1.4826 * (g - med).abs().median())
    f['z'] = z
    w = f.sort_values(['z', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id')
    d = w.sort_values(['z', 'is_decoy'], ascending=[False, False]).drop_duplicates('pair_key').sort_values('z', ascending=False)
    dd = d.is_decoy.values == 1; T = np.cumsum(~dd)
    wd = np.where(dd, WB[binof(d.precursor_id.map(COM).fillna(0).values)], 0.0)
    out = []
    for nm, Dn in (('plain', np.cumsum(dd)), ('commonness-weighted', np.cumsum(wd))):
        q = np.minimum.accumulate(np.minimum(1, (Dn + 1) / np.maximum(T, 1))[::-1])[::-1]
        for qq in (0.01, 0.05):
            a = d[(~dd) & (q <= qq)]; out.append(f'{nm} q<={qq}: {len(a)} FDP_c {a.e.sum() * (1 + 1 / R) / max(1, len(a)):.3f}')
    print(f'{run}: ' + ' | '.join(out), flush=True)
