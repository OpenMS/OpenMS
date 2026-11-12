"""raw01, q <= 1%: DIA-NN peaks inside our RT window that none of our <= 5 reported peaks contains (cut by the
co-elution ranking?), and DIA-NN peaks we contain but do not rank first. Our co-elution is recomputed on our traces
over DIA-NN's [RT.Start, RT.Stop] (at least +-2 scans around the apex), as pickXicPeaks defines it:
sum over fragments of max(0, Pearson r with the sum of the other fragments). Computed with all 6 of our fragments
and with only the fragments DIA-NN kept (it drops b1/b2/y1/y2 etc.). Evaluation only.
"""
import sqlite3, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; X = B + '/extraction-parity-20261112'
out = open(X + '/rtmono/PEAK_CHOICE_raw01.txt', 'w')
def say(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')
rows = pd.read_csv(X + '/rtmono/peak_parity_rows_raw01.tsv', sep='\t')
q1 = rows[rows.grp == 'q<=1%'].copy()
idm = pd.read_csv(X + '/idmap.tsv', sep='\t')
q1 = q1.merge(idm[['precursor_id', 'Precursor.Id']], on='precursor_id')
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
tr = pd.read_sql("""select tp.PRECURSOR_ID precursor_id, t.PRODUCT_MZ product_mz, t.ANNOTATION annotation
  from TRANSITION t join TRANSITION_PRECURSOR_MAPPING tp on tp.TRANSITION_ID=t.ID""", c)
tr = tr[tr.precursor_id.isin(set(q1.precursor_id))].sort_values(['precursor_id', 'product_mz'])
tr['slot'] = tr.groupby('precursor_id').cumcount()
dn = pd.read_parquet(X + '/diann/report_raw01_xic/raw01.xic.parquet', columns=['pr', 'feature'])
kept = dn[dn.pr.isin(set(q1['Precursor.Id'])) & ~dn.feature.isin(['index', 'ms1'])].drop_duplicates()
kept = set(zip(kept.pr, kept.feature))
us = pd.read_csv(X + '/rtmono/traces.tsv', sep='\t', usecols=['precursor_id', 'rt', 'slot', 'intensity'])
us = us[us.precursor_id.isin(set(q1.precursor_id))].merge(tr[['precursor_id', 'slot', 'annotation']], on=['precursor_id', 'slot'])
feat = pd.read_csv(X + '/rtmono/search.features.tsv', sep='\t', usecols=['precursor_id', 'is_decoy', 'rt_apex', 'rt_start', 'rt_end', 'coelution', 'peak_rank', 'apex_intensity'])
feat = feat[(feat.is_decoy == 0) & feat.precursor_id.isin(set(q1.precursor_id))]
fg = dict(tuple(feat.groupby('precursor_id')))
ug = dict(tuple(us.groupby('precursor_id')))

def coel(M, a, b, ia):
    a = min(a, ia - 2); b = max(b, ia + 2); a = max(a, 0); b = min(b, M.shape[0] - 1)
    W = M[a:b + 1]; tot = W.sum(1); s = 0.0
    for k in range(W.shape[1]):
        x = W[:, k]; y = tot - x
        if x.std() > 0 and y.std() > 0:
            s += max(0.0, float(np.corrcoef(x, y)[0, 1]))
    return s
res = []
for r in q1.itertuples():
    g = ug[r.precursor_id]
    P = g.pivot_table(index='rt', columns='annotation', values='intensity', aggfunc='sum').fillna(0)
    t = P.index.values
    if r.RT < t[0] or r.RT > t[-1]:
        continue
    ia = int(np.argmin(np.abs(t - r.RT))); a = int(np.argmin(np.abs(t - r._5))); b = int(np.argmin(np.abs(t - r._6)))
    M = P.values; cols = list(P.columns)
    sub = [i for i, cname in enumerate(cols) if (r._21, cname) in kept]
    f = fg[r.precursor_id]
    top = f[f.peak_rank == 0].iloc[0]
    res.append(dict(precursor_id=r.precursor_id, e=r.e, contains=r.contains, is_top=r.is_top, n_kept=len(sub),
                    coel_dn_all=coel(M, a, b, ia), coel_dn_kept=coel(M[:, sub], a, b, ia) if len(sub) >= 2 else np.nan,
                    coel_min_reported=f.coelution.min(), coel_top=top.coelution,
                    top_dist_scans=abs(top.rt_apex - r.RT) / 1.646,
                    dn_apex_height=M[ia].sum(), top_apex_height=top.apex_intensity))
R = pd.DataFrame(res)
nc = R[~R.contains]
say(f'q<=1% with DIA-NN apex inside our window: {len(R)}; not contained in any of our 5 peaks: {len(nc)} (entrapments {int(nc.e.sum())})')
say(f'  DIA-NN peak co-elution (our 6 fragments) below the 5th reported peak (cut by ranking): {(nc.coel_dn_all < nc.coel_min_reported).sum()}')
say('  co-elution DIA-NN peak, all ours / kept only:', nc.coel_dn_all.round(2).tolist(), '/', nc.coel_dn_kept.round(2).tolist())
say('  co-elution of our rank-0 peak:', nc.coel_top.round(2).tolist())
say('  our rank-0 apex height / DIA-NN apex height (our consensus):', (nc.top_apex_height / nc.dn_apex_height).round(2).tolist())
nt = R[R.contains & ~R.is_top]
say(f'\ncontained but not our rank-0: {len(nt)} (entrapments {int(nt.e.sum())})')
for k, s in [('all', 'coel_dn_all'), ('kept', 'coel_dn_kept')]:
    # normalise by fragment count so 6-fragment and 4-fragment scores compare
    n = 6 if k == 'all' else nt.n_kept
    say(f'  co-elution per fragment of DIA-NN peak ({k}) quantiles', (nt[s] / n).quantile([.1, .5, .9]).round(3).to_dict())
say('  rank-0 co-elution per fragment', (nt.coel_top / 6).quantile([.1, .5, .9]).round(3).to_dict())
say('  rank-0 apex distance from DIA-NN apex, scans', nt.top_dist_scans.quantile([.1, .5, .9]).round(1).to_dict())
# would ranking on DIA-NN's fragment subset pick DIA-NN's peak? (needs per-peak recompute on kept fragments)
better = 0; tot = 0
for r in nt.itertuples():
    g = ug[r.precursor_id]; P = g.pivot_table(index='rt', columns='annotation', values='intensity', aggfunc='sum').fillna(0)
    t = P.index.values; M = P.values
    pr = q1.set_index('precursor_id').loc[r.precursor_id]
    sub = [i for i, cname in enumerate(P.columns) if (pr['Precursor.Id'], cname) in kept]
    if len(sub) < 2:
        continue
    f = fg[r.precursor_id]; scores = []
    for p in f.itertuples():
        ia = int(np.argmin(np.abs(t - p.rt_apex))); a = int(np.argmin(np.abs(t - p.rt_start))); b = int(np.argmin(np.abs(t - p.rt_end)))
        scores.append((coel(M[:, sub], a, b, ia), p.rt_start <= pr.RT <= p.rt_end))
    tot += 1; better += max(scores)[1]
say(f'  re-ranking our 5 peaks by co-elution on DIA-NN\'s kept fragments only picks DIA-NN\'s peak: {better} of {tot}')
R.to_csv(X + '/rtmono/peak_choice_rows_raw01.tsv', sep='\t', index=False)
