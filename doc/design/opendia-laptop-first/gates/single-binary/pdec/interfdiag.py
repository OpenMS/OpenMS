"""PRECHECK_INTERF (diagnostic, labels for evaluation only). Input: runs/D_cur_s0/<run>.features.tsv.gz (cut columns) and the
predicted library. For accepted precursor winners (q<=0.01), classes:
  ES  entrapment target sharing >=3 fragments with a co-eluting identified HeLa target (|dRT|<=10 s, |dmz|<=12.5, 20 ppm)
  EO  other accepted entrapment targets
  DE  accepted decoys of entrapment pairs
  HT  accepted HeLa targets (sample)
Per interference feature: class medians, and AUC(ES > DE), AUC(EO > DE), AUC(HT > DE). The precursor winner is the row
with the highest score per precursor."""
import sqlite3, glob, bisect, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
LIB = 'hela_r5.predicted_decoys.pqp'
c = sqlite3.connect(f'file:{LIB}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID precursor_id, p.TRAML_ID t, p.DECOY d, p.PRECURSOR_MZ mz,
 min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
 from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID
 join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c).set_index('precursor_id')
pq['key'] = np.where(pq.d == 1, pq.t.str.replace('^DECOY_', '', regex=True), pq.t)
pq['pe'] = pq.key.map(pq[pq.d == 0].set_index('key').e)
fall = sorted(glob.glob('runs/[DH]_cur_s0/*.features.tsv.gz'))
frames = []
for f in fall:
    F = pd.read_csv(f, sep='\t'); F['run'] = f.split('/')[1][0]; frames.append(F)
F = pd.concat(frames)
feats = [x for x in F.columns if x.startswith('i_')]
W = F[F.score > -1e29].sort_values('score', ascending=False).drop_duplicates(['run', 'precursor_id'])
W = W[W.q_value <= 0.01].join(pq[['d', 'mz', 'e', 'pe']], on='precursor_id')
ids = W.precursor_id.tolist()
fr = pd.read_sql(f"select tm.PRECURSOR_ID p, tr.PRODUCT_MZ m from TRANSITION_PRECURSOR_MAPPING tm join TRANSITION tr on tr.ID=tm.TRANSITION_ID "
                 f"where tm.PRECURSOR_ID in ({','.join(map(str, ids))})", c)
frag = fr.groupby('p').m.apply(lambda s: np.sort(s.values)).to_dict()
def shared_run(Wr):
  hela = Wr[(Wr.d == 0) & (Wr.e == 0)].sort_values('rt_apex'); hr = hela.rt_apex.values; hm = hela.mz.values; hp = hela.precursor_id.values
  def shared(row):
    a = frag.get(row.precursor_id)
    if a is None: return 0
    lo = bisect.bisect_left(hr, row.rt_apex - 10); hi = bisect.bisect_right(hr, row.rt_apex + 10); best = 0
    for k in range(lo, hi):
        if hp[k] == row.precursor_id or abs(hm[k] - row.mz) > 12.5: continue
        b = frag.get(hp[k])
        if b is None: continue
        idx = np.searchsorted(b, a); n = 0
        for v, i in zip(a, idx):
            if any(0 <= j < len(b) and abs(b[j] - v) <= v * 20e-6 for j in (i - 1, i)): n += 1
        best = max(best, n)
    return best
  E = Wr[(Wr.d == 0) & (Wr.e == 1)].copy(); E['sh'] = E.apply(shared, axis=1); return E
E = pd.concat([shared_run(W[W.run == r]) for r in sorted(W.run.unique())])
print('runs', sorted(W.run.unique()))
cls = {'ES': E[E.sh >= 3], 'EO': E[E.sh < 3], 'DE': W[(W.d == 1) & (W.pe == 1)],
       'HT': W[(W.d == 0) & (W.e == 0)].sample(min(2000, int(((W.d == 0) & (W.e == 0)).sum())), random_state=0)}
print('accepted winners', len(W), {k: len(v) for k, v in cls.items()})
def auc(a, b):
    a = a.dropna(); b = b.dropna()
    return mannwhitneyu(a, b).statistic / (len(a) * len(b)) if len(a) and len(b) else float('nan')
for x in feats + ['score']:
    med = '  '.join(f'{k} {v[x].median():.3f}' for k, v in cls.items())
    print(f'{x:16s} median: {med}   AUC(ES>DE) {auc(cls["ES"][x], cls["DE"][x]):.3f}  AUC(EO>DE) {auc(cls["EO"][x], cls["DE"][x]):.3f}  AUC(HT>DE) {auc(cls["HT"][x], cls["DE"][x]):.3f}')
print('-- fraction with i_explained >= t, and with i_libcos_resid == 0, per class')
for t in (0.95, 0.99, 0.999):
    print(f'i_explained>={t}: ' + '  '.join(f'{k} {(v.i_explained>=t).mean():.3f}' for k, v in cls.items()))
print('i_libcos_resid==0: ' + '  '.join(f'{k} {(v.i_libcos_resid<=1e-6).mean():.3f}' for k, v in cls.items()))
# what a label-free rule (drop windows fully explained by co-eluting claimants) would remove among all accepted winners
allW = W.copy(); m = allW.i_explained >= 0.999
print(f'all accepted winners: {len(allW)}  fully explained {m.sum()} ({m.mean():.3%}); of those decoys {int((allW.d[m]==1).sum())}, entrapments {int(((allW.d[m]==0)&(allW.e[m]==1)).sum())}, HeLa {int(((allW.d[m]==0)&(allW.e[m]==0)).sum())}')
