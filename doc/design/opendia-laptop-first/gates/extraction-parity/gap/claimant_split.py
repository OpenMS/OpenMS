"""Precheck for a presence-weighted (joint) explaining-away fit. Diagnostic only, labels only in the grouping.
Recomputes features_I's per-link contributions (same rules) and splits each window's explained area by the claimant
window's own current score (C++ CrossFitScorer seed 0): claimants above the run's 1 % target cut (likely present) vs
below (likely absent). Compares DIA-NN-supported rejected IDs with the decoy winners matched in score.
Usage: claimant_split.py <raw01|raw02>"""
import json, sqlite3, sys, numpy as np, pandas as pd
RAW = sys.argv[1]
G = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/gap-20261114'; B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
D = f'{G}/armE_{RAW}'
C0, C1, K, BIN, CH = 12, 20, 12, 30.0 / 32, 2_000_000
TOL = float(json.load(open(f'{D}/cal/cal.json'))['ms2_mz']['tolerance_half_width_ppm']) * 1e-6
iso = pd.read_csv(f'{G}/iso_windows.tsv', sep='\t')[['lo', 'hi']].sort_values('lo').values
pop = pd.read_csv(f'{D}/population.tsv', sep='\t')
n = len(pop)
f = pd.read_csv(f'{D}/search.features.tsv', sep='\t', usecols=['candidate_row', 'rt_start', 'rt_end']).set_index('candidate_row')
pop['rt_start'] = f.rt_start.reindex(pop.sample_id).values; pop['rt_end'] = f.rt_end.reindex(pop.sample_id).values
c = sqlite3.connect(f'file:{D}/library12_pop.pqp?mode=ro', uri=True)
pmz = pd.read_sql('select ID pid, PRECURSOR_MZ mz from PRECURSOR', c).set_index('pid').mz
t = pd.read_sql("""select m.PRECURSOR_ID pid, t.PRODUCT_MZ mz from TRANSITION t join TRANSITION_PRECURSOR_MAPPING m
  on m.TRANSITION_ID=t.ID order by m.PRECURSOR_ID, t.ID""", c)
t['slot'] = t.groupby('pid').cumcount(); t = t[t.slot < K]
uniq, inv = np.unique(pop.precursor_id.values, return_inverse=True)
FMu = np.full((len(uniq), K), np.nan)
i = np.searchsorted(uniq, t.pid.values); ok = (i < len(uniq)) & (uniq[np.minimum(i, len(uniq) - 1)] == t.pid.values)
FMu[i[ok], t.slot.values[ok]] = t.mz.values[ok]
mz = pmz.reindex(pop.precursor_id).values
k = np.searchsorted(iso[:, 0], mz, side='right') - 1
pop['iso'] = np.where((k >= 0) & (mz < iso[np.clip(k, 0, None), 1]), k, -1)
X = np.load(f'{D}/xic12.npy', mmap_mode='r'); L = np.asarray(np.load(f'{D}/lib12.npy', mmap_mode='r'), np.float64)
A = np.zeros((n, K)); CS = np.zeros((n, 33))
for b in range(0, n, 100000):
    x = np.asarray(X[b:b + 100000], np.float64); A[b:b + 100000] = x[:, :, C0:C1].sum(2); CS[b:b + 100000, 1:] = np.cumsum(x.sum(1), 1)
pid = pop.precursor_id.values; apex = pop.rt_apex.values
ra, rb, ma = [], [], []
for _, g in pop[pop.iso >= 0].groupby('iso'):
    g = g.sort_values('rt_start'); idx = g.index.values; st, en = g.rt_start.values, g.rt_end.values
    last = np.searchsorted(st, en, side='right'); cnt = np.maximum(last - np.arange(len(g)) - 1, 0)
    I = np.repeat(np.arange(len(g)), cnt); J = I + 1 + (np.arange(len(I)) - np.repeat(np.cumsum(cnt) - cnt, cnt))
    I, J = idx[I], idx[J]; keep = pid[I] != pid[J]; I, J = I[keep], J[keep]
    for s0 in range(0, len(I), CH):
        a, b = I[s0:s0 + CH], J[s0:s0 + CH]
        M = np.abs(FMu[inv[a]][:, :, None] - FMu[inv[b]][:, None, :]) <= TOL * FMu[inv[a]][:, :, None]
        hit = M.any((1, 2))
        if not hit.any(): continue
        M = M[hit]; a, b = a[hit], b[hit]
        m_ab = np.where(M.any(2), M.argmax(2), -1).astype(np.int8); m_ba = np.where(M.any(1), M.argmax(1), -1).astype(np.int8)
        ra += [a, b]; rb += [b, a]; ma += [(m_ab, m_ba), (m_ba, m_ab)]
mab = np.concatenate([x[0] for x in ma]); mba = np.concatenate([x[1] for x in ma]); ra, rb = np.concatenate(ra), np.concatenate(rb)
d = np.abs(apex[ra] - apex[rb]); o = np.lexsort((d, pid[rb], ra)); ra, rb, mab, mba = ra[o], rb[o], mab[o], mba[o]
first = np.r_[True, (ra[1:] != ra[:-1]) | (pid[rb][1:] != pid[rb][:-1])][:len(ra)]
ra, rb, mab, mba = ra[first], rb[first], mab[first], mba[first]
# Claimant window score class.
s = np.load(f'{G}/cppscore/{RAW}_s0.npy').astype(float)
sys.argv = ['x']
code = open(B + '/pairlib-20261113/scripts/evaluate_v2.py').read().split('rng = np.random.default_rng(5)')[0].replace("B + '/pairlib-20261113/lib/library.pqp", f"'{G}/armE_raw01/library.pqp")
exec(code)
pop['s'] = s; pop['row'] = np.arange(n)
w, best = tdc(pop, s)
tw = w[w.is_decoy == 0].sort_values('s', ascending=False); e = tw.precursor_id.map(ent).values
fdp = np.cumsum(e) * (1 + 1 / R) / np.arange(1, len(tw) + 1); n1 = int(np.where(fdp <= 0.01)[0].max() + 1)
cut_raw = s[tw.row.values[:n1].astype(int)].min()  # raw-score cut (scores are tail-calibrated per fold already)
present = s >= cut_raw
# Is the claimant window its precursor's winner?
is_winner = np.zeros(n, bool); is_winner[w.row.values.astype(int)] = True
contrib = {k_: np.zeros((n, K)) for k_ in ('present', 'absent_winner', 'absent_nonwinner')}
for s0 in range(0, len(ra), CH):
    a, b, m_ab, m_ba = ra[s0:s0 + CH], rb[s0:s0 + CH], mab[s0:s0 + CH].astype(int), mba[s0:s0 + CH].astype(int)
    Ab, Lb = A[b], L[b]; use = (m_ba < 0) & (Lb > 0) & (Ab > 0)
    scale = np.where(use.any(1), np.nanmedian(np.where(use, Ab / np.maximum(Lb, 1e-12), np.nan), 1), 0.0)
    sh = np.rint((apex[a] - apex[b]) / BIN).astype(int); lo, hi = np.clip(C0 + sh, 0, 32), np.clip(C1 + sh, 0, 32)
    own = CS[b, C1] - CS[b, C0]; r_ = np.arange(len(b))
    frac = np.where(own > 0, (CS[b, hi] - CS[b, lo]) / np.maximum(own, 1e-12), 0.0)
    ev = np.where(m_ab >= 0, scale[:, None] * Lb[r_[:, None], np.maximum(m_ab, 0)] * frac[:, None], 0.0)
    cls = np.where(present[b], 0, np.where(is_winner[b], 1, 2))
    for ci, k_ in enumerate(('present', 'absent_winner', 'absent_nonwinner')):
        m = cls == ci
        np.add.at(contrib[k_], a[m], ev[m])
tot = A.sum(1)
total = sum(contrib.values())
expl_all = np.minimum(A, total).sum(1)
share = {k_: np.where(tot > 0, np.minimum(A, v).sum(1) / np.maximum(tot, 1e-12), 0) for k_, v in contrib.items()}
expl_present = np.where(tot > 0, np.minimum(A, contrib['present']).sum(1) / np.maximum(tot, 1e-12), 0)
# Groups.
c2 = sqlite3.connect(f'file:{G}/armE_raw01/library.pqp?mode=ro', uri=True)
pre = pd.read_sql("select p.ID precursor_id,p.CHARGE z,pe.MODIFIED_SEQUENCE ms,p.DECOY decoy from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID", c2)
pre['key'] = pre.ms + '/' + pre.z.astype(str); tg = pre[pre.decoy == 0].set_index('key')
REF = {'raw01': B + '/diann-pairlib-20261113/report_raw01.parquet', 'raw02': B + '/diann-decoys-20261114/own/raw02/report.parquet'}[RAW]
ref = pd.read_parquet(REF, columns=['Modified.Sequence', 'Precursor.Charge', 'RT', 'Q.Value'])
ref['key'] = ref['Modified.Sequence'] + '/' + ref['Precursor.Charge'].astype(str); ref = ref[ref['Q.Value'] <= 0.01].drop_duplicates('key'); ref = ref[ref.key.isin(tg.index)]
ref['pid'] = tg.loc[ref.key, 'precursor_id'].values; ref = ref[ref.pid.map(ent) == 0].set_index('pid')
acc = set(tw.precursor_id.iloc[:n1]); W = w.set_index('precursor_id')
right = [p for p in ref.index if p in W.index and abs(apex[int(W.row[p])] - ref.RT[p] * 60) <= 6]
rR = W.loc[[p for p in right if p not in acc], 'row'].astype(int).values; rA = W.loc[[p for p in right if p in acc], 'row'].astype(int).values
dw = w[w.is_decoy == 1].sort_values('s'); ds = dw.s.values
ix = np.clip(np.searchsorted(ds, W.loc[[p for p in right if p not in acc], 's'].values), 1, len(ds) - 2)
rD = np.unique(np.r_[dw.row.values[ix - 1], dw.row.values[ix], dw.row.values[ix + 1]].astype(int))
print(f'{RAW}: cut keeps {n1}; windows above the cut (likely present) {present.sum()} of {n}; groups: rejected {len(rR)}, matched decoys {len(rD)}, accepted {len(rA)}')
print('median (mean) explained share of the window central area, by claimant class:')
for nm, rows in (('rejected (DIA-NN right peak)', rR), ('matched decoy winners', rD), ('accepted', rA)):
    parts = ' | '.join(f'{k_} {np.median(v[rows]):.3f} ({v[rows].mean():.3f})' for k_, v in share.items())
    print(f'  {nm:30s} all {np.median(np.where(tot[rows] > 0, expl_all[rows] / np.maximum(tot[rows], 1e-12), 0)):.3f} | {parts}')
from sklearn.metrics import roc_auc_score
y = np.r_[np.ones(len(rR)), np.zeros(len(rD))]
for nm, v in (('explained, all claimants (current)', np.where(tot > 0, expl_all / np.maximum(tot, 1e-12), 0)), ('explained, present claimants only', expl_present)):
    print(f'  AUC rejected > matched decoys, -{nm}: {roc_auc_score(y, -np.r_[v[rR], v[rD]]):.3f}')
