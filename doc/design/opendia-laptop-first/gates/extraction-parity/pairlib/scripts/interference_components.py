"""Interference components (step 1 of DESIGN_INTERFERENCE_GRAPH.md): how often do candidate peaks of different library
precursors claim the same fragment signal?

Graph (label-free): node = search candidate window (every row of search.features.tsv, targets and decoys).
Edge between windows a, b of different precursors when
  - their precursors are isolated in the same DIA window (m/z inside the same [lo, hi) of the run's MS2 scans),
  - their peak ranges [rt_start, rt_end] overlap,
  - at least one library fragment of a matches one of b within the run's fitted fragment tolerance (cal.json
    ms2_mz.tolerance_half_width_ppm). The number of matched fragment pairs is the edge weight.
Reported:
  - connected components over all windows (edges with weight >= 1 and >= 2);
  - for the arm's precursor winners at q <= 0.01 / 0.05 (paired TDC as evaluate_v2.py): share with a claimant
    (an edge of weight >= 1, >= 2, >= 3), the claimant's class (target/decoy), and whether the claimant scores higher.
  - --labels (descriptive, evaluation only): the same by human / entrapment / decoy and, for entrapments, whether a
    claimant has the same composition and charge (its FDRBench parent).
Usage: interference_components.py <root> <raw> <arm> <library.pqp> <search.features.tsv> <cal.json> <scan_index or mzML> [--labels]"""
import json, sqlite3, sys, numpy as np, pandas as pd

ROOT, RAW, ARM, LIB, FEAT, CAL, SCANS = sys.argv[1:8]
LABELS = '--labels' in sys.argv
cal = json.load(open(CAL))
TOL = float(cal['ms2_mz']['tolerance_half_width_ppm']) * 1e-6

# DIA isolation windows of the run.
if SCANS.endswith('.parquet'):
    sc = pd.read_parquet(SCANS)
    iso = sc[sc.ms == 2][['lo', 'hi']].drop_duplicates().sort_values('lo').values
else:
    import pyopenms as oms
    e = oms.MSExperiment(); oms.MzMLFile().load(SCANS, e)
    iso = sorted({(s.getPrecursors()[0].getMZ() - s.getPrecursors()[0].getIsolationWindowLowerOffset(),
                   s.getPrecursors()[0].getMZ() + s.getPrecursors()[0].getIsolationWindowUpperOffset())
                  for s in e if s.getMSLevel() == 2})
    iso = np.array(iso)
print(f'{RAW}: {len(iso)} isolation windows, fragment tolerance {TOL * 1e6:.2f} ppm')

c = sqlite3.connect(f'file:{LIB}?mode=ro', uri=True)
pre = pd.read_sql("""select p.ID precursor_id, p.PRECURSOR_MZ mz, p.CHARGE z, p.DECOY decoy, pe.UNMODIFIED_SEQUENCE seq from PRECURSOR p
  join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID""", c).drop_duplicates('precursor_id').set_index('precursor_id')
tr = pd.read_sql("""select m.PRECURSOR_ID precursor_id, t.PRODUCT_MZ fmz from TRANSITION t
  join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID""", c)
fr = tr.assign(r=tr.groupby('precursor_id').cumcount())
KF = int(fr.r.max()) + 1
FM = np.full((int(pre.index.max()) + 1, KF), np.nan)
FM[fr.precursor_id.values, fr.r.values] = fr.fmz.values  # fragment m/z matrix by precursor id (NaN = none)

# Isolation window of each precursor (the window whose [lo, hi) holds its m/z; -1 outside).
k = np.searchsorted(iso[:, 0], pre.mz.values, side='right') - 1
ok = (k >= 0) & (pre.mz.values < iso[np.clip(k, 0, None), 1])
pre['iso'] = np.where(ok, k, -1)

cand = pd.read_csv(FEAT, sep='\t', usecols=['candidate_row', 'precursor_id', 'is_decoy', 'rt_start', 'rt_apex', 'rt_end'])
cand['iso'] = cand.precursor_id.map(pre.iso)
cand = cand[cand.iso >= 0].reset_index(drop=True)
print(f'{len(cand)} candidate windows of {cand.precursor_id.nunique()} precursors')


def shared(pa, pb):
    """Per pair: number of fragments of a with a partner in b within the tolerance."""
    fa, fb = FM[pa][:, :, None], FM[pb][:, None, :]
    return (np.abs(fa - fb) <= TOL * fa).any(axis=2).sum(axis=1)


# Edges: within each isolation window, sorted by rt_start, window i pairs with the following windows j whose start is
# at or before i's end (so the ranges overlap). Vectorised over pairs, in chunks.
ea, eb, ew = [], [], []
for iso_k, g in cand.groupby('iso'):
    g = g.sort_values('rt_start')
    idx, st, en, pid = g.index.values, g.rt_start.values, g.rt_end.values, g.precursor_id.values
    last = np.searchsorted(st, en, side='right')  # first j with st[j] > en[i]
    cnt = np.maximum(last - np.arange(len(g)) - 1, 0)
    I = np.repeat(np.arange(len(g)), cnt)
    J = I + 1 + (np.arange(len(I)) - np.repeat(np.cumsum(cnt) - cnt, cnt))
    keep = pid[I] != pid[J]
    I, J = I[keep], J[keep]
    for s0 in range(0, len(I), 2_000_000):
        a, b = I[s0:s0 + 2_000_000], J[s0:s0 + 2_000_000]
        wv = shared(pid[a], pid[b])
        m = wv > 0
        ea.append(idx[a[m]]); eb.append(idx[b[m]]); ew.append(wv[m])
ea, eb, ew = np.concatenate(ea), np.concatenate(eb), np.concatenate(ew)
print(f'edges (>= 1 shared fragment) {len(ew)}; weight counts', dict(zip(*np.unique(ew, return_counts=True))))
EDGES_OUT = next((a.split('=', 1)[1] for a in sys.argv if a.startswith('--edges-out=')), None)
if EDGES_OUT:  # candidate_row pairs and weights, for later diagnostics
    np.savez_compressed(EDGES_OUT, a=cand.candidate_row.values[ea], b=cand.candidate_row.values[eb], w=ew.astype(np.int8),
                        pa=cand.precursor_id.values[ea], pb=cand.precursor_id.values[eb])


def components(n, a, b):
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    return connected_components(coo_matrix((np.ones(len(a)), (a, b)), shape=(n, n)), directed=False)[1]


for wmin in (1, 2):
    m = ew >= wmin
    comp = components(len(cand), ea[m], eb[m])
    size = pd.Series(comp).value_counts()
    nprec = cand.assign(c=comp).groupby('c').precursor_id.nunique()
    print(f'components (weight >= {wmin}): windows in non-singletons {np.mean(size.reindex(comp).values > 1):.3f}; '
          f'precursors per component p50/p90/p99/max {np.percentile(nprec, [50, 90, 99]).round(1).tolist()} {nprec.max()}')

# Winners of the arm (paired TDC as evaluate_v2.py) on the population, mapped to candidate rows.
pop = pd.read_csv(f'{ROOT}/{RAW}/population.tsv', sep='\t')
s = np.load(f'{ROOT}/{RAW}/scores_{ARM}.npy')
w = pop.assign(s=s).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').copy()
for f_, gq in w[w.is_decoy == 1].groupby('fold'):
    med = gq.s.median(); mad = 1.4826 * (gq.s - med).abs().median()
    w.loc[w.fold == f_, 's'] = (w.loc[w.fold == f_, 's'] - med) / mad
best = w.sort_values(['s', 'is_decoy'], ascending=[False, False]).drop_duplicates('pair_key').sort_values('s', ascending=False)
T = np.cumsum(1 - best.is_decoy.values); D = np.cumsum(best.is_decoy.values)
best['q'] = np.minimum.accumulate(np.minimum(1, (D + 1) / np.maximum(T, 1))[::-1])[::-1]
score_of = dict(zip(w.sample_id, w.s))  # standardised score of every precursor's winning window
row_of = pd.Series(cand.index.values, index=cand.candidate_row.values)

# Per candidate row: its claimants (other precursors' windows sharing signal), with weight.
nb = pd.DataFrame({'a': np.r_[ea, eb], 'b': np.r_[eb, ea], 'w': np.r_[ew, ew]})
nb['b_pid'] = cand.precursor_id.values[nb.b.values]; nb['b_dec'] = cand.is_decoy.values[nb.b.values]
nb['b_row'] = cand.candidate_row.values[nb.b.values]
# A claimant counts as a competing winner when that window is its precursor's winning window in the arm.
nb['b_win'] = nb.b_row.isin(set(w.sample_id))
nb['b_s'] = nb.b_row.map(score_of)

ent = None
if LABELS:
    ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
      from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID
      group by pp.PRECURSOR_ID""", c).set_index('precursor_id').e
    comp_key = (pre.seq.map(lambda x: ''.join(sorted(x))) + '/' + pre.z.astype(str))

for qcut in (0.01, 0.05):
    acc = best[best.q <= qcut]
    rows = acc.sample_id.map(row_of).dropna().astype(int)
    sub = nb[nb.a.isin(set(rows))]
    print(f'\n{RAW} {ARM} q <= {qcut}: {len(acc)} winners ({int((acc.is_decoy == 0).sum())} target, {int(acc.is_decoy.sum())} decoy)')
    cls = acc.assign(row=acc.sample_id.map(row_of))
    if LABELS:
        cls['cls'] = np.where(cls.is_decoy == 1, 'decoy', np.where(cls.precursor_id.map(ent) == 1, 'entrap', 'human'))
    else:
        cls['cls'] = np.where(cls.is_decoy == 1, 'decoy', 'target')
    for name, g in cls.groupby('cls'):
        r = set(g.row.dropna().astype(int))
        sg = sub[sub.a.isin(r)]
        line = f'  {name:7s} n {len(g):5d}'
        for wmin in (1, 2, 3):
            line += f' | claim>={wmin} {sg[sg.w >= wmin].a.nunique() / len(g):.3f}'
        win2 = sg[(sg.w >= 2) & sg.b_win]
        line += f' | claimant-winner>=2 {win2.a.nunique() / len(g):.3f}'
        stronger = win2[win2.b_s > win2.a.map(dict(zip(cls.row, cls.s)))]
        line += f' (stronger {stronger.a.nunique() / len(g):.3f})'
        if LABELS and name == 'entrap':
            ka = cand.precursor_id.values[sg.a.values]
            same = comp_key.reindex(ka).values == comp_key.reindex(sg.b_pid.values).values
            line += f' | claimant with same composition {sg[same].a.nunique() / len(g):.3f}'
        print(line)
