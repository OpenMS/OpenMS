"""PRECHECK_CANDREC part A (diagnostic, 21 Nov 2026): why do DIA-NN HeLa IDs get no OpenDIA window within 6 s of DIA-NN's RT?
Input: the PRECHECK_CEIL0 -features_out runs (current default, seed 0), the run's calibration.json and log, the DIA-NN report
of the same library level. DIA-NN labels are used for evaluation only.

The search runs on the refined RT prediction (second pass, PLAN_RTR). Each window carries abs_library_rt_error = |rt_apex - pred|,
so pred is recovered per precursor as the value rt_apex +- err that its windows agree on (one window: two possibilities).
Classes for a C-class ID (no window within 6 s of DIA-NN RT):
  rt_out    DIA-NN RT outside pred +- refined tolerance (both possibilities when ambiguous)
  near      a window apex 6-15 s from DIA-NN RT (apex disagreement, peak is probably found)
  cap       in window, the precursor has max_peaks (5) windows, none near DIA-NN RT
  picker    in window, fewer than 5 windows, none near (valley merge or no picked maximum)
Also: C share by quartile of DIA-NN quantity and of the candidate count, plus the nearest-window distance distribution.
Usage: precheck_candrec_a.py <run dir> <lib.pqp> <diann report.parquet> <tag>"""
import sys, glob, re, json, sqlite3, numpy as np, pandas as pd
RD, L, DN, TAG = sys.argv[1:5]
feat = glob.glob(f'{RD}/*_feat.parquet')[0]
log = open(RD.rstrip('/') + '.log').read()
tol = float(re.findall(r'rt_tolerance_s [0-9.]+ -> ([0-9.]+)', log)[-1])
F = pd.read_parquet(feat, columns=['precursor_id', 'is_decoy', 'rt_apex', 'abs_library_rt_error', 'score', 'q_value'])
F = F[F.is_decoy == 0]
c = sqlite3.connect(f'file:{L}?mode=ro', uri=True)
tg = pd.read_sql("""select p.ID pid, p.CHARGE z, pe.MODIFIED_SEQUENCE ms,
  min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID where p.DECOY=0 group by p.ID""", c)
tg = tg.assign(sk=tg.ms + '/' + tg.z.astype(str)).drop_duplicates('sk').set_index('sk')
d = pd.read_parquet(DN, columns=['Modified.Sequence', 'Precursor.Charge', 'Q.Value', 'RT', 'Precursor.Quantity', 'RT.Start', 'RT.Stop'])
d = d[d['Q.Value'] <= 0.01]
k = d['Modified.Sequence'] + '/' + d['Precursor.Charge'].astype(str)
d['pid'] = tg.pid.reindex(k).values; d['e'] = tg.e.reindex(k).values
d = d.dropna(subset=['pid']); d = d[d.e == 0]; d['pid'] = d.pid.astype(np.int64)
d = d.sort_values('Q.Value').drop_duplicates('pid').set_index('pid')
d['rt_s'] = d.RT * 60; d['dn_width'] = (d['RT.Stop'] - d['RT.Start']) * 60
tool = set(F.precursor_id[F.q_value <= 0.01])  # window q equals the precursor q of its best window in this output
T = F[F.precursor_id.isin(d.index)].copy()
T['drt'] = np.abs(T.rt_apex.values - d.rt_s.reindex(T.precursor_id).values)
# recover the predicted RT
T['p1'] = (T.rt_apex + T.abs_library_rt_error).round(1); T['p2'] = (T.rt_apex - T.abs_library_rt_error).round(1)
def pred_of(g):
    if len(g) == 1: return pd.Series({'pa': g.p1.iloc[0], 'pb': g.p2.iloc[0]})
    v = pd.concat([g.p1, g.p2]).round(0).value_counts()
    best = v.index[0]; cand = pd.concat([g.p1, g.p2]); m = cand[(cand.round(0) - best).abs() <= 1].median()
    return pd.Series({'pa': m, 'pb': m})
g = T.groupby('precursor_id')
P = g[['p1', 'p2']].apply(pred_of)
d['n'] = g.size().reindex(d.index).fillna(0).astype(int)
d['near'] = g.drt.min().reindex(d.index)
d['pa'] = P.pa.reindex(d.index); d['pb'] = P.pb.reindex(d.index)
d['err_a'] = d.rt_s - d.pa; d['err_b'] = d.rt_s - d.pb
d['err'] = np.where(d.err_a.abs() <= d.err_b.abs(), d.err_a, d.err_b)
d['ambig'] = (d.pa - d.pb).abs() > 1
d['missed'] = ~d.index.isin(list(tool))
d['C'] = ~(d.near <= 6)
d['qq'] = pd.qcut(np.log10(d['Precursor.Quantity'].clip(lower=1)), 4, labels=['Q1', 'Q2', 'Q3', 'Q4'])
M = d[d.missed & d.C].copy()
rt_out = M.err.abs() > tol
near = ~rt_out & (M.near <= 15)
cap = ~rt_out & ~near & (M.n >= 5)
M['cat'] = np.select([M.n == 0, rt_out, near, cap], ['no window at all', 'rt_out', 'near (6-15 s)', 'cap (5 windows)'], 'picker (<5 windows)')
print(f'{TAG}: refined RT tolerance +-{tol:.1f} s; DIA-NN HeLa IDs {len(d)}, missed by tool {d.missed.sum()}, of which C {len(M)}')
print(f'  all DIA-NN IDs: |RT - refined pred| p50 {d.err.abs().median():.1f} p95 {d.err.abs().quantile(.95):.1f} p99 {d.err.abs().quantile(.99):.1f} s; '
      f'outside the window {(d.err.abs() > tol).mean():.3f}; single-window (ambiguous pred) {d.ambig.mean():.3f}')
tab = pd.crosstab(M.cat, M.qq); tab['all'] = M.cat.value_counts(); tab['share'] = (tab['all'] / len(M)).round(3)
print(tab.to_string())
print(f'  windows per precursor: C-missed p50 {M.n.median():.0f} (share with 5: {(M.n >= 5).mean():.3f}); '
      f'found-and-accepted p50 {d[~d.missed].n.median():.0f} ({(d[~d.missed].n >= 5).mean():.3f})')
print('  nearest-window distance for C-missed (s): ' + ' '.join(f'p{q} {M.near.quantile(q / 100):.1f}' for q in (10, 25, 50, 75, 90)))
print(f'  DIA-NN peak width (RT.Stop-RT.Start) p50: C-missed {M.dn_width.median():.1f} s, accepted {d[~d.missed].dn_width.median():.1f} s')
ro = M[M.cat == 'rt_out']
if len(ro):
    print(f'  rt_out: signed error p5/p50/p95 {ro.err.quantile(.05):.1f}/{ro.err.median():.1f}/{ro.err.quantile(.95):.1f} s; |err|/tol p50 {(ro.err.abs() / tol).median():.2f}')
# C share vs windows per precursor among all DIA-NN IDs
print('  C share (all DIA-NN IDs, accepted or not) by windows per precursor: ' +
      ' '.join(f'n={n}: {d.C[d.n == n].mean():.3f} ({(d.n == n).sum()})' for n in range(1, 6)))
print('  missed share by windows per precursor: ' + ' '.join(f'n={n}: {d.missed[d.n == n].mean():.3f}' for n in range(1, 6)))
M[['cat', 'n', 'near', 'err', 'qq', 'rt_s']].to_csv(f'{RD}/candrec_a_{TAG}.tsv', sep='\t')
