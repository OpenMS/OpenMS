"""GAPDIAG (diagnostic): what does DIA-NN see in the IDs we miss? r~7 HeLa D seed 0, current default.
DIA-NN IDs (q<=0.01, non-entrapment) are split by what OpenDIA did with them:
  A   we identify it (our run q<=0.01)
  B1  not identified; our best window is the right peak (|rt - DIA-NN RT| <= 6 s): right peak, scored too low
  B2  not identified; a right window exists but our best window is another peak: wrong peak chosen
  C   no window within 6 s of DIA-NN's RT: no candidate
Then: (1) DIA-NN's own report fields by class (is B weak for DIA-NN too, or clean?),
(2) our features on the right window, per class, against the decoys' best windows (which of our features fail),
(3) agreement of our extraction with DIA-NN (quantity, MS1, RT, IM) per class,
(4) where B1 sits in our q ranking. Labels (DIA-NN, entrapment) are for evaluation only.
Usage: gap_diag.py <features.tsv> <run.precursors.tsv> <calibration.json> <lib.pqp> <diann report.parquet>"""
import sys, json, sqlite3, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr
FT, PT, CJ, L, DN = sys.argv[1:6]
W = 6.0
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40)

c = sqlite3.connect(f'file:{L}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID pid, p.DECOY d, p.CHARGE z, p.PRECURSOR_MZ mz, p.LIBRARY_RT lrt, p.LIBRARY_DRIFT_TIME lim,
  pe.MODIFIED_SEQUENCE ms, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c).set_index('pid')
tg = pq[pq.d == 0]
d = pd.read_parquet(DN)
d = d[(d['Q.Value'] <= 0.01) & (d.Decoy == 0)].copy()
sk = (tg.ms + '/' + tg.z.astype(str)).reset_index().drop_duplicates(0).set_index(0).pid
d['pid'] = sk.reindex(d['Modified.Sequence'] + '/' + d['Precursor.Charge'].astype(str)).values
print(f'DIA-NN q<=0.01 targets {len(d)}, matched to library {d.pid.notna().sum()}')
d = d.dropna(subset=['pid']); d['pid'] = d.pid.astype(np.int64)
d = d[pq.e.reindex(d.pid).values == 0].drop_duplicates('pid').set_index('pid')
d['rt_s'] = d.RT * 60.0
print(f'DIA-NN non-entrapment IDs used: {len(d)}')

meta = {'candidate_row': np.int64, 'precursor_id': np.int64, 'is_decoy': np.int8, 'fold': np.int8, 'rt_apex': np.float64, 'q_value': np.float64}
head = pd.read_csv(FT, sep='\t', nrows=0).columns
F = pd.read_csv(FT, sep='\t', engine='c', dtype={k: (meta[k] if k in meta else np.float32) for k in head})
feat = [x for x in F.columns if x not in meta and x != 'score']
print(f'windows {len(F)}, features {len(feat)}', flush=True)
P = pd.read_csv(PT, sep='\t', usecols=['candidate_row', 'precursor_id', 'is_decoy', 'rt_apex', 'score', 'q_value'])
P = P[P.score > -1e29].sort_values('score', ascending=False).drop_duplicates('precursor_id').set_index('precursor_id')
ours = P[(P.is_decoy == 0) & (P.q_value <= 0.01)].index
print(f'our q<=0.01 targets {len(ours)}; of them in DIA-NN set {np.isin(ours, d.index).sum()}')

# Right window per DIA-NN precursor: the window of that precursor closest to DIA-NN's RT, if within W.
Ft = F[F.precursor_id.isin(d.index)].copy()
Ft['drt'] = (Ft.rt_apex - d.rt_s.reindex(Ft.precursor_id).values).abs()
right = Ft.sort_values('drt').drop_duplicates('precursor_id').set_index('precursor_id')
right = right[right.drt <= W]
best_row = P.candidate_row
cls = pd.Series('C', index=d.index)
has_right = d.index.isin(right.index)
cls[has_right] = 'B2'
br = best_row.reindex(d.index)
cls[has_right & (br.values == right.candidate_row.reindex(d.index).values)] = 'B1'
cls[d.index.isin(ours)] = 'A'
d['cls'] = cls
print('\n== classes of DIA-NN IDs')
print(cls.value_counts().reindex(['A', 'B1', 'B2', 'C']).to_frame('n').assign(share=lambda x: (x.n / x.n.sum()).round(3)))
# A split: does our ID sit on DIA-NN's peak?
a = d.index[cls == 'A']
print(f'A on the same peak (|our rt - DIA-NN rt| <= {W} s): {(np.abs(P.rt_apex.reindex(a) - d.rt_s.reindex(a)) <= W).mean():.3f}')

# (1) DIA-NN's own fields by class
cal = json.load(open(CJ))
im_a, im_b = cal['im']['intercept_vssc'], cal['im']['slope']
kr, ks = np.array(cal['rt']['knots_library_rt']), np.array(cal['rt']['knots_run_rt_s'])
lrt = pq.lrt.reindex(d.index).values; lim = pq.lim.reindex(d.index).values
d['our_pred_rt_s'] = np.interp(lrt, kr, ks)
d['dn_rt_err_vs_our_map_s'] = (d.rt_s - d.our_pred_rt_s).abs()
d['dn_im_err_vs_our_cal'] = (d.IM - (im_a + im_b * lim)).abs()
d['log10_quantity'] = np.log10(d['Precursor.Quantity'].clip(lower=1))
d['log10_ms1_area'] = np.log10(d['Ms1.Area'].clip(lower=1))
d['fwhm_s'] = d.FWHM * 60.0
d['dn_rt_err_own_s'] = (d.RT - d['Predicted.RT']).abs() * 60.0
d['dn_im_err_own'] = (d.IM - d['Predicted.IM']).abs()
d['length'] = d['Stripped.Sequence'].str.len()
dnf = ['log10_quantity', 'log10_ms1_area', 'Evidence', 'Ms1.Profile.Corr', 'Quantity.Quality', 'PEP', 'fwhm_s',
       'dn_rt_err_own_s', 'dn_rt_err_vs_our_map_s', 'dn_im_err_own', 'dn_im_err_vs_our_cal', 'Precursor.Mz', 'Precursor.Charge', 'length']
print('\n== (1) DIA-NN report fields by class: median [AUC of A vs B1, A vs C; 0.5 = no difference]')
rows = []
for f_ in dnf:
    r = {'field': f_}
    for k in ['A', 'B1', 'B2', 'C']: r[k] = d.loc[d.cls == k, f_].median()
    for k in ['B1', 'C']:
        x = d.loc[d.cls.isin(['A', k]), [f_, 'cls']].dropna()
        r[f'AUC A>{k}'] = roc_auc_score(x.cls == 'A', x[f_]) if x.cls.nunique() == 2 else np.nan
    rows.append(r)
print(pd.DataFrame(rows).set_index('field').round(4).to_string())
print('\n== share of each class by DIA-NN quantity decile (1 = faintest)')
d['qdec'] = pd.qcut(d.log10_quantity, 10, labels=False) + 1
print(pd.crosstab(d.qdec, d.cls, normalize='index').round(3).to_string())

# (2) our features on the right window per class, against the decoys' best windows
dec_rows = P[P.is_decoy == 1].candidate_row
DW = F.set_index('candidate_row').loc[dec_rows.values, feat + ['score']]
R = right.join(d.cls)
print(f'\n== (2) our features: right window per class vs decoy best windows (n decoys {len(DW)})')
print('   AUC(class vs decoys); low AUC for B1 where A is high = feature that fails on the missed IDs')
rows = []
for f_ in feat + ['score']:
    r = {'feature': f_}
    dv = DW[f_].values.astype(np.float64); dv = dv[np.isfinite(dv)]
    for k in ['A', 'B1', 'B2']:
        v = R.loc[R.cls == k, f_].values.astype(np.float64); v = v[np.isfinite(v)]
        r[f'med {k}'] = np.median(v) if len(v) else np.nan
        if len(v) > 20:
            rng = np.random.default_rng(0); dd = rng.choice(dv, min(len(dv), 200000), replace=False)
            r[f'AUC {k}'] = roc_auc_score(np.r_[np.ones(len(v)), np.zeros(len(dd))], np.r_[v, dd])
    r['med dec'] = np.median(dv)
    rows.append(r)
T = pd.DataFrame(rows).set_index('feature')
T['A-B1 gap'] = (T['AUC A'] - 0.5).abs() - (T['AUC B1'] - 0.5).abs()
print(T.round(3).sort_values('A-B1 gap', ascending=False).to_string())

# (3) agreement of our extraction with DIA-NN, per class, on the right window
print('\n== (3) agreement with DIA-NN on the right window (Spearman per class)')
R = R.join(d[['log10_quantity', 'log10_ms1_area', 'rt_s', 'IM', 'Ms1.Profile.Corr', 'Evidence']])
for ours_f, dn_f in [('x_ms2_log_area', 'log10_quantity'), ('x_ms1_log_area', 'log10_ms1_area'), ('x_ms1_coelution', 'Ms1.Profile.Corr'),
                     ('x_fragment_coelution', 'Evidence'), ('transition_support', 'Evidence')]:
    if ours_f not in R: continue
    s = {k: spearmanr(R.loc[R.cls == k, ours_f], R.loc[R.cls == k, dn_f], nan_policy='omit')[0] for k in ['A', 'B1', 'B2']}
    print(f'  {ours_f:22s} vs {dn_f:18s} ' + '  '.join(f'{k} {v:.3f}' for k, v in s.items()))
R['rt_diff'] = R.rt_apex - R.rt_s
print('  our right-window apex - DIA-NN RT (s), median / p90 |.|: ' + '  '.join(
    f'{k} {R.loc[R.cls == k, "rt_diff"].median():+.2f} / {R.loc[R.cls == k, "rt_diff"].abs().quantile(0.9):.2f}' for k in ['A', 'B1', 'B2']))
print(f'  our IM tolerance used: +-{cal["im"]["tolerance_half_width_vssc"] * 0.6:.4f} (0.6 x calibrated)')
for k in ['A', 'B1', 'B2', 'C']:
    e = d.loc[d.cls == k, 'dn_im_err_vs_our_cal']
    print(f'  {k}: DIA-NN IM vs our calibrated library IM, median {e.median():.4f}, p90 {e.quantile(0.9):.4f}, share outside our window {(e > cal["im"]["tolerance_half_width_vssc"] * 0.6).mean():.3f}')
for k in ['A', 'B1', 'B2', 'C']:
    e = d.loc[d.cls == k, 'dn_rt_err_vs_our_map_s']
    print(f'  {k}: DIA-NN RT vs our RT map, median {e.median():.1f} s, p90 {e.quantile(0.9):.1f} s, share outside our {cal["rt"]["tolerance_half_width_s"]:.0f} s window {(e > cal["rt"]["tolerance_half_width_s"]).mean():.3f}')

# (4) where B1 sits in our ranking
print('\n== (4) our q of B1 (best window = right peak)')
qb = P.q_value.reindex(d.index[d.cls == 'B1'])
print(pd.cut(qb, [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0001]).value_counts(sort=False).to_string())
print('\n== (5) B2: our best window vs the right one')
b2 = d.index[d.cls == 'B2']
bw = F.set_index('candidate_row').loc[best_row.reindex(b2).values]
rw = right.loc[b2]
print(f'  score best (wrong) median {bw.score.median():.3f}  right median {rw.score.median():.3f}; |rt best - DIA-NN| median {np.median(np.abs(bw.rt_apex.values - d.rt_s.reindex(b2).values)):.1f} s')
for f_ in ['transition_support', 'x_fragment_coelution', 'x_library_correlation', 'abs_library_rt_error', 'x_ms2_log_area', 'x_im_centroid_offset']:
    if f_ in F: print(f'  {f_:24s} wrong {np.nanmedian(bw[f_].values):.3f}  right {np.nanmedian(rw[f_].values):.3f}')
d[['Precursor.Id', 'cls']].to_csv('classes_D.tsv', sep='\t', index=False)  # classes.tsv
print('done')
