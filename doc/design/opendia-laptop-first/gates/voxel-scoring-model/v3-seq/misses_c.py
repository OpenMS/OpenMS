"""D11 misses, part c (descriptive): (1) direction of the D11 score per feature (per fold, least-squares fit of the score on the
z-scored 11 features over all windows, and the median score by transition_support); (2) the entrapment target winners near the
FDP_c 1 % cut (they set the cut)."""
import sqlite3, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; O = B + '/astral-voxel-pilot-20261110/v2'; R = 148983 / 151239
FEATS = ['transition_support', 'log_matched_scans', 'log_peak_width', 'abs_library_rt_error', 'boundary_log_fwhm', 'boundary_log_width',
         'boundary_coelution', 'boundary_library_correlation', 'boundary_apex_agreement', 'byion_b', 'byion_y']
BASEF = {'raw01': B + '/astral-design-raw01-20261107/stage2-byion-design/byion/features.tsv', 'raw02': B + '/astral-raw02-byion-design-20261109/byion/features.tsv'}
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c).set_index('precursor_id').e
for raw in ('raw01', 'raw02'):
    pop = pd.read_csv(f'{O}/{raw}/population.tsv', sep='\t'); H = np.load(f'{O}/{raw}/heuristics.npy').astype(np.float64)
    s = np.load(f'{O}/{raw}/scores_D11.npy').astype(np.float64)
    base = pd.read_csv(BASEF[raw], sep='\t', usecols=['sample_id', 'stage2_score']).set_index('sample_id').stage2_score.loc[pop.sample_id].values
    print(f'\n=== {raw}: linear fit of the D11 score on z-scored features (per fold; R^2; coefficient per 1 SD)')
    for k in range(3):
        m = pop.fold.values == k; Z = (H[m] - H[m].mean(0)) / H[m].std(0); X = np.c_[np.ones(m.sum()), Z]
        coef, *_ = np.linalg.lstsq(X, s[m], rcond=None); r2 = 1 - ((s[m] - X @ coef) ** 2).sum() / ((s[m] - s[m].mean()) ** 2).sum()
        cb, *_ = np.linalg.lstsq(X, base[m], rcond=None)
        print(f'  fold {k} R2 {r2:.3f}  D11: ' + ' '.join(f'{f[:10]}={v:+.1f}' for f, v in zip(FEATS, coef[1:])))
        print(f'         BASE (same fit):  ' + ' '.join(f'{f[:10]}={v:+.2f}' for f, v in zip(FEATS, cb[1:])))
    d = pd.DataFrame({'ts': H[:, 0], 'rt': H[:, 3], 's': s, 'dec': pop.is_decoy.values})
    print('  median D11 score by transition_support (targets):', d[d.dec == 0].groupby('ts').s.median().round(1).to_dict())
    d['rtbin'] = pd.qcut(d.rt, 10, labels=False, duplicates='drop')
    print('  median D11 score by abs_library_rt_error decile (targets):', d[d.dec == 0].groupby('rtbin').s.median().round(1).tolist())
    # Entrapment near the cut.
    w = pop.assign(s=s, row=np.arange(len(pop))).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').copy()
    for f_, g in w[w.is_decoy == 1].groupby('fold'):
        med = g.s.median(); mad = 1.4826 * (g.s - med).abs().median(); w.loc[w.fold == f_, 's'] = (w.loc[w.fold == f_, 's'] - med) / mad
    t = w[w.is_decoy == 0].assign(e=lambda x: x.precursor_id.map(ent).values).sort_values('s', ascending=False).reset_index(drop=True)
    t['N'] = np.arange(1, len(t) + 1); t['ne'] = t.e.cumsum(); t['fdp'] = t.ne * (1 + 1 / R) / t.N
    er = t[t.e == 1].head(12)
    print('  first 12 entrapment target winners: rank N, FDP_c there, fold, score, features (ts, rt_err, coel, libcorr, apex, b, y)')
    for _, r in er.iterrows():
        h = H[int(r.row)]
        print(f'   N {int(r.N):5d} fdp {r.fdp:.4f} fold {int(r.fold)} s {r.s:8.1f}  ts {h[0]:.0f} rt {h[3]:5.1f} coel {h[6]:.2f} lib {h[7]:.2f} apex {h[8]:.2f} b {h[9]:.0f} y {h[10]:.0f}')
