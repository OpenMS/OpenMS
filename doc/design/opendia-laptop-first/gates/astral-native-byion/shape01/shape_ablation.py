"""h14b: which shape subgroup brings the IDs and which breaks 5% calibration? (raw 01, declared before running)
Exact replica of the tool LDA (validated: reproduces stage2_score on the byion set). Arms add to the byion base:
  xcorr   = shape_xcorr_* (5)
  library = shape_library_* (4)
  signal  = shape_log_sn, shape_intensity_share, shape_signal_fragments, shape_min_fragment_share
  all     = all 13 (the tool run)
Evaluation: eval_arm metric (decoy-standardised score), entrapment-matched N at FDP_c 1%, plain and bulk-null (seed 0)
q<=0.05 FDP_c. Reading: a subgroup is a candidate if N >= 1.03 x base (1010) and bulk-null q<=0.05 FDP_c <= 0.06.
Diagnostic for design only; any candidate must then be run natively on both raws."""
import sqlite3, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; R = 148983 / 151239
ARM = B + '/astral-design-raw01-20261107/stage2-shape-design/shape'
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
p = pd.read_sql('select ID precursor_id, TRAML_ID, DECOY from PRECURSOR', c)
p['pair_key'] = np.where(p.DECOY == 1, p.TRAML_ID.str.replace('^DECOY_', '', regex=True), p.TRAML_ID)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c)
f = pd.read_csv(ARM + '/features.tsv', sep='\t').merge(p[['precursor_id', 'pair_key']], on='precursor_id').merge(ent, on='precursor_id', how='left')
BASE = ['cnn_z', 'transition_support', 'log_matched_scans', 'log_peak_width', 'abs_library_rt_error', 'boundary_log_fwhm',
        'boundary_log_width', 'boundary_coelution', 'boundary_library_correlation', 'boundary_apex_agreement', 'byion_b', 'byion_y']
SH = [x for x in f.columns if x.startswith('shape_')]
GROUPS = {'xcorr': [x for x in SH if 'xcorr' in x], 'library': [x for x in SH if 'library' in x],
          'signal': ['shape_log_sn', 'shape_intensity_share', 'shape_signal_fragments', 'shape_min_fragment_share'], 'all': SH}
prec = pd.factorize(f.precursor_id)[0]; pair = pd.factorize(f.pair_key)[0]; dec = f.is_decoy.values == 1; fold = f.fold.values
entr = f.e.fillna(0).values.astype(bool)
def winners(idx, s):
    o = idx[np.lexsort((idx, -s[idx]))]; _, first = np.unique(prec[o], return_index=True); return o[first]
def lda(cols):
    X = f[cols].values.astype(np.float64); out = np.zeros(len(f))
    for k in np.unique(fold):
        tr = np.where(fold != k)[0]; te = np.where(fold == k)[0]
        mu = X[tr].mean(0); sd = X[tr].std(0); sd[~(sd > 0)] = 1; Z = (X - mu) / sd; s = X[:, 0].copy()
        for fdr in (0.15, 0.05, 0.05):
            w = winners(tr, s); w = w[np.lexsort((w, -s[w]))]; d = dec[w]
            q = np.minimum.accumulate(np.minimum(1, (np.cumsum(d) + 1) / np.maximum(np.cumsum(~d), 1))[::-1])[::-1]
            a, b = Z[w[(~d) & (q <= fdr)]], Z[w[d]]
            S = (np.cov(a.T) * (len(a) - 1) + np.cov(b.T) * (len(b) - 1)) / (len(a) + len(b) - 2) + 1e-3 * np.eye(len(cols))
            s = Z @ np.linalg.solve(S, a.mean(0) - b.mean(0))
        out[te] = s[te]
    return out, X
def evaluate(name, s, X):
    w = winners(np.arange(len(s)), s); z = s.copy()
    for k in np.unique(fold):
        dk = w[dec[w] & (fold[w] == k)]; med = np.median(s[dk]); z[fold == k] = (s[fold == k] - med) / (1.4826 * np.median(np.abs(s[dk] - med)))
    w = winners(np.arange(len(z)), z); tg = w[~dec[w]]; o = tg[np.argsort(-z[tg], kind='stable')]
    fd = np.cumsum(entr[o]) * (1 + 1 / R) / np.arange(1, len(o) + 1)
    mt = [int(np.where(fd <= t)[0].max() + 1) if (fd <= t).any() else 0 for t in (0.005, 0.01, 0.015, 0.02)]
    wdv = w[dec[w]]; twv = w[~dec[w]]; thr = np.median(z[wdv]); bt, bd = twv[z[twv] < thr], wdv[z[wdv] < thr]
    r = np.ones(len(z)); pr = np.zeros(len(z))
    for k in np.unique(fold):
        trn = np.r_[bt[fold[bt] != k], bd[fold[bd] != k]]
        m = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.05, max_leaf_nodes=31, random_state=0).fit(X[trn], (~dec[trn]).astype(int))
        te = wdv[fold[wdv] == k]; pr[te] = m.predict_proba(X[te])[:, 1]
    r[wdv] = np.clip(pr[wdv] / np.maximum(1 - pr[wdv], 1e-9) * len(bd) / len(bt), 0.2, 5); r[wdv] /= r[bd].mean()
    qs = []
    for nm, wv in (('plain', np.ones(len(z))), ('bulk-null', r)):
        d = pd.DataFrame({'i': w, 'p': pair[w], 's': z[w], 'd': dec[w]}).sort_values(['s', 'd'], ascending=[False, False]).drop_duplicates('p').sort_values('s', ascending=False)
        dd = d.d.values; Dn = np.cumsum(np.where(dd, wv[d.i.values], 0.0))
        q = np.minimum.accumulate(np.minimum(1, (Dn + 1) / np.maximum(np.cumsum(~dd), 1))[::-1])[::-1]; kept = d.i.values
        for qq in (0.01, 0.05):
            a = kept[(~dec[kept]) & (q <= qq)]; qs.append(f'{nm} q<={qq}: {len(a)} ({entr[a].sum() * (1 + 1 / R) / max(1, len(a)):.3f})')
    print(f'{name:16s} N@FDP_c 0.5/1/1.5/2%: {mt} | ' + ' | '.join(qs), flush=True)
evaluate('tool +shape', f.stage2_score.values.astype(np.float64), f[BASE + SH].values)
s, X = lda(BASE); evaluate('base (replica)', s, X)
for g, cols in GROUPS.items():
    s, X = lda(BASE + cols); evaluate('+' + g, s, X)
