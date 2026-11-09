"""k7: do library-agreement features drive the 5% decoy-q excess? (declared before running)
Exact tool-LDA replica on the native search+boundary+byion features (raws 01/02). Arms:
  base           : all 12 features
  -libcorr       : minus boundary_library_correlation
  -libcorr-apex  : minus boundary_library_correlation and boundary_apex_agreement
  -lib-any       : minus the above and abs_library_rt_error (control: an unrelated library property)
Metric: entrapment-matched N at FDP_c 1% (decoy-standardised score) and plain pooled paired q<=0.01/0.05 FDP_c.
Rule: SUPPORTED if -libcorr (or -libcorr-apex) gives plain q<=0.05 FDP_c <= 0.06 on BOTH raws, or lowers it by >= 0.01 on
both, while N at FDP_c 1% falls by < 10%. REJECTED if q<=0.05 FDP_c changes by < 0.005 on either raw."""
import sqlite3, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; R = 148983 / 151239
ARMS = {'raw01': B + '/astral-design-raw01-20261107/stage2-byion-design/byion', 'raw02': B + '/astral-raw02-byion-design-20261109/byion'}
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
p = pd.read_sql('select ID precursor_id, TRAML_ID, DECOY from PRECURSOR', c)
p['pair_key'] = np.where(p.DECOY == 1, p.TRAML_ID.str.replace('^DECOY_', '', regex=True), p.TRAML_ID)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c)
BASE = ['cnn_z', 'transition_support', 'log_matched_scans', 'log_peak_width', 'abs_library_rt_error', 'boundary_log_fwhm',
        'boundary_log_width', 'boundary_coelution', 'boundary_library_correlation', 'boundary_apex_agreement', 'byion_b', 'byion_y']
DROPS = {'base': [], '-libcorr': ['boundary_library_correlation'], '-libcorr-apex': ['boundary_library_correlation', 'boundary_apex_agreement'],
         '-lib-any': ['boundary_library_correlation', 'boundary_apex_agreement', 'abs_library_rt_error']}
for run, arm in ARMS.items():
    f = pd.read_csv(arm + '/features.tsv', sep='\t').merge(p[['precursor_id', 'pair_key']], on='precursor_id').merge(ent, on='precursor_id', how='left')
    prec = pd.factorize(f.precursor_id)[0]; pair = pd.factorize(f.pair_key)[0]; dec = f.is_decoy.values == 1; fold = f.fold.values
    entr = f.e.fillna(0).values.astype(bool)
    def winners(idx, s):
        o = idx[np.lexsort((idx, -s[idx]))]; _, first = np.unique(prec[o], return_index=True); return o[first]
    for name, drop in DROPS.items():
        cols = [x for x in BASE if x not in drop]; X = f[cols].values.astype(np.float64); s_out = np.zeros(len(f))
        for k in np.unique(fold):
            tr = np.where(fold != k)[0]; te = np.where(fold == k)[0]
            mu = X[tr].mean(0); sd = X[tr].std(0); sd[~(sd > 0)] = 1; Z = (X - mu) / sd; s = X[:, 0].copy()
            for fdr in (0.15, 0.05, 0.05):
                w = winners(tr, s); w = w[np.lexsort((w, -s[w]))]; d = dec[w]
                q = np.minimum.accumulate(np.minimum(1, (np.cumsum(d) + 1) / np.maximum(np.cumsum(~d), 1))[::-1])[::-1]
                a, b = Z[w[(~d) & (q <= fdr)]], Z[w[d]]
                S = (np.cov(a.T) * (len(a) - 1) + np.cov(b.T) * (len(b) - 1)) / (len(a) + len(b) - 2) + 1e-3 * np.eye(len(cols))
                s = Z @ np.linalg.solve(S, a.mean(0) - b.mean(0))
            s_out[te] = s[te]
        s = s_out; w = winners(np.arange(len(s)), s); z = s.copy()
        for k in np.unique(fold):
            dk = w[dec[w] & (fold[w] == k)]; med = np.median(s[dk]); z[fold == k] = (s[fold == k] - med) / (1.4826 * np.median(np.abs(s[dk] - med)))
        w = winners(np.arange(len(z)), z); tg = w[~dec[w]]; o = tg[np.argsort(-z[tg], kind='stable')]
        fd = np.cumsum(entr[o]) * (1 + 1 / R) / np.arange(1, len(o) + 1)
        n1 = int(np.where(fd <= 0.01)[0].max() + 1)
        d0 = pd.DataFrame({'i': w, 'p': pair[w], 's': z[w], 'd': dec[w]}).sort_values(['s', 'd'], ascending=[False, False]).drop_duplicates('p').sort_values('s', ascending=False)
        dd = d0.d.values; q = np.minimum.accumulate(np.minimum(1, (np.cumsum(dd) + 1) / np.maximum(np.cumsum(~dd), 1))[::-1])[::-1]; kept = d0.i.values
        qs = []
        for qq in (0.01, 0.05):
            a = kept[(~dec[kept]) & (q <= qq)]; qs.append(f'q<={qq}: {len(a)} FDP_c {entr[a].sum() * (1 + 1 / R) / max(1, len(a)):.3f}')
        print(f'{run} {name:14s} N@FDP_c1% {n1} | ' + ' | '.join(qs), flush=True)
