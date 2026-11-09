"""h14a: which already-native feature groups add to search+boundary+byion? (declared before running)
Replica of crossfitStage2LDA (OpenDIAStage2.cpp): per fold, train-row standardisation, 3 rounds, round 1 ranked by cnn_z,
positives = target precursor winners at (D+1)/T <= 0.15 then 0.05 (unpaired, winners over all train windows),
negatives = all decoy winners, ridge 1e-3 LDA. Sanity: the replica on the byion set must reproduce the tool's N
(1029 raw 01 / 1298 raw 02 at FDP_c 1%) within 1%.
Arms: byion (replica), byion+extraction (7 tensor/MS1 features of the -stage2 default set, from the earlier runs with
identical windows), byion+extraction minus ms1_isotope/fragment_fraction (the two near-constant ones).
Rule: an arm PASSES if entrapment-matched N at FDP_c 1% >= 1.03 x replica on BOTH raws and bulk-null q<=0.05 FDP_c
<= 0.06 on both. Usage: python feature_union.py R01|R02"""
import sys, sqlite3, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; R = 148983 / 151239; RUN = sys.argv[1]
BY = {'R01': B + '/astral-design-raw01-20261107/stage2-byion-design/byion', 'R02': B + '/astral-raw02-byion-design-20261109/byion'}[RUN]
EX = {'R01': B + '/astral-design-raw01-20261107/stage2-design/features.tsv', 'R02': B + '/astral-raw02-stage2-check-20261108/crossfit/features.tsv'}[RUN]
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
p = pd.read_sql('select ID precursor_id, TRAML_ID, DECOY from PRECURSOR', c)
p['pair_key'] = np.where(p.DECOY == 1, p.TRAML_ID.str.replace('^DECOY_', '', regex=True), p.TRAML_ID)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c)
EXT = ['fragment_fraction', 'fragment_coelution', 'library_correlation', 'ms2_abs_ppm', 'ms1_abs_ppm', 'ms1_coelution', 'ms1_isotope']
f = pd.read_csv(BY + '/features.tsv', sep='\t')
x = pd.read_csv(EX, sep='\t', usecols=['sample_id', 'precursor_id'] + EXT)
n0 = len(f); f = f.merge(x, on=['sample_id', 'precursor_id']); assert len(f) == n0, 'window mismatch'
f = f.merge(p[['precursor_id', 'pair_key']], on='precursor_id').merge(ent, on='precursor_id', how='left')
BASE = ['cnn_z', 'transition_support', 'log_matched_scans', 'log_peak_width', 'abs_library_rt_error', 'boundary_log_fwhm',
        'boundary_log_width', 'boundary_coelution', 'boundary_library_correlation', 'boundary_apex_agreement', 'byion_b', 'byion_y']
prec = pd.factorize(f.precursor_id)[0]; pair = pd.factorize(f.pair_key)[0]; dec = f.is_decoy.values == 1; fold = f.fold.values
entr = f.e.fillna(0).values.astype(bool)
def winners(idx, s):
    o = idx[np.lexsort((idx, -s[idx]))]; _, first = np.unique(prec[o], return_index=True); return o[first]
def lda_crossfit(cols):
    X = f[cols].values.astype(np.float64); out = np.zeros(len(f))
    for k in np.unique(fold):
        tr = np.where(fold != k)[0]; te = np.where(fold == k)[0]
        mu = X[tr].mean(0); sd = X[tr].std(0); sd[~(sd > 0)] = 1; Z = (X - mu) / sd
        s = X[:, 0].copy()
        for it, fdr in enumerate((0.15, 0.05, 0.05)):
            w = winners(tr, s); w = w[np.lexsort((w, -s[w]))]; d = dec[w]
            T = np.cumsum(~d); Dn = np.cumsum(d); q = np.minimum.accumulate(np.minimum(1, (Dn + 1) / np.maximum(T, 1))[::-1])[::-1]
            pos, neg = w[(~d) & (q <= fdr)], w[d]
            a, b = Z[pos], Z[neg]
            S = (np.cov(a.T) * (len(a) - 1) + np.cov(b.T) * (len(b) - 1)) / (len(a) + len(b) - 2) + 1e-3 * np.eye(len(cols))
            wt = np.linalg.solve(S, a.mean(0) - b.mean(0)); s = Z @ wt
        out[te] = s[te]
    return out, X
def report(name, s, X):
    # pooled: standardise per fold by decoy precursor winners (tool pooled TDC)
    w = winners(np.arange(len(s)), s); z = s.copy()
    for k in np.unique(fold):
        dk = w[dec[w] & (fold[w] == k)]; med = np.median(s[dk]); mad = 1.4826 * np.median(np.abs(s[dk] - med))
        z[fold == k] = (s[fold == k] - med) / mad
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
        dd = d.d.values; T = np.cumsum(~dd); Dn = np.cumsum(np.where(dd, wv[d.i.values], 0.0))
        q = np.minimum.accumulate(np.minimum(1, (Dn + 1) / np.maximum(T, 1))[::-1])[::-1]; kept = d.i.values
        for qq in (0.01, 0.05):
            a = kept[(~dec[kept]) & (q <= qq)]; qs.append(f'{nm} q<={qq}: {len(a)} FDP_c {entr[a].sum() * (1 + 1 / R) / max(1, len(a)):.3f}')
    print(f'{name:34s} N at FDP_c 0.5/1/1.5/2%: {mt} | ' + ' | '.join(qs), flush=True)
print(RUN, len(f), flush=True)
report('tool stage2_score', f.stage2_score.values, f[BASE].values)
for name, cols in (('replica byion', BASE), ('byion+extraction', BASE + EXT),
                   ('byion+extraction-const', BASE + [e for e in EXT if e not in ('ms1_isotope', 'fragment_fraction')])):
    s, X = lda_crossfit(cols); report(name, s, X)
