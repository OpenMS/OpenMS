"""Decoy-quirk check (evaluation only; reads entrapment labels; no score, model or rule changes).

Question: do the scorer's inputs let it tell decoys from false targets, so that decoys stop sampling the false-target
null in the tail (rho ~ 1.5)? If decoys carry construction artifacts that entrapments lack, a feature or a model will
separate decoy winners from entrapment winners at the same score.

Sets per raw (PAIRLIB population, arm's per-precursor winning windows, as evaluate_v2.py):
  E = entrapment target winners, D = decoy winners, H = human target winners.
Measures:
 1. Per input feature (the 20 D11FD columns: 11 heuristics, 6 F12, 3 featD): AUC(D > E) among winners in the tail band
    (standardised score between the q <= 0.20 and q <= 0.01 cuts), and AUC(D > E) score-matched within 10 score bins of
    that band (Mann-Whitney per bin, pooled by bin weights). 0.5 = no tell.
 2. Classifier test: gradient boosting on the 20 features, D vs E winners, 5-fold cross-validated AUC (grouped by
    pair_key family so a decoy and its own target never straddle folds). Repeated with the arm score added (should not
    help if the score is already used) and score-matched (sample weights equalising the score distributions).
    Reference: same classifier for H vs E (human vs entrapment) to show the feature space can separate true from false.
 3. Null control: the same classifier D vs target among null windows (|RT error| > 3 sigma, no peptide expected).
    AUC > 0.5 there is a pure construction tell (no signal involved).
 4. Shared signal (if interf_pre/edges_<raw>.npz exists): share of D, E, H winners with >= 2 shared fragments.
Usage: decoy_quirk.py <root> <arm>"""
import json, sqlite3, sys, os, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score

ROOT, ARM = sys.argv[1:3]
NAMES = ['transition_support', 'log_matched_scans', 'log_peak_width', 'abs_library_rt_error', 'boundary_log_fwhm',
         'boundary_log_width', 'boundary_coelution', 'boundary_library_correlation', 'boundary_apex_agreement', 'byion_b',
         'byion_y', 'f12_support', 'f12_coelution', 'f12_libcos', 'f12_extra_support', 'f12_extra_coelution',
         'f12_extra_share', 'd_max', 'd6_max', 'd_lib']
c = sqlite3.connect(f'file:{ROOT}/lib/library.pqp?mode=ro', uri=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID
  group by pp.PRECURSOR_ID""", c).set_index('precursor_id').e
# The decoy partner class: a decoy of an entrapment vs a decoy of a human target.
pre = pd.read_sql('select ID precursor_id, TRAML_ID traml, DECOY decoy from PRECURSOR', c)
tid = pre[pre.decoy == 0].set_index('traml').precursor_id
pre['partner'] = np.where(pre.decoy == 1, pre.traml.str.replace('^DECOY_', '', regex=True).map(tid), pre.precursor_id)
partner_ent = pre.set_index('precursor_id').partner.map(ent)


def auc(p, n):
    return mannwhitneyu(p, n).statistic / (len(p) * len(n)) if len(p) and len(n) else np.nan


def matched_auc(xp, sp, xn, sn, bins):
    """AUC(p > n) within score bins, weighted by the bin's pair count."""
    num = den = 0.0
    for lo, hi in zip(bins[:-1], bins[1:]):
        a, b = xp[(sp >= lo) & (sp < hi)], xn[(sn >= lo) & (sn < hi)]
        if len(a) and len(b):
            num += auc(a, b) * len(a) * len(b); den += len(a) * len(b)
    return num / den if den else np.nan


def cv_auc(X, y, groups, w=None, seed=0):
    p = np.zeros(len(y))
    for tr, te in GroupKFold(5).split(X, y, groups):
        m = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.05, max_leaf_nodes=15, random_state=seed)
        m.fit(X[tr], y[tr], sample_weight=None if w is None else w[tr]); p[te] = m.predict_proba(X[te])[:, 1]
    return roc_auc_score(y, p, sample_weight=w)


for raw in ('raw01', 'raw02'):
    pop = pd.read_csv(f'{ROOT}/{raw}/population.tsv', sep='\t')
    X = np.c_[np.load(f'{ROOT}/{raw}/heuristics.npy'), np.load(f'{ROOT}/{raw}/feat12.npy'), np.load(f'{ROOT}/{raw}/featD.npy')].astype(np.float64)
    s = np.load(f'{ROOT}/{raw}/scores_{ARM}.npy').astype(np.float64)
    # Winners and q as evaluate_v2.py.
    w = pop.assign(s=s, r=np.arange(len(pop))).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').copy()
    for f_, g in w[w.is_decoy == 1].groupby('fold'):
        med = g.s.median(); mad = 1.4826 * (g.s - med).abs().median()
        w.loc[w.fold == f_, 's'] = (w.loc[w.fold == f_, 's'] - med) / mad
    best = w.sort_values(['s', 'is_decoy'], ascending=[False, False]).drop_duplicates('pair_key').sort_values('s', ascending=False)
    T = np.cumsum(1 - best.is_decoy.values); D = np.cumsum(best.is_decoy.values)
    best['q'] = np.minimum.accumulate(np.minimum(1, (D + 1) / np.maximum(T, 1))[::-1])[::-1]
    cut = {q: best[best.q <= q].s.min() for q in (0.01, 0.05, 0.20)}
    # Classes over all precursor winners (not only pair winners: a decoy that lost to its target is still a decoy winner).
    w['cls'] = np.where(w.is_decoy == 1, 'D', np.where(w.precursor_id.map(ent) == 1, 'E', 'H'))
    w['dec_of_ent'] = w.precursor_id.map(partner_ent)
    band = w[(w.s >= cut[0.20])]
    print(f'\n==== {raw} {ARM}: cuts q0.01 {cut[0.01]:.2f} q0.05 {cut[0.05]:.2f} q0.20 {cut[0.20]:.2f}; winners above the q<=0.20 cut:',
          band.cls.value_counts().to_dict())
    Dw, Ew, Hw = band[band.cls == 'D'], band[band.cls == 'E'], band[band.cls == 'H']
    bins = np.quantile(pd.concat([Dw.s, Ew.s]), np.linspace(0, 1, 11)); bins[-1] += 1e-9
    print('1. per feature AUC(decoy > entrapment): band | score-matched (10 bins) | for reference AUC(human > entrapment) band')
    rows = []
    for k, nm in enumerate(NAMES):
        xd, xe, xh = X[Dw.r.values, k], X[Ew.r.values, k], X[Hw.r.values, k]
        rows.append((nm, auc(xd, xe), matched_auc(xd, Dw.s.values, xe, Ew.s.values, bins), auc(xh, xe)))
    for nm, a, am, ah in sorted(rows, key=lambda t: -abs(t[2] - 0.5)):
        print(f'   {nm:30s} {a:.3f} | {am:.3f} | {ah:.3f}')
    # 2. Classifier D vs E (and H vs E) on the band, grouped folds by pair family.
    sub = pd.concat([Dw, Ew]); y = (sub.cls == 'D').values.astype(int); grp = sub.pair_key.values
    Xs = X[sub.r.values]
    a0 = cv_auc(Xs, y, grp)
    a1 = cv_auc(np.c_[Xs, sub.s.values], y, grp)
    # Score-matched weights: inverse density of each class within the 10 score bins.
    b = np.clip(np.searchsorted(bins, sub.s.values, side='right') - 1, 0, 9)
    cnt = pd.crosstab(b, y); wts = np.array([1.0 / cnt.loc[bi, yi] if cnt.loc[bi, yi] else 0 for bi, yi in zip(b, y)])
    keep = np.array([cnt.loc[bi].min() > 0 for bi in b]); a2 = cv_auc(Xs[keep], y[keep], grp[keep], wts[keep])
    subh = pd.concat([Hw.sample(min(len(Hw), 4 * len(Ew)), random_state=0), Ew]); yh = (subh.cls == 'H').values.astype(int)
    ah = cv_auc(X[subh.r.values], yh, subh.pair_key.values)
    print(f'2. classifier CV AUC decoy vs entrapment winners (n D {len(Dw)}, E {len(Ew)}): features {a0:.3f} | + score {a1:.3f} | '
          f'score-matched {a2:.3f};  human vs entrapment (reference) {ah:.3f}')
    # Split decoys by partner class: decoys of entrapments vs decoys of human targets.
    for lab, m in (('decoys of human targets', Dw.dec_of_ent == 0), ('decoys of entrapments', Dw.dec_of_ent == 1)):
        dd = Dw[m.values]; ss = pd.concat([dd, Ew]); yy = (ss.cls == 'D').values.astype(int)
        print(f'   {lab}: n {len(dd)}, CV AUC vs entrapment winners {cv_auc(X[ss.r.values], yy, ss.pair_key.values):.3f}')
    # 3. Null control: windows far from library RT, target (any) vs decoy.
    sig = json.load(open(f'{ROOT}/{raw}/cal/cal.json'))['rt']['sigma_s']
    null = np.where(X[:, 3] > 3 * sig)[0]
    rng = np.random.default_rng(0); null = rng.choice(null, min(len(null), 60000), replace=False)
    yn = pop.is_decoy.values[null]; gn = pop.pair_key.values[null]
    print(f'3. null windows (|RT err| > 3 sigma, n {len(null)}): classifier CV AUC decoy vs target {cv_auc(X[null], yn, gn):.3f}')
    for k in np.argsort([-abs(auc(X[null][yn == 1, k], X[null][yn == 0, k]) - 0.5) for k in range(len(NAMES))])[:5]:
        print(f'   {NAMES[k]:30s} AUC(decoy > target) {auc(X[null][yn == 1, k], X[null][yn == 0, k]):.3f}')
    # 4. Shared signal.
    ef = f'{ROOT}/interf_pre/edges_{raw}.npz'
    if os.path.exists(ef):
        e = np.load(ef); m2 = e['w'] >= 2
        shared = set(e['a'][m2]) | set(e['b'][m2])
        for nm, g in (('D', Dw), ('E', Ew), ('H', Hw)):
            print(f'4. shared signal (>= 2 fragments) {nm}: {np.mean(g.sample_id.isin(shared)):.3f}')
