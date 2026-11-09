"""h13: does an MLP head on the NATIVE stage-2 features fix 5% calibration and keep N? (declared before running)
Input: <arm>/features.tsv of the native byion run (all in-window windows, 12 features, tool fold, tool LDA stage2_score).
Head: per tool fold k, train on folds != k. Start from the tool's out-of-fold LDA score; 2 rounds of: positives = target
precursor winners at paired-TDC q<=0.02 (>= 200), negatives = all decoy precursor winners; MLP 32-16 GELU (oswlike_score
mlp_fit, standardised features). Held-out fold scored and standardised by its decoy winners (median / 1.4826 MAD).
Ensemble = mean of the standardised scores over seeds 0, 10, 20, 30, 40. Then bulk-null weighting (bulknull_native method).
Rule: PASS if on BOTH raws (a) entrapment-matched N at FDP_c 1% >= the tool LDA's (1029 / 1298) and (b) bulk-null
q<=0.05 FDP_c <= 0.06 and q<=0.01 FDP_c <= 0.015. Usage: python mlp_head_native.py <arm dir>"""
import sys, sqlite3, numpy as np, pandas as pd, torch, torch.nn as nn
from sklearn.ensemble import HistGradientBoostingClassifier
torch.set_num_threads(12)
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; R = 148983 / 151239; ARM = sys.argv[1]
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
p = pd.read_sql('select ID precursor_id, TRAML_ID, DECOY from PRECURSOR', c)
p['pair_key'] = np.where(p.DECOY == 1, p.TRAML_ID.str.replace('^DECOY_', '', regex=True), p.TRAML_ID)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c)
f = pd.read_csv(ARM + '/features.tsv', sep='\t').merge(p[['precursor_id', 'pair_key']], on='precursor_id').merge(ent, on='precursor_id', how='left')
FEAT = [x for x in f.columns if x not in ('sample_id', 'precursor_id', 'is_decoy', 'fold', 'stage2_score', 'pair_key', 'e')]
X = f[FEAT].values.astype(np.float64); Z = np.nan_to_num((X - X.mean(0)) / (X.std(0) + 1e-9)); Zt = torch.tensor(Z, dtype=torch.float32)
prec = pd.factorize(f.precursor_id)[0]; pair = pd.factorize(f.pair_key)[0]; dec = f.is_decoy.values == 1; fold = f.fold.values
entr = f.e.fillna(0).values.astype(bool); sid = f.sample_id.values
def winners(idx, s):
    o = idx[np.lexsort((sid[idx], -s[idx]))]; _, first = np.unique(prec[o], return_index=True); return o[first]
def tdc(w, s, wdec=None):
    wd = np.ones(len(s)) if wdec is None else wdec
    d = pd.DataFrame({'i': w, 'p': pair[w], 's': s[w], 'd': dec[w]}).sort_values(['s', 'd'], ascending=[False, False]).drop_duplicates('p').sort_values('s', ascending=False)
    dd = d.d.values; T = np.cumsum(~dd); Dn = np.cumsum(np.where(dd, wd[d.i.values], 0.0))
    return d.i.values, np.minimum.accumulate(np.minimum(1, (Dn + 1) / np.maximum(T, 1))[::-1])[::-1]
def mlp_fit(pos, neg, seed):
    torch.manual_seed(seed); idx = np.r_[pos, neg]; y = torch.tensor(np.r_[np.ones(len(pos)), np.zeros(len(neg))], dtype=torch.float32)
    net = nn.Sequential(nn.Linear(Z.shape[1], 32), nn.GELU(), nn.Linear(32, 16), nn.GELU(), nn.Linear(16, 1))
    opt = torch.optim.AdamW(net.parameters(), lr=3e-3, weight_decay=1e-3); lf = nn.BCEWithLogitsLoss(pos_weight=torch.tensor(len(neg) / max(1, len(pos))))
    for _ in range(5):
        perm = np.random.default_rng(seed).permutation(len(idx))
        for b in range(0, len(idx), 4096):
            j = perm[b:b + 4096]; loss = lf(net(Zt[idx[j]]).squeeze(1), y[j]); opt.zero_grad(); loss.backward(); opt.step()
    with torch.no_grad():
        return np.concatenate([net(Zt[i:i + 500000]).squeeze(1).numpy() for i in range(0, len(Zt), 500000)]).astype(np.float64)
def head(seed):
    out = np.zeros(len(Z))
    for k in np.unique(fold):
        tr = np.where(fold != k)[0]; s = f.stage2_score.values.copy()
        for it in range(2):
            w = winners(tr, s); kept, q = tdc(w, s); pos = kept[(~dec[kept]) & (q <= 0.02)]
            if len(pos) < 200: tw = w[~dec[w]]; pos = tw[np.argsort(-s[tw], kind='stable')[:200]]
            s = mlp_fit(pos, w[dec[w]], seed + it)
        te = np.where(fold == k)[0]; wd = winners(te, s); wd = wd[dec[wd]]
        med = np.median(s[wd]); out[te] = (s[te] - med) / (1.4826 * np.median(np.abs(s[wd] - med)) + 1e-12)
    return out
def report(name, s):
    w = winners(np.arange(len(s)), s); tg = w[~dec[w]]; o = tg[np.argsort(-s[tg], kind='stable')]
    fd = np.cumsum(entr[o]) * (1 + 1 / R) / np.arange(1, len(o) + 1)
    mt = [int(np.where(fd <= t)[0].max() + 1) if (fd <= t).any() else 0 for t in (0.005, 0.01, 0.015, 0.02)]
    # bulk-null weights from the same score
    wdv = w[dec[w]]; twv = w[~dec[w]]; thr = np.median(s[wdv]); bt, bd = twv[s[twv] < thr], wdv[s[wdv] < thr]
    r = np.ones(len(s)); pr = np.zeros(len(s))
    for k in np.unique(fold):
        tr = np.r_[bt[fold[bt] != k], bd[fold[bd] != k]]
        m = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.05, max_leaf_nodes=31, random_state=0).fit(X[tr], (~dec[tr]).astype(int))
        te = wdv[fold[wdv] == k]; pr[te] = m.predict_proba(X[te])[:, 1]
    r[wdv] = np.clip(pr[wdv] / np.maximum(1 - pr[wdv], 1e-9) * len(bd) / len(bt), 0.2, 5); r[wdv] /= r[bd].mean()
    qs = []
    for nm, wv in (('plain', None), ('bulk-null', r)):
        kept, q = tdc(w, s, wv)
        for qq in (0.01, 0.05):
            a = kept[(~dec[kept]) & (q <= qq)]; qs.append(f'{nm} q<={qq}: {len(a)} FDP_c {entr[a].sum() * (1 + 1 / R) / max(1, len(a)):.3f}')
    print(f'{name:16s} N at FDP_c 0.5/1/1.5/2%: {mt} | ' + ' | '.join(qs), flush=True)
print(ARM, len(f), FEAT, flush=True)
report('tool LDA', f.stage2_score.values.astype(np.float64))
ens = np.zeros(len(Z))
for i, seed in enumerate((0, 10, 20, 30, 40)):
    h = head(seed); ens += h
    if i == 0: report('MLP seed 0', h)
report('MLP ens5', ens / 5)
