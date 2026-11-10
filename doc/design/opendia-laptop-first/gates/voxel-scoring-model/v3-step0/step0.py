"""Voxel v3 step 0 (PLAN_STEP0.md): check A (brightness dependence) and check B (voxel score as an extra base LDA feature).
Usage: step0.py. Writes to stdout. Entrapment labels are read only for the check B metrics."""
import sqlite3, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
O = B + '/astral-voxel-pilot-20261110/v2'
R = 148983 / 151239
BASEF = {'raw01': B + '/astral-design-raw01-20261107/stage2-byion-design/byion/features.tsv',
         'raw02': B + '/astral-raw02-byion-design-20261109/byion/features.tsv'}
FEATS = ['cnn_z', 'transition_support', 'log_matched_scans', 'log_peak_width', 'abs_library_rt_error', 'boundary_log_fwhm', 'boundary_log_width',
         'boundary_coelution', 'boundary_library_correlation', 'boundary_apex_agreement', 'byion_b', 'byion_y']
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c).set_index('precursor_id').e

def auc(p, n): return mannwhitneyu(p, n).statistic / (len(p) * len(n))

# ---------- shared metric code (as evaluate_v2.py) ----------
def tdc(pop, s):
    w = pop.assign(s=s).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').copy()
    for f_, g in w[w.is_decoy == 1].groupby('fold'):
        med = g.s.median(); mad = 1.4826 * (g.s - med).abs().median()
        w.loc[w.fold == f_, 's'] = (w.loc[w.fold == f_, 's'] - med) / mad
    best = w.sort_values(['s', 'is_decoy'], ascending=[False, False]).drop_duplicates('pair_key').sort_values('s', ascending=False)
    T = np.cumsum(1 - best.is_decoy.values); D = np.cumsum(best.is_decoy.values)
    best['q'] = np.minimum.accumulate(np.minimum(1, (D + 1) / np.maximum(T, 1))[::-1])[::-1]
    return w, best

def n_at(s, e, t):
    o = np.argsort(-s, kind='stable'); e = e[o]
    f = np.cumsum(e) * (1 + 1 / R) / np.arange(1, len(e) + 1); ok = np.where(f <= t)[0]
    return int(ok.max() + 1) if len(ok) else 0

def metrics(pop, s):
    w, best = tdc(pop, s)
    tw = w[w.is_decoy == 0].assign(e=lambda d: d.precursor_id.map(ent).values); dw = w[w.is_decoy == 1]
    curve = [n_at(tw.s.values, tw.e.values, t) for t in (0.005, 0.01, 0.015, 0.02)]
    tb = best[best.is_decoy == 0].assign(e=lambda d: d.precursor_id.map(ent).values)
    qn = {q: (int((tb.q <= q).sum()), tb[tb.q <= q].e.sum() * (1 + 1 / R) / max((tb.q <= q).sum(), 1)) for q in (0.01, 0.05)}
    a = auc(tw[tw.e == 1].s.values, dw.s.values)
    return curve, qn, a, tw.set_index('precursor_id')[['s', 'e']]

# ---------- Python copy of crossfitStage2LDA ----------
def lda_crossfit(X, pid, dec, fold, iters=3, f0=0.15, fi=0.05, ridge=1e-3):
    out = np.zeros(len(X))
    for k in range(3):
        tr = np.where(fold != k)[0]; te = np.where(fold == k)[0]
        mu = X[tr].mean(0); sd = X[tr].std(0); sd[sd == 0] = 1
        Z = (X - mu) / sd
        score = X[:, 0].copy()
        for it in range(iters):
            d = pd.DataFrame({'i': tr, 's': score[tr], 'p': pid[tr]})
            win = d.sort_values(['s', 'i'], ascending=[False, True]).drop_duplicates('p').i.values
            o = win[np.lexsort((win, -score[win]))]
            isd = dec[o]; T = np.cumsum(~isd); D = np.cumsum(isd)
            q = np.minimum(1, (D + 1) / np.maximum(T, 1)); q = np.where(T == 0, 1, q)
            q = np.minimum.accumulate(q[::-1])[::-1]
            pos = o[(~isd) & (q <= (f0 if it == 0 else fi))]; neg = o[isd]
            mp, mn = Z[pos].mean(0), Z[neg].mean(0)
            S = ((Z[pos] - mp).T @ (Z[pos] - mp) + (Z[neg] - mn).T @ (Z[neg] - mn)) / (len(pos) + len(neg) - 2) + ridge * np.eye(X.shape[1])
            w = np.linalg.solve(S, mp - mn)
            score = Z @ w
        out[te] = Z[te] @ w
        print(f'      fold {k}: positives {len(pos)} negatives {len(neg)} weights ' + ' '.join(f'{v:+.2f}' for v in w))
    return out

rng = np.random.default_rng(5)
for raw in ('raw01', 'raw02'):
    pop = pd.read_csv(f'{O}/{raw}/population.tsv', sep='\t')
    bf = pd.read_csv(BASEF[raw], sep='\t', usecols=['sample_id', 'fold'] + FEATS + ['stage2_score']).set_index('sample_id').loc[pop.sample_id]
    V0 = np.load(f'{O}/{raw}/scores_V0.npy'); FS = np.load(f'{O}/{raw}/scores_FS.npy')
    xic = np.load(f'{O}/{raw}/xic.npy', mmap_mode='r')
    bright = np.empty(len(pop), np.float32)
    for b in range(0, len(pop), 200000): bright[b:b + 200000] = np.log1p(np.asarray(xic[b:b + 200000]).sum((1, 2)))
    print(f'\n===================== {raw}')
    # ---------- check A ----------
    tw = pd.read_csv(f'{O}/{raw}/train_windows.tsv', sep='\t')
    row = pd.Series(np.arange(len(pop)), index=pop.sample_id)
    anc = row.loc[tw[tw.group == 'anchor'].sample_id].values
    decw = np.where(pop.is_decoy.values == 1)[0]
    Rn = rng.choice(decw, 20 * len(anc), replace=False)
    edges = np.quantile(bright[anc], np.linspace(0, 1, 21)); edges[-1] += 1e-6
    ab = np.clip(np.searchsorted(edges, bright[anc], side='right') - 1, 0, 19); db = np.searchsorted(edges, bright[decw], side='right') - 1
    Mn = np.concatenate([rng.choice(decw[db == k], 20 * (ab == k).sum(), replace=(db == k).sum() < 20 * (ab == k).sum()) for k in range(20)])
    base_hits = None
    sc = {'BASE': bf.stage2_score.values, 'cnn_z': bf.cnn_z.values, 'V0': V0, 'FS': FS, 'brightness': bright}
    print(f'  check A: {len(anc)} anchors, R {len(Rn)}, M {len(Mn)} decoy windows')
    drops = {}
    for k, s in sc.items():
        a_r, a_m = auc(s[anc], s[Rn]), auc(s[anc], s[Mn]); drops[k] = a_r - a_m
        print(f'    {k:10s} AUC anchor>R {a_r:.4f}  anchor>M {a_m:.4f}  drop {a_r - a_m:+.4f}')
    # weaker view: non-anchor base hits (q<=0.01 base) vs M-matched decoys
    _, best = tdc(pop, bf.stage2_score.values)
    hits = best[(best.is_decoy == 0) & (best.q <= 0.01)]
    w_b, _ = tdc(pop, bf.stage2_score.values)
    hit_rows = pop.reset_index().merge(w_b[w_b.precursor_id.isin(hits.precursor_id)][['sample_id']], on='sample_id')['index'].values
    hit_rows = hit_rows[~np.isin(hit_rows, anc)]
    hb = np.clip(np.searchsorted(edges, bright[hit_rows], side='right') - 1, 0, 19)
    Mh = np.concatenate([rng.choice(decw[db == k], 20 * (hb == k).sum(), replace=True) for k in range(20) if (hb == k).sum()])
    print(f'    weaker view: {len(hit_rows)} non-anchor base hits vs {len(Mh)} matched decoys: ' +
          ' '.join(f'{k} {auc(s[hit_rows], s[Mh]):.4f}' for k, s in sc.items()))
    for k in ('V0', 'FS'):
        print(f'    shortcut rule {k}: drop {drops[k]:+.4f} >= 0.01 and >= 2x BASE drop ({2 * drops["BASE"]:+.4f}): {drops[k] >= 0.01 and drops[k] >= 2 * drops["BASE"]}')
    # ---------- check B ----------
    X0 = bf[FEATS].values.astype(float); pid = pop.precursor_id.values; dec = pop.is_decoy.values.astype(bool); fold = bf.fold.values
    res = {}
    for name, X in (('B0', X0), ('B1', np.c_[X0, V0]), ('B2', np.c_[X0, FS])):
        print(f'    {name}:')
        s = lda_crossfit(X, pid, dec, fold)
        # metrics use the base folds for per-fold decoy standardisation (as the base LDA output)
        p2 = pop.assign(fold=fold)
        curve, qn, a, tab = metrics(p2, s); res[name] = (tab, a)
        line = f'    {name} N@FDP_c 0.5/1/1.5/2%: {curve} | q<=0.01 {qn[0.01][0]} FDP_c {qn[0.01][1]:.3f} | q<=0.05 {qn[0.05][0]} FDP_c {qn[0.05][1]:.3f} | AUC {a:.4f}'
        if name != 'B0':
            m = tab.join(res['B0'][0].s.rename('b'), how='inner'); a_, b_, e_ = m.s.values, m.b.values, m.e.values.astype(int); d1, d2 = [], []
            for _ in range(1000):
                i = rng.integers(0, len(m), len(m))
                d1.append(n_at(a_[i], e_[i], 0.01) - n_at(b_[i], e_[i], 0.01)); d2.append(n_at(a_[i], e_[i], 0.02) - n_at(b_[i], e_[i], 0.02))
            line += f' | vs B0 90%: @1% {np.percentile(d1, 5):+.0f}..{np.percentile(d1, 95):+.0f} @2% {np.percentile(d2, 5):+.0f}..{np.percentile(d2, 95):+.0f}'
        print(line, flush=True)
