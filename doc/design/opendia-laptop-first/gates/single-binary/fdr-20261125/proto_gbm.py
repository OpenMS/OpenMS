#!/usr/bin/env python3
"""Python stand-in for gbmRescore (sklearn HGB, same labels/negatives/weights/folds/validation split/rounds), applied to
the exported network scores and features, for screening variants. Not bit-identical to C++. Diagnosis only.

Variants (argv[1], comma list):
  base         C++ settings (hard 2 %, pos q 0.01, 3 rounds, pos weight n_neg/n_pos, min_samples_leaf 20)
  leaf200      min_samples_leaf 200
  hess         min_samples_leaf 20 but positives weighted 1 (class balance through no weight)
  stack        base + network score as an extra feature
  shrink       final score = network + GBM margin standardised (half weight)
Output per run: N at own q<=0.01, E, combined bound, N@FDP_c1%."""
import sqlite3, sys, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
LIB = '/home/singjc/.jcode/scratch/perf/lib/search.pqp'
c = sqlite3.connect(f'file:{LIB}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID precursor_id, p.TRAML_ID t, p.DECOY d,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
pq['key'] = np.where(pq.d == 1, pq.t.str.replace('^DECOY_', '', regex=True), pq.t)
P = pq.set_index('precursor_id')
META = ('candidate_row', 'precursor_id', 'is_decoy', 'fold', 'rt_apex', 'score', 'q_value', 'network_score', 'network_q_value')

def tdc(prec, pair, dec, fold, s):
    """pooledPairedTDC: winners per precursor, per-fold decoy median/MAD standardisation, better of pair, q."""
    df = pd.DataFrame({'prec': prec, 'pair': pair, 'dec': dec, 'fold': fold, 's': s, 'row': np.arange(len(s))})
    w = df.sort_values(['s', 'row'], ascending=[False, True]).drop_duplicates('prec')
    st = np.empty(len(w))
    for f_ in np.unique(w.fold):
        m = w.fold.values == f_
        d = w.s.values[m & (w.dec.values == 1)]
        med = np.median(d); mad = 1.4826 * np.median(np.abs(d - med))
        if not mad > 0:
            print(f'   [warn] fold {f_}: decoy MAD 0 (n {len(d)}, unique {len(np.unique(d))}); using SD', file=sys.stderr)
            mad = np.std(d) or 1.0
        st[m] = (w.s.values[m] - med) / mad
    w = w.assign(st=st).sort_values(['st', 'dec'], ascending=[False, False])
    pb = w.drop_duplicates('pair')
    dd = np.cumsum(pb.dec.values); tt = np.cumsum(1 - pb.dec.values)
    qv = np.minimum(1, (dd + 1) / np.maximum(tt, 1))
    qv = np.minimum.accumulate(qv[::-1])[::-1]
    return w.row.values, pb.row.values, qv

def rescore(f, X, variant, seed):
    prec = f.precursor_id.values; pair = f.key.values; dec = f.is_decoy.values; fold = f.fold.values
    ph = f.pair_hash.values
    s = f.network_score_raw.values.copy()
    leaf = 200 if variant == 'leaf200' else 20
    for rnd in range(3):
        winners, ranked, qv = tdc(prec, pair, dec, fold, s)
        pos = np.sort(ranked[(qv <= 0.01) & (dec[ranked] == 0)])
        dw = winners[dec[winners] == 1]
        dw = dw[np.argsort(-s[dw], kind='stable')]
        hard = set(prec[dw[:max(1, int(np.ceil(0.02 * len(dw))))]])
        neg = np.where((dec == 1) & np.isin(prec, list(hard)))[0]
        nxt = np.zeros(len(s))
        for k in range(3):
            fp, fn = pos[fold[pos] != k], neg[fold[neg] != k]
            rows = np.concatenate([fp, fn]); y = np.r_[np.ones(len(fp)), np.zeros(len(fn))]
            wpos = 1.0 if variant == 'hess' else len(fn) / max(1, len(fp))
            w = np.where(y == 1, wpos, 1.0)
            m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=31, min_samples_leaf=leaf,
                                               early_stopping=True, validation_fraction=0.1, n_iter_no_change=10, random_state=seed + 1000 * (rnd + 1) + k)
            m.fit(X[rows], y, sample_weight=w)
            te = fold == k
            nxt[te] = m.decision_function(X[te])
        s = nxt
    return s

def evaluate(f, s):
    """Keep best window per precursor, competition, q; report entrapment metrics."""
    winners, ranked, qv = tdc(f.precursor_id.values, f.key.values, f.is_decoy.values, f.fold.values, s)
    r = pd.DataFrame({'row': ranked, 'q': qv})
    r['dec'] = f.is_decoy.values[ranked]; r['e'] = f.e.values[ranked]
    t = r[r.dec == 0]
    acc = t[t.q <= 0.01]
    fdp = np.cumsum(t.e.values) * 2 / np.arange(1, len(t) + 1)
    ok = np.where(fdp <= 0.01)[0]
    return len(acc), int(acc.e.sum()), 2 * acc.e.sum() / max(1, len(acc)), int(ok.max() + 1) if len(ok) else 0

variants = sys.argv[1].split(',')
for path in sys.argv[2:]:
    f = pd.read_csv(path, sep='\t')
    f['key'] = f.precursor_id.map(P.key); f['e'] = f.precursor_id.map(P.e).fillna(0).astype(int)
    f['pair_hash'] = pd.util.hash_pandas_object(f.key, index=False).values % 100
    feats = [x for x in f.columns if x not in META + ('key', 'e', 'pair_hash')]
    # Network raw score is not exported; the calibrated kept score is monotone per fold for kept windows only.
    # Use the tool's own GBM input instead: network score on kept windows, other windows get their calibrated value too
    # (NOT_KEPT rows set below every kept row) -- matches the C++ input `scored.score`.
    f['network_score_raw'] = np.where(f.network_score > -1e29, f.network_score, -1e3)
    X = f[feats].values.astype(np.float32)
    seed = int(path.split('_s')[1].split('/')[0])
    line = [path.split('/')[0]]
    line.append('net ' + '/'.join(map(lambda v: f'{v:.4f}' if isinstance(v, float) else str(v), evaluate(f, f.network_score_raw.values))))
    line.append('tool ' + '/'.join(map(lambda v: f'{v:.4f}' if isinstance(v, float) else str(v), evaluate(f, np.where(f.score > -1e29, f.score, -1e3)))))
    for v in variants:
        Xv = np.c_[X, f.network_score_raw.values] if v == 'stack' else X
        s = rescore(f, Xv, 'base' if v == 'stack' else v, seed)
        line.append(f'{v} ' + '/'.join(map(lambda x: f'{x:.4f}' if isinstance(x, float) else str(x), evaluate(f, s))))
    print('  '.join(line), flush=True)
