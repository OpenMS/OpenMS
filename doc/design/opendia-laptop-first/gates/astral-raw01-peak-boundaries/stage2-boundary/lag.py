"""Why do we still lag? Raw 01 design diagnostic (target/decoy labels only for fitting; OSW used only to locate misses).

1. Funnel of same-library OSW+Percolator IDs through the search+boundary result, and which features still separate the
   missed right-peak targets from the decoys that compete with them.
2. Combiner test on identical inputs: cross-fitted semi-supervised linear LDA vs a small MLP, each on (a) the search+boundary
   features the tool used and (b) the union with the seven extraction features. Fold k is scored only by models fitted on
   the other folds (PyProphet-style: winners, decoy winners as negatives, target winners at q<=0.15 then 0.05 as positives).
   If the MLP accepts many more at the same decoy-count estimate, the combiner is the limit; if not, the information is.
"""
import sys, time, sqlite3
import numpy as np, pandas as pd, torch, torch.nn as nn

sys.path.insert(0, '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-design-raw01-20261107/cnn-topaz-design')
from evaluate import paired_pooled

torch.manual_seed(0)
torch.set_num_threads(5)
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
D = B + '/astral-design-raw01-20261107'
LIB = B + '/astral-entrapment-20261107/lib/library.pqp'
OSW = '/home/singjc/Documents/datasets/proteobench-astral-entrapment/100K_sampled_fasta/opendia_results/OpenDIA.results.tsv'

f = pd.read_csv(D + '/stage2-boundary-design/boundary/features.tsv', sep='\t')
x = pd.read_csv(D + '/stage2-design/features.tsv', sep='\t',
                usecols=['sample_id', 'fragment_fraction', 'fragment_coelution', 'library_correlation', 'ms2_abs_ppm', 'ms1_abs_ppm', 'ms1_coelution', 'ms1_isotope'])
cr = pd.read_csv(D + '/stage2-boundary-design/boundary/crossfit.tsv', sep='\t', usecols=['sample_id', 'in_rt_window', 'precursor_best', 'q_value'])
s = pd.read_csv(D + '/pipeline/search.features.tsv', sep='\t', usecols=['candidate_row', 'rt_apex', 'apex_intensity']).rename(columns={'candidate_row': 'sample_id'})
df = f.merge(x, on='sample_id').merge(cr, on='sample_id').merge(s, on='sample_id')
c = sqlite3.connect('file:' + LIB + '?mode=ro', uri=True)
p = pd.read_sql('select ID precursor_id, TRAML_ID, DECOY from PRECURSOR', c)
p['pair_key'] = np.where(p.DECOY == 1, p.TRAML_ID.str.replace('^DECOY_', '', regex=True), p.TRAML_ID)
df = df.merge(p[['precursor_id', 'pair_key']], on='precursor_id')
df['q'] = pd.to_numeric(df.q_value, errors='coerce')
SEARCH = ['cnn_z', 'transition_support', 'log_matched_scans', 'log_peak_width', 'abs_library_rt_error']
BOUND = [c for c in f.columns if c.startswith('boundary_')]
EXTR = ['fragment_fraction', 'fragment_coelution', 'library_correlation', 'ms2_abs_ppm', 'ms1_abs_ppm', 'ms1_coelution', 'ms1_isotope']
print(f'windows in RT window: {len(df)}')

# ---- 1. Funnel and separability of what is still missed.
lib = pd.read_parquet(B + '/astral-diann-gap-20261108/lib_targets.parquet')
lib['mkey'] = lib.modseq + '/' + lib.charge.astype(str)
o = pd.read_csv(OSW, sep='\t', usecols=['run_name', 'FullPeptideName', 'Charge', 'RT', 'd_score'])
o = o[o.run_name.str.endswith('_01')]
o['mkey'] = o.FullPeptideName + '/' + o.Charge.astype(str)
o = o.merge(lib[['mkey', 'precursor_id']], on='mkey').rename(columns={'RT': 'ref_rt'})
win = df[df.precursor_best == 1]
m = win.merge(o[['precursor_id', 'ref_rt', 'd_score']], on='precursor_id')
m['right'] = (m.rt_apex - m.ref_rt).abs() <= 10
acc = m.q <= 0.01
print(f'\nOSW 1% IDs {len(o)} | our winner at OSW RT {int(m.right.sum())} | accepted {int(acc.sum())} | right peak not accepted {int((m.right & ~acc).sum())}')
miss = m[m.right & ~acc]
print('missed right-peak q quantiles', miss.q.quantile([.25, .5, .75, .9]).round(3).to_dict())
dec = win[win.is_decoy == 1]
cut = win[(win.is_decoy == 0) & (win.q <= 0.01)].stage2_score.min()
# Competitors: decoy winners in the same stage-2 score range as the missed targets.
lo, hi = miss.stage2_score.quantile(0.1), cut
comp = dec[(dec.stage2_score >= lo) & (dec.stage2_score < hi)]
tgt = miss[(miss.stage2_score >= lo)]


def auc(a, b):
    a, b = a.dropna().values, b.dropna().values
    r = pd.Series(np.r_[a, b]).rank().values
    return (r[:len(a)].sum() - len(a) * (len(a) + 1) / 2) / (len(a) * len(b))


print(f'\nwithin the score band where the misses sit (stage2 {lo:.1f}..{hi:.1f}): missed OSW targets {len(tgt)}, decoy winners {len(comp)}')
print('  per-feature AUC missed-target vs decoy (0.5 = no information left; <0.5 means decoys look better):')
for k in SEARCH + BOUND + EXTR + ['apex_intensity']:
    print(f'    {k:30s} {auc(tgt[k], comp[k]):.2f}   target med {tgt[k].median():9.3f}  decoy med {comp[k].median():9.3f}')

# ---- 2. Combiner test.
folds = sorted(df.fold.unique())
fold_of = df.fold.values
is_dec = df.is_decoy.values.astype(bool)
prec = df.precursor_id.values


def winners_idx(idx, score):
    sub = pd.DataFrame({'i': idx, 'p': prec[idx], 's': score[idx]}).sort_values(['s', 'i'], ascending=[False, True])
    return sub.drop_duplicates('p').i.values


def select(idx, score, fdr):
    w = winners_idx(idx, score)
    w = w[np.argsort(-score[w], kind='stable')]
    d = is_dec[w]
    T, Dc = np.cumsum(~d), np.cumsum(d)
    q = np.minimum.accumulate(np.minimum(1, (Dc + 1) / np.maximum(T, 1))[::-1])[::-1]
    return w[(~d) & (q <= fdr)], w[d]


def lda(X, pos, neg, ridge=1e-3):
    mu, sd = X[np.r_[pos, neg]].mean(0), X[np.r_[pos, neg]].std(0) + 1e-9
    Z = (X - mu) / sd
    a, b = Z[pos], Z[neg]
    S = np.cov(a.T) * (len(a) - 1) + np.cov(b.T) * (len(b) - 1)
    S = S / (len(a) + len(b) - 2) + ridge * np.eye(X.shape[1])
    w = np.linalg.solve(S, a.mean(0) - b.mean(0))
    return Z @ w


def mlp_fit(X, pos, neg, epochs=6):
    mu, sd = X[np.r_[pos, neg]].mean(0), X[np.r_[pos, neg]].std(0) + 1e-9
    Xt = torch.tensor((X - mu) / sd, dtype=torch.float32)
    net = nn.Sequential(nn.Linear(X.shape[1], 32), nn.GELU(), nn.Linear(32, 16), nn.GELU(), nn.Linear(16, 1))
    opt = torch.optim.AdamW(net.parameters(), lr=3e-3, weight_decay=1e-3)
    idx = np.r_[pos, neg]
    y = torch.tensor(np.r_[np.ones(len(pos)), np.zeros(len(neg))], dtype=torch.float32)
    wpos = len(neg) / max(1, len(pos))
    loss_fn = nn.BCEWithLogitsLoss(pos_weight=torch.tensor(wpos))
    for _ in range(epochs):
        perm = np.random.permutation(len(idx))
        for i in range(0, len(idx), 2048):
            b = perm[i:i + 2048]
            out = net(Xt[idx[b]]).squeeze(1)
            loss = loss_fn(out, y[b])
            opt.zero_grad(); loss.backward(); opt.step()
    with torch.no_grad():
        return np.concatenate([net(Xt[i:i + 200000]).squeeze(1).numpy() for i in range(0, len(Xt), 200000)])


def crossfit(cols, kind):
    X = df[cols].to_numpy(np.float64)
    out = np.full(len(df), np.nan)
    for k in folds:
        tr = np.where(fold_of != k)[0]
        score = X[:, 0].copy()  # rank by cnn_z first
        for it, fdr in enumerate([0.15, 0.05, 0.05]):
            pos, neg = select(tr, score, fdr)
            score = lda(X, pos, neg)
        if kind == 'mlp':
            for it in range(2):
                pos, neg = select(tr, score, 0.05)
                score = mlp_fit(X, pos, neg)
        out[fold_of == k] = score[fold_of == k]
    return out


res = {}
for name, cols in [('search+boundary', SEARCH + BOUND), ('search+boundary+extraction', SEARCH + BOUND + EXTR)]:
    for kind in ['lda', 'mlp']:
        t0 = time.time()
        sc = crossfit(cols, kind)
        tmp = df[['sample_id', 'precursor_id', 'is_decoy', 'fold', 'in_rt_window', 'pair_key']].assign(score=sc)
        best, w = paired_pooled(tmp, 'score')
        t = best[best.is_decoy == 0]
        counts = [int((t.q <= q).sum()) for q in (0.005, 0.01, 0.02, 0.05)]
        # Null symmetry: targets/decoys below the decoy-winner median, relative to totals.
        wz = w.assign(zz=w.z)
        nt, nd = (wz.is_decoy == 0).sum(), (wz.is_decoy == 1).sum()
        low = wz[wz.zz < wz[wz.is_decoy == 1].zz.median()]
        R = ((low.is_decoy == 0).sum() / (low.is_decoy == 1).sum()) / (nt / nd)
        osw_rec = t[t.q <= 0.01].precursor_id.isin(o.precursor_id).sum()
        res[(name, kind)] = counts
        print(f'{name:28s} {kind}: q<=0.005/0.01/0.02/0.05 = {counts} | OSW IDs among 1% {osw_rec} | R(decoy median) {R:.3f} | {time.time() - t0:.0f} s')
