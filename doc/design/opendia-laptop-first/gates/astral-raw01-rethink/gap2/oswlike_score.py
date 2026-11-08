"""Teacher-free PyProphet-style scoring of the OSW-like candidates (raw 01). Labels: target/decoy only.

Learner: cross-fitted (pair-hash 3 folds) semi-supervised LDA. Iteration 0 ranks by one main subscore; then, for
several rounds, positives = target top-ranked candidates at q<=0.15 then 0.02 (PyProphet ss_initial_fdr / ss_iteration_fdr
defaults), negatives = decoy top-ranked candidates; refit LDA on all subscores. Final score = out-of-fold LDA. Peak
group per precursor = rank 1 by that score. Paired TDC (target vs its paired decoy) for q; also reported: N at
entrapment FDP_c <= 0.5/1/1.5/2% (combined method), and the confirmed-right count, to compare against OSW+Percolator
(1015 / 1151 / 1251 / 1286 at those FDP_c) and our tool (641 / 757 / 923 / 950).
Ablations drop one subscore group at a time. Also: an MLP variant (same SS loop, small MLP in the last two rounds) and
a 'no rank-1' variant (pool all candidates, best per precursor only at the end) to test the within-precursor prior.
"""
import os, sys, numpy as np, pandas as pd, torch, torch.nn as nn
S = os.environ['JCODE_SCRATCH_DIR'] + '/xic' + os.environ.get('OPENDIA_RUN', '01')
X0 = np.load(S + '/oswlike_X.npy'); m = pd.read_parquet(S + '/oswlike_meta.parquet')
NAMES = ['xcorr_coel', 'xcorr_shape', 'xcorr_coel_w', 'xcorr_shape_w', 'lib_corr', 'lib_dotprod', 'lib_manhattan',
         'lib_rmsd', 'lib_sangle', 'intensity', 'log_sn', 'n_sig', 'log_width', 'norm_rt', 'ms1_xcorr_coel',
         'ms1_xcorr_shape', 'ms1_iso_corr', 'ms1_log_sn', 'frag_min_share', 'frag_corr_min']
prec = m.precursor_id.values; dec = m.decoy.values.astype(bool); fold = m.fold.values; pair = m.pair.values
entr = m.ent.fillna(0).values.astype(bool); right = m.right.values
R_ENT = 148983 / 151239
torch.set_num_threads(12)

def winners(idx, score):
    o = idx[np.lexsort((idx, -score[idx]))]
    _, first = np.unique(prec[o], return_index=True)
    return o[first]

def tdc_q(w, score):
    """Paired TDC on precursor winners w: keep the better of target/decoy per pair; returns kept idx sorted and q."""
    d = pd.DataFrame({'i': w, 'p': pair[w], 's': score[w], 'd': dec[w]})
    d = d.sort_values(['s', 'd'], ascending=[False, False]).drop_duplicates('p').sort_values('s', ascending=False)
    dd = d.d.values; T = np.cumsum(~dd); Dn = np.cumsum(dd)
    q = np.minimum.accumulate(np.minimum(1, (Dn + 1) / np.maximum(T, 1))[::-1])[::-1]
    return d.i.values, q

def select(idx, score, fdr):
    w = winners(idx, score); kept, q = tdc_q(w, score)
    pos = kept[(~dec[kept]) & (q <= fdr)]
    neg = w[dec[w]]
    return pos, neg

def lda_fit(Z, pos, neg):
    a, b = Z[pos], Z[neg]
    Sw = (np.cov(a.T) * (len(a) - 1) + np.cov(b.T) * (len(b) - 1)) / (len(a) + len(b) - 2) + 1e-4 * np.eye(Z.shape[1])
    return np.linalg.solve(Sw, a.mean(0) - b.mean(0))

def mlp_fit(Z, pos, neg, seed):
    torch.manual_seed(seed)
    idx = np.r_[pos, neg]; y = torch.tensor(np.r_[np.ones(len(pos)), np.zeros(len(neg))], dtype=torch.float32)
    Zt = torch.tensor(Z, dtype=torch.float32)
    net = nn.Sequential(nn.Linear(Z.shape[1], 32), nn.GELU(), nn.Linear(32, 16), nn.GELU(), nn.Linear(16, 1))
    opt = torch.optim.AdamW(net.parameters(), lr=3e-3, weight_decay=1e-3)
    lf = nn.BCEWithLogitsLoss(pos_weight=torch.tensor(len(neg) / max(1, len(pos))))
    for _ in range(5):
        perm = np.random.default_rng(seed).permutation(len(idx))
        for b in range(0, len(idx), 4096):
            j = perm[b:b + 4096]; loss = lf(net(Zt[idx[j]]).squeeze(1), y[j]); opt.zero_grad(); loss.backward(); opt.step()
    with torch.no_grad():
        return np.concatenate([net(Zt[i:i + 500000]).squeeze(1).numpy() for i in range(0, len(Zt), 500000)])

def crossfit(cols, kind='lda', main='xcorr_shape_w', seed=0):
    X = X0[:, [NAMES.index(c) for c in cols]]
    Z = np.nan_to_num((X - np.nanmean(X, 0)) / (np.nanstd(X, 0) + 1e-9))
    out = np.zeros(len(Z))
    for k in range(3):
        tr = np.where(fold != k)[0]
        sc = Z[:, cols.index(main)].copy() if main in cols else Z[:, 0].copy()
        for it, fdr in enumerate([0.15, 0.02, 0.02, 0.02, 0.02]):
            pos, neg = select(tr, sc, fdr)
            sc = Z @ lda_fit(Z, pos, neg)
        if kind == 'mlp':
            for it in range(2):
                pos, neg = select(tr, sc, 0.02)
                sc = mlp_fit(Z, pos, neg, seed + it)
        # put this fold's score on a common scale: decoy-winner median / MAD of the held-out fold
        te = np.where(fold == k)[0]; wd = winners(te, sc); wd = wd[dec[wd]]
        med = np.median(sc[wd]); mad = 1.4826 * np.median(np.abs(sc[wd] - med)) + 1e-12
        out[te] = (sc[te] - med) / mad
    return out

def report(name, score):
    w = winners(np.arange(len(score)), score); kept, q = tdc_q(w, score)
    tg = kept[~dec[kept]]
    qs = []
    for qq in (0.01, 0.05):
        a = kept[(~dec[kept]) & (q <= qq)]
        fdp = entr[a].sum() * (1 + 1 / R_ENT) / max(1, len(a))
        qs.append(f'q<={qq}: {len(a)} (conf-right {int(right[a].sum())}, FDP_c {fdp:.3f})')
    # entrapment-matched counts over target winners ranked by score
    o = tg[np.argsort(-score[tg], kind='stable')]
    e = entr[o]; nn_ = np.arange(1, len(o) + 1); fdpc = np.cumsum(e) * (1 + 1 / R_ENT) / nn_
    mt = []
    for t in (0.005, 0.01, 0.015, 0.02):
        okk = np.where(fdpc <= t)[0]; mt.append(int(okk.max() + 1) if len(okk) else 0)
    print(f'{name:44s} {" | ".join(qs)} | N at FDP_c<=0.5/1/1.5/2%: {mt}', flush=True)

ALL = list(NAMES)
GROUPS = {'xcorr': ['xcorr_coel', 'xcorr_shape', 'xcorr_coel_w', 'xcorr_shape_w', 'frag_corr_min'],
          'library': ['lib_corr', 'lib_dotprod', 'lib_manhattan', 'lib_rmsd', 'lib_sangle'],
          'signal': ['intensity', 'log_sn', 'n_sig', 'frag_min_share'],
          'shape/rt': ['log_width', 'norm_rt'],
          'ms1': ['ms1_xcorr_coel', 'ms1_xcorr_shape', 'ms1_iso_corr', 'ms1_log_sn']}
if __name__ == '__main__':
    print('reference (raw 01, entrapment-matched N at FDP_c<=0.5/1/1.5/2%): OSW+Percolator [1015, 1151, 1251, 1286]; OpenDIA stage2 [641, 757, 923, 950]')
    report('OSW-like LDA, all 20 subscores', crossfit(ALL))
    for g, cols in GROUPS.items():
        report(f'  minus {g}', crossfit([c for c in ALL if c not in cols]))
    report('OSW-like MLP (SS loop), all 20', crossfit(ALL, 'mlp'))
    report('OSW-like MLP (SS loop), all 20, seed 1', crossfit(ALL, 'mlp', seed=1))
