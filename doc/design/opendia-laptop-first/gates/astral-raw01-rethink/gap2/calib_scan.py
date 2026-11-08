"""How much of the decoy-q excess does the out-of-scan-range artefact explain? Design diagnostic.
Fragments < 150 Th are never acquired (the lowest non-zero m/z is 150.06), but the library keeps them in the top 6, and decoys more
often (9.7% of decoys, 5.7% of entrapment, 5.2% of human precursors), including decoy transitions at 55-100 Th, below
the 100 Th assay limit. Simulated fix: drop every target-decoy pair in which either member has a top-6 slot < 150 Th
(a symmetric library filter), then redo TDC q, FDP_c and the tail ratio on the cached scores. Not a full fix (it
removes pairs instead of re-picking fragments), but it bounds the effect."""
import os, sys, sqlite3, numpy as np, pandas as pd
sys.path.insert(0, '.')
RUN = os.environ.get('OPENDIA_RUN', '01')
import oswlike_score as S
c = sqlite3.connect('/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-entrapment-20261107/lib/library.pqp')
tr = pd.read_sql('select m.PRECURSOR_ID pid, t.LIBRARY_INTENSITY li, t.ID tid, t.PRODUCT_MZ fmz from TRANSITION t join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID', c)
tr = tr.sort_values(['pid', 'li', 'tid'], ascending=[True, False, True]); tr['slot'] = tr.groupby('pid').cumcount()
bad_p = np.array(sorted(set(tr[(tr.slot < 6) & (tr.fmz < 150)].pid)))
bad_pair = pd.unique(S.pair[np.isin(S.prec, bad_p)])
keep = ~pd.Series(S.pair).isin(set(bad_pair)).values
print(f'raw {RUN}: pairs dropped {len(bad_pair)}; candidates kept {keep.mean():.3f}')
z = np.load(f'{S.S}/calib_scores.npz'); ens = np.load(f'{S.S}/calib_ens.npy')
for name, sc in (('LDA', z['lda']), ('MLP seed 0', z['mlp']), ('MLP 5-seed ensemble', ens)):
    for lab, idx in (('all pairs', np.arange(len(sc))), ('in-range pairs only', np.where(keep)[0])):
        w = S.winners(idx, sc); kept, q = S.tdc_q(w, sc)
        a = kept[(~S.dec[kept]) & (q <= 0.01)]; b = kept[(~S.dec[kept]) & (q <= 0.05)]
        f = lambda x: S.entr[x].sum() * (1 + 1 / S.R_ENT) / max(1, len(x))
        dw = w[S.dec[w]]; ew = w[(~S.dec[w]) & S.entr[w]]
        rr = []
        for qq in (0.01, 0.05):
            t = sc[kept[(~S.dec[kept]) & (q <= qq)]].min(); rr.append(f'{((sc[ew] >= t).mean()) / max(1e-12, (sc[dw] >= t).mean()):.2f}')
        tg = kept[~S.dec[kept]]; o = tg[np.argsort(-sc[tg], kind='stable')]; fd = np.cumsum(S.entr[o]) * (1 + 1 / S.R_ENT) / np.arange(1, len(o) + 1)
        n1 = int(np.where(fd <= 0.01)[0].max() + 1)
        print(f'  {name:20s} {lab:20s} q<=0.01 {len(a):5d} (FDP_c {f(a):.3f}, lower {S.entr[a].mean():.3f}) | q<=0.05 {len(b):5d} (FDP_c {f(b):.3f}) | N at FDP_c 1% {n1} | S_e/S_d at 1%/5% cuts {"/".join(rr)}')
