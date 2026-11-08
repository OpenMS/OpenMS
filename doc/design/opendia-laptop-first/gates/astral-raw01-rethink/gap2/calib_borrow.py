"""Do the entrapment hits that break decoy calibration borrow signal from co-eluting human precursors? Design diagnostic.

For a precursor P, its neighbour overlap = the largest number of P's top-6 library fragments (within 20 ppm) shared with
any one HUMAN target precursor whose precursor m/z lies within +-1 Th of P's and whose calibrated library RT lies within
+-30 s. Compared for: entrapment winners above the q<=0.01 / 0.05 score cuts; decoy winners above the same cuts; and
random samples of all entrapment and decoy precursors. If the hits are enriched for overlap relative to decoys, false
targets are scoring on borrowed real signal that decoys cannot reach (pseudo-reversed, mass-shifted), which is a
decoy-model mismatch rather than learner overfitting.
"""
import os, sys, sqlite3, numpy as np, pandas as pd
sys.path.insert(0, '.')
RUN = os.environ.get('OPENDIA_RUN', '01')
import oswlike_score as S
LIB = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-entrapment-20261107/lib/library.pqp'
cache = f'{S.S}/calib_scores.npz'
if os.path.exists(cache):
    z = np.load(cache); lda, mlp = z['lda'], z['mlp']
else:
    lda = S.crossfit(S.ALL); mlp = S.crossfit(S.ALL, kind='mlp', seed=0); np.savez(cache, lda=lda, mlp=mlp)

c = sqlite3.connect(LIB)
pr = pd.read_sql('select ID pid, PRECURSOR_MZ mz, DECOY dec from PRECURSOR', c).set_index('pid')
tr = pd.read_sql('select m.PRECURSOR_ID pid, t.PRODUCT_MZ fmz, t.LIBRARY_INTENSITY li from TRANSITION t join '
                 'TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID where t.DETECTING=1', c)
tr = tr.sort_values(['pid', 'li'], ascending=[True, False]).groupby('pid').head(6)
frags = tr.groupby('pid').fmz.apply(np.array).to_dict()
rt = S.m.groupby('precursor_id').cal_rt.first()
pr['rt'] = rt.reindex(pr.index).values
ent_p = S.m.groupby('precursor_id').ent.first().reindex(pr.index).fillna(0).astype(bool)
hum = pr[(pr.dec == 0) & (~ent_p) & pr.rt.notna()].sort_values('mz')
hmz = hum.mz.values; hrt = hum.rt.values; hid = hum.index.values

def overlap(pid):
    if pid not in frags or not np.isfinite(pr.rt.get(pid, np.nan)): return np.nan
    a, b = np.searchsorted(hmz, [pr.mz[pid] - 1, pr.mz[pid] + 1])
    f = frags[pid]; best = 0
    for k in range(a, b):
        if hid[k] == pid or abs(hrt[k] - pr.rt[pid]) > 30: continue
        g = frags.get(hid[k])
        if g is None: continue
        n = int((np.abs(f[:, None] - g[None, :]) <= 20e-6 * f[:, None]).any(1).sum())
        best = max(best, n)
    return best

def summ(label, pids):
    o = np.array([overlap(p) for p in pids], float); o = o[np.isfinite(o)]
    print(f'  {label:42s} n {len(o):5d}  mean {o.mean():.2f}  share>=2 {(o >= 2).mean():.3f}  share>=3 {(o >= 3).mean():.3f}', flush=True)

rng = np.random.default_rng(0)
allent = pr.index[(pr.dec == 0) & ent_p.values]; alldec = pr.index[pr.dec == 1]
print(f'raw {RUN}')
summ('random entrapment precursors', rng.choice(allent, 3000, replace=False))
summ('random decoy precursors', rng.choice(alldec, 3000, replace=False))
for name, sc in (('LDA', lda), ('MLP', mlp)):
    w = S.winners(np.arange(len(sc)), sc); kept, q = S.tdc_q(w, sc)
    for qq in (0.01, 0.05):
        t = sc[kept[(~S.dec[kept]) & (q <= qq)]].min()
        e = w[(~S.dec[w]) & S.entr[w] & (sc[w] >= t)]; d = w[S.dec[w] & (sc[w] >= t)]
        summ(f'{name} entrapment winners above q<={qq} cut', S.prec[e])
        summ(f'{name} decoy winners above q<={qq} cut', S.prec[d])
