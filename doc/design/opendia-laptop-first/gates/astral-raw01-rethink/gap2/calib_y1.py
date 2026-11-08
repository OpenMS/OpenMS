"""Mechanism check: is the decoy y1 deficit a scan-range cut-off? Zero-trace rate of top-6 fragment slots by product
m/z (1 Th bins, 100-200 Th) and group; plus per-group share of precursors with a top-6 slot below 150 Th, and that
share among tail winners. Design diagnostic."""
import os, sys, sqlite3, numpy as np, pandas as pd
sys.path.insert(0, '.')
RUN = os.environ.get('OPENDIA_RUN', '01')
import oswlike_score as S
Sd = S.S
meta = pd.read_parquet(Sd + '/meta.parquet'); lib = np.load(Sd + '/lib.npy'); frag = np.load(Sd + '/frag.npy', mmap_mode='r')
pid2row = pd.Series(np.arange(len(meta)), index=meta.precursor_id.values)
fmax = np.zeros(lib.shape, np.float32)
for b in range(0, len(meta), 50000):
    fmax[b:b + 50000] = np.asarray(frag[b:b + 50000]).max(2)
ent_p = S.m.groupby('precursor_id').ent.first().reindex(meta.precursor_id).fillna(0).astype(bool).values
grp = np.where(meta.decoy.values.astype(bool), 'decoy', np.where(ent_p, 'entrapment', 'human'))
c = sqlite3.connect('/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-entrapment-20261107/lib/library.pqp')
tr = pd.read_sql('select m.PRECURSOR_ID pid, t.LIBRARY_INTENSITY li, t.ID tid, t.PRODUCT_MZ fmz from TRANSITION t join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID', c)
tr = tr.sort_values(['pid', 'li', 'tid'], ascending=[True, False, True]); tr['slot'] = tr.groupby('pid').cumcount()
tr = tr[(tr.slot < 6) & tr.pid.isin(pid2row.index)].copy()
r = pid2row[tr.pid.values].values; s = tr.slot.values
tr['grp'] = grp[r]; tr = tr[lib[r, s] > 0]; r = pid2row[tr.pid.values].values; s = tr.slot.values
tr['zero'] = fmax[r, s] <= 0
lo = tr[tr.fmz < 200].copy(); lo['bin'] = (lo.fmz // 5 * 5).astype(int)
t = lo.pivot_table(index='bin', columns='grp', values='zero', aggfunc=['mean', 'size'], observed=True)
print(f'raw {RUN}: zero-trace rate (and slot count) by fragment m/z, 5 Th bins below 200')
print(t.round(3).to_string())
print('min nonzero fragment m/z seen:', round(float(tr[~tr.zero].fmz.min()), 3))
below = tr[tr.fmz < 150].groupby('pid').size()
print('share of precursors with a top-6 slot < 150 Th:', {g_: round(float(below.reindex(meta.precursor_id[grp == g_]).fillna(0).gt(0).mean()), 4) for g_ in ('human', 'entrapment', 'decoy')})
z = np.load(f'{Sd}/calib_scores.npz')
for name in ('lda', 'mlp'):
    sc = z[name]; w = S.winners(np.arange(len(sc)), sc); kept, q = S.tdc_q(w, sc)
    dw = w[S.dec[w]]; ew = w[(~S.dec[w]) & S.entr[w]]
    tq = sc[kept[(~S.dec[kept]) & (q <= 0.05)]].min()
    f = lambda idx: below.reindex(S.prec[idx]).fillna(0).gt(0)
    print(f'{name.upper()}: share with slot<150 Th - all decoy winners {f(dw).mean():.3f}, decoys above 5% cut {f(dw[sc[dw] >= tq]).mean():.3f}; all entrapment winners {f(ew).mean():.3f}, above cut {f(ew[sc[ew] >= tq]).mean():.3f}')
    for lab, m_ in (('decoys WITHOUT slot<150', ~f(dw).values), ('decoys WITH slot<150', f(dw).values)):
        print(f'   tail survival {lab}: {(sc[dw[m_]] >= tq).mean() * 1e4:.2f}e-4 (n {m_.sum()})')
    print(f'   tail survival entrapment: {(sc[ew] >= tq).mean() * 1e4:.2f}e-4')
