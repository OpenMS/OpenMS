"""Is the entrapment, not the decoy, the odd one out? Design diagnostic.

1. Per fragment: zero-trace rate by fragment m/z bin, for human / entrapment / decoy. Both decoys and entrapment
   peptides are false, so for valid entrapment they should match each other (human targets mix true and false).
2. Stratified tail ratio: re-weight entrapment winners so that their distribution over structural strata
   (number of empty fragment traces x number of peak candidates) matches the decoy winners. Then recompute S_e/S_d
   at the q<=0.01/0.02/0.05 cuts. If the excess shrinks towards 1, part of it is an entrapment construction artefact,
   not decoy failure.
"""
import os, sys, sqlite3, numpy as np, pandas as pd
sys.path.insert(0, '.')
RUN = os.environ.get('OPENDIA_RUN', '01')
import oswlike_score as S
Sd = S.S
meta = pd.read_parquet(Sd + '/meta.parquet'); lib = np.load(Sd + '/lib.npy'); frag = np.load(Sd + '/frag.npy', mmap_mode='r')
pid2row = pd.Series(np.arange(len(meta)), index=meta.precursor_id.values)
nz = np.zeros(len(meta), int)
fmax = np.zeros(lib.shape, np.float32)
for b in range(0, len(meta), 50000):
    f = np.asarray(frag[b:b + 50000]); fmax[b:b + 50000] = f.max(2)
    nz[b:b + 50000] = ((fmax[b:b + 50000] <= 0) & (lib[b:b + 50000] > 0)).sum(1)
ent_p = S.m.groupby('precursor_id').ent.first().reindex(meta.precursor_id).fillna(0).astype(bool).values
dec_p = meta.decoy.values.astype(bool)
grp = np.where(dec_p, 'decoy', np.where(ent_p, 'entrapment', 'human'))

# 1. per-fragment zero rate by fragment m/z (slot m/z from the PQP, top-6 by library intensity as in build_xic)
c = sqlite3.connect('/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-entrapment-20261107/lib/library.pqp')
tr = pd.read_sql('select m.PRECURSOR_ID pid, t.LIBRARY_INTENSITY li, t.ID tid, t.PRODUCT_MZ fmz from TRANSITION t join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID', c)
tr = tr.sort_values(['pid', 'li', 'tid'], ascending=[True, False, True]); tr['slot'] = tr.groupby('pid').cumcount()
tr = tr[(tr.slot < 6) & tr.pid.isin(pid2row.index)]
r = pid2row[tr.pid.values].values; s = tr.slot.values
tr = tr.assign(grp=grp[r], zero=(fmax[r, s] <= 0) & (lib[r, s] > 0), mapped=lib[r, s] > 0)
tr = tr[tr.mapped]
tr['fbin'] = pd.cut(tr.fmz, [100, 300, 500, 700, 900, 1200, 1800])
print(f'raw {RUN}: per-fragment zero-trace rate by fragment m/z')
print(tr.pivot_table(index='fbin', columns='grp', values='zero', observed=True).round(3).to_string())
print('fragment count share by bin:'); print(pd.crosstab(tr.fbin, tr.grp, normalize='columns').round(3).to_string())

# 2. stratified tail ratio
z = np.load(f'{Sd}/calib_scores.npz')
ncand = S.m.groupby('precursor_id').size()
for name in ('lda', 'mlp'):
    sc = z[name]; w = S.winners(np.arange(len(sc)), sc); kept, q = S.tdc_q(w, sc)
    dw = w[S.dec[w]]; ew = w[(~S.dec[w]) & S.entr[w]]
    def strata(idx):
        p = S.prec[idx]; return np.minimum(nz[pid2row[p].values], 3) * 10 + np.minimum(ncand.reindex(p).values, 5)
    sd, se = strata(dw), strata(ew)
    pd_ = pd.Series(sd).value_counts(normalize=True); pe_ = pd.Series(se).value_counts(normalize=True)
    wt = (pd_.reindex(se).fillna(0).values / pe_.reindex(se).values)
    out = []
    for qq in (0.01, 0.02, 0.05):
        t = sc[kept[(~S.dec[kept]) & (q <= qq)]].min()
        sd_ = (sc[dw] >= t).mean(); raw = (sc[ew] >= t).mean() / sd_; adj = (wt * (sc[ew] >= t)).sum() / wt.sum() / sd_
        out.append(f'q<={qq}: raw {raw:.2f} -> reweighted {adj:.2f}')
    print(f'{name.upper()} S_e/S_d ' + ' | '.join(out))
