"""Mechanism test: tryptic-background low-mass fragments. Design diagnostic.

Target and entrapment peptides end in K/R, so their y1/y2 (and some b2) ions match fragments that any co-isolated
tryptic peptide also produces, for example y1 = 147.113 (K) or 175.119 (R). Decoys here are pseudo-reversed with the
terminal K/R switched to another residue, so their low-ordinal y ions do not get this shared background. That gives
false targets free co-eluting signal that decoys cannot get, which is the anti-conservative direction.
Reports: decoy C-terminal residue distribution; zero-trace rate and median max intensity by (type, ordinal <= 3) for
each group; the share of top-6 slots that are y1/y2; and the tail ratio S_e/S_d after removing all y1/y2/b1/b2 slots
from the subscores, which is not possible with cached subscores, so the share of tail winners carrying such slots is shown instead.
"""
import os, sys, sqlite3, numpy as np, pandas as pd
sys.path.insert(0, '.')
RUN = os.environ.get('OPENDIA_RUN', '01')
import oswlike_score as S
Sd = S.S
meta = pd.read_parquet(Sd + '/meta.parquet'); lib = np.load(Sd + '/lib.npy'); frag = np.load(Sd + '/frag.npy', mmap_mode='r')
pid2row = pd.Series(np.arange(len(meta)), index=meta.precursor_id.values)
fmax = np.zeros(lib.shape, np.float32); fsum = np.zeros(lib.shape, np.float32)
for b in range(0, len(meta), 50000):
    f = np.asarray(frag[b:b + 50000]); fmax[b:b + 50000] = f.max(2); fsum[b:b + 50000] = f.sum(2)
ent_p = S.m.groupby('precursor_id').ent.first().reindex(meta.precursor_id).fillna(0).astype(bool).values
grp = np.where(meta.decoy.values.astype(bool), 'decoy', np.where(ent_p, 'entrapment', 'human'))
c = sqlite3.connect('/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-entrapment-20261107/lib/library.pqp')
pe = pd.read_sql('select m.PRECURSOR_ID pid, pe.UNMODIFIED_SEQUENCE seq from PRECURSOR_PEPTIDE_MAPPING m join PEPTIDE pe on pe.ID=m.PEPTIDE_ID', c).drop_duplicates('pid')
pe = pe[pe.pid.isin(pid2row.index)]; pe['grp'] = grp[pid2row[pe.pid].values]
print(f'raw {RUN}: C-terminal residue share (top 6)')
for g_, x in pe.groupby('grp'):
    print(f'  {g_:10s}', x.seq.str[-1].value_counts(normalize=True).head(6).round(3).to_dict())
tr = pd.read_sql('select m.PRECURSOR_ID pid, t.LIBRARY_INTENSITY li, t.ID tid, t.PRODUCT_MZ fmz, t.TYPE ty, t.ORDINAL ord, t.CHARGE fz from TRANSITION t join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID', c)
tr = tr.sort_values(['pid', 'li', 'tid'], ascending=[True, False, True]); tr['slot'] = tr.groupby('pid').cumcount()
tr = tr[(tr.slot < 6) & tr.pid.isin(pid2row.index)].copy()
r = pid2row[tr.pid.values].values; s = tr.slot.values
tr['grp'] = grp[r]; tr['mapped'] = lib[r, s] > 0; tr['zero'] = fmax[r, s] <= 0; tr['area'] = fsum[r, s]
tr = tr[tr.mapped]
tr['ion'] = np.where(tr.ord <= 3, tr.ty + tr.ord.astype(str) + '^' + tr.fz.astype(str), 'other')
print('share of top-6 slots by ion (ordinal<=3):')
print(pd.crosstab(tr.ion, tr.grp, normalize='columns').round(4).loc[lambda x: x.max(1) > 0.003].to_string())
print('zero-trace rate by ion:')
print(tr.pivot_table(index='ion', columns='grp', values='zero', observed=True).round(3).loc[lambda x: x.index.isin(pd.crosstab(tr.ion, tr.grp).loc[lambda y: y.min(1) > 2000].index)].to_string())
print('median log10 area (nonzero) by ion:')
nzt = tr[~tr.zero].assign(la=lambda x: np.log10(x.area + 1))
print(nzt.pivot_table(index='ion', columns='grp', values='la', aggfunc='median', observed=True).round(2).loc[lambda x: x.index.isin(pd.crosstab(tr.ion, tr.grp).loc[lambda y: y.min(1) > 2000].index)].to_string())
low = tr[tr.ion.str.match(r'^y[12]\^1$')].groupby('pid').size()
z = np.load(f'{Sd}/calib_scores.npz')
for name in ('lda', 'mlp'):
    sc = z[name]; w = S.winners(np.arange(len(sc)), sc); kept, q = S.tdc_q(w, sc)
    dw = w[S.dec[w]]; ew = w[(~S.dec[w]) & S.entr[w]]
    for qq in (0.05,):
        t = sc[kept[(~S.dec[kept]) & (q <= qq)]].min()
        e = S.prec[ew[sc[ew] >= t]]; d = S.prec[dw[sc[dw] >= t]]
        print(f'{name.upper()} above q<={qq} cut: share carrying a y1/y2 slot: entrapment {low.reindex(e).fillna(0).gt(0).mean():.3f} (n {len(e)}), decoys {low.reindex(d).fillna(0).gt(0).mean():.3f} (n {len(d)}); all entrapment {low.reindex(S.prec[ew]).fillna(0).gt(0).mean():.3f}, all decoys {low.reindex(S.prec[dw]).fillna(0).gt(0).mean():.3f}')
