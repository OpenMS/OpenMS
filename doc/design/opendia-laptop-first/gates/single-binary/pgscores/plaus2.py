"""Fraction true among IDs NEW adds that DIA-NN does not report at 1 % (HeLa D, seed 0), by mixture:
rate(group) = pi * rate(true) + (1 - pi) * rate(false).
- rate(true): the old (CUR) accepted IDs.
- rate(false): the same group with our apex shifted by +-30 s, so any coincidence is chance at the same library and RT density.
Indicators use a tight +-3 s apex match: DIA-NN reports the precursor (any q, q<=0.05) and OpenSWATH SVM q<=0.01."""
import sqlite3, glob, pandas as pd, numpy as np
exec(open('plaus.py').read().split("groups={")[0])
groups={'added, in DIA-NN':add[add.key.isin(D1.index)],'added, NOT in DIA-NN':add[~add.key.isin(D1.index)],'old accepted':acc[acc.precursor_id.isin(C)]}
def ind(g, shift=0.0, tol=3.0):
    d=DD.reindex(g.key); o=osw.reindex(g.precursor_id); rt=g.rt_apex.values+shift
    dn_any=np.abs(d.RT.values*60-rt)<=tol; dn05=dn_any&(d['Q.Value'].values<=0.05)
    os_=(np.abs(o.rt.values-rt)<=tol)&(o.q.values<=0.01)
    return np.array([dn_any.mean(),dn05.mean(),os_.mean()])
names=['DIA-NN any q, +-3 s','DIA-NN q<=0.05, +-3 s','OpenSWATH q<=0.01, +-3 s']
T=ind(groups['old accepted']); 
for gname in ('added, in DIA-NN','added, NOT in DIA-NN'):
    g=groups[gname]; r=ind(g); f=(ind(g,30)+ind(g,-30))/2
    print(f'\n{gname} (n={len(g)})')
    for i,n in enumerate(names):
        pi=(r[i]-f[i])/(T[i]-f[i])
        print(f'  {n:26s} group {r[i]:.3f} | true (old accepted) {T[i]:.3f} | chance (apex +-30 s) {f[i]:.3f} | est. fraction true {pi:.2f}')
