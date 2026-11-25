#!/usr/bin/env python3
"""Step A, HeLa: do pseudo-reverse predicted decoys differ from their own entrapment in LIBRARY properties, and does the
pair outcome P (entrapment scores above own decoy, best q <= 0.05) follow that difference? Diagnosis only.
Per pair: delta = decoy - entrapment for library RT, IM, #transitions, top-fragment intensity share, mean fragment m/z,
fraction of y ions, and near-neighbour fragment sharing nb_win (diag_nb definition, W 12.5 Th, RT 0.015, IM 0.03, 15 ppm),
computed only for pairs that are competitive in at least one run. P by tercile of each delta, pooled over runs.
Also P by best-member q band. Usage: pair_lib.py <lib.pqp> <arm>=<glob> ..."""
import sqlite3, sys, glob, numpy as np, pandas as pd
lib = sys.argv[1]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
P = pd.read_sql("""select p.ID id, p.TRAML_ID t, p.DECOY d, p.CHARGE z, p.PRECURSOR_MZ mz, p.LIBRARY_RT rt, p.LIBRARY_DRIFT_TIME im,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
P['key'] = np.where(P.d == 1, P.t.str.replace('^DECOY_', '', regex=True), P.t)
T = P[P.d == 0]
r = (T.e == 1).sum() / (T.e == 0).sum()
pairs = T[T.e == 1].set_index('key').join(P[P.d == 1].set_index('key')[['id', 'rt', 'im']].rename(
    columns={'id': 'did', 'rt': 'drt', 'im': 'dim'}), how='inner').reset_index()
runs = []
for spec in sys.argv[2:]:
    arm, g = spec.split('=', 1)
    for f in sorted(glob.glob(g)):
        t = pd.read_csv(f, sep='\t', usecols=['precursor_id', 'score', 'q_value'])
        t = t[t.score > -1e29].sort_values('score', ascending=False).drop_duplicates('precursor_id').set_index('precursor_id')
        st, sd = t.score.reindex(pairs.id).values, t.score.reindex(pairs.did).values
        qt, qd = t.q_value.reindex(pairs.id).values, t.q_value.reindex(pairs.did).values
        runs.append((arm, st, sd, np.fmin(qt, qd)))
comp = np.zeros(len(pairs), bool)
for _, st, sd, qb in runs: comp |= ~np.isnan(st) & ~np.isnan(sd) & (qb <= 0.05)
x = pairs[comp].copy()
print(f'library r {r:.4f}; entrapment pairs {len(pairs)}, competitive in >= 1 run {len(x)}', flush=True)
ids = np.concatenate([x.id.values, x.did.values])
q = ','.join(map(str, ids))
tr = pd.read_sql(f"""select tpm.PRECURSOR_ID id, tr.PRODUCT_MZ f, tr.LIBRARY_INTENSITY i, tr.TYPE ty
  from TRANSITION tr join TRANSITION_PRECURSOR_MAPPING tpm on tpm.TRANSITION_ID=tr.ID where tpm.PRECURSOR_ID in ({q})""", c)
g = tr.groupby('id')
L = pd.DataFrame({'nt': g.size(), 'top': g.i.max() / g.i.sum(), 'fmz': g.f.mean(), 'yfrac': g.ty.apply(lambda s: (s == 'y').mean())})
for k in L.columns:
    x['d_' + k] = L[k].reindex(x.did).values - L[k].reindex(x.id).values
x['d_rt'] = x.drt - x.rt; x['d_im'] = x.dim - x.im
x['d_absrt'] = np.abs(x.d_rt)
# library RT density of all targets around each member (crowding proxy): #targets within +-0.005 RT and +-0.02 IM
tg = T.sort_values('rt'); trt, tim = tg.rt.values, tg.im.values
def dens(rt, im):
    lo, hi = np.searchsorted(trt, rt - 0.005), np.searchsorted(trt, rt + 0.005)
    return int((np.abs(tim[lo:hi] - im) <= 0.02).sum())
x['dens_e'] = [dens(a, b) for a, b in zip(x.rt, x.im)]
x['dens_d'] = [dens(a, b) for a, b in zip(x.drt, x.dim)]
x['d_dens'] = np.log((x.dens_d + 1) / (x.dens_e + 1))
idx = np.where(comp)[0]
print('library deltas (decoy - entrapment), competitive pairs: ' + ', '.join(
    f'{k} median {x[k].median():.4g} mean {x[k].mean():.4g}' for k in ['d_rt', 'd_im', 'd_nt', 'd_top', 'd_fmz', 'd_yfrac', 'd_dens']), flush=True)
for arm in dict.fromkeys(a for a, *_ in runs):
    W = np.zeros(len(x)); N = np.zeros(len(x))
    band = {b: [0, 0] for b in ('q<=0.001', '0.001-0.01', '0.01-0.05')}
    for a, st, sd, qb in runs:
        if a != arm: continue
        st, sd, qb = st[idx], sd[idx], qb[idx]
        ok = ~np.isnan(st) & ~np.isnan(sd) & (qb <= 0.05)
        N += ok; W += ok & (st > sd)
        for b, m in (('q<=0.001', qb <= 0.001), ('0.001-0.01', (qb > 0.001) & (qb <= 0.01)), ('0.01-0.05', (qb > 0.01) & (qb <= 0.05))):
            mm = ok & m; band[b][0] += int((mm & (st > sd)).sum()); band[b][1] += int(mm.sum())
    print(f'== {arm}: P {W.sum() / N.sum():.3f} (n {int(N.sum())}); by best q: ' + ', '.join(
        f'{b} {w / max(n, 1):.3f} (n {n})' for b, (w, n) in band.items()), flush=True)
    for k in ['d_rt', 'd_absrt', 'd_im', 'd_nt', 'd_top', 'd_fmz', 'd_yfrac', 'd_dens']:
        v = x[k].values
        cut = np.nanquantile(v, [1 / 3, 2 / 3])
        parts = []
        for lab, m in (('low', v <= cut[0]), ('mid', (v > cut[0]) & (v <= cut[1])), ('high', v > cut[1])):
            parts.append(f'{lab} P {W[m].sum() / max(N[m].sum(), 1):.3f} (n {int(N[m].sum())})')
        print(f'  {k:8s} cuts {cut[0]:.4g}/{cut[1]:.4g}: ' + ' | '.join(parts), flush=True)
