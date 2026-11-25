#!/usr/bin/env python3
"""Step A: does borrowing explain the entrapment-vs-own-decoy pair asymmetry? Diagnosis only.
Per run: donors = accepted human IDs (non-entrapment target, q <= 0.01). For each entrapment pair whose better member has
q <= 0.05, each member gets share = max #library fragments (FTOL ppm) shared with a donor (not itself) whose apex is within
DT s of the member's apex and precursor m/z within W Th (-1: no donor). Then P(entrapment scores above own decoy) by
(a) entrapment share class, (b) max(share_E, share_D) <= 1 (neither borrows) vs not, and by best-q band.
Also the q <= 0.01 count accounting: accepted E and counted decoys D_mt split by share class.
Usage: pair_borrow.py <lib.pqp> <W Th> <DT s> <FTOL ppm> <arm>=<glob> ..."""
import sqlite3, sys, glob, numpy as np, pandas as pd
lib, W, DT, FT = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]) * 1e-6
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
P = pd.read_sql("""select p.ID id, p.TRAML_ID t, p.DECOY d, p.PRECURSOR_MZ mz,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
P['key'] = np.where(P.d == 1, P.t.str.replace('^DECOY_', '', regex=True), P.t)
import os
if os.environ.get('RELABEL'):  # entrapment-labelled exact human substrings count as human (step A label audit)
    seqs = dict(c.execute("select p.ID, pe.UNMODIFIED_SEQUENCE from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID").fetchall())
    rl = set(pd.read_csv(os.environ['RELABEL'], sep='\t').seq)
    m = (P.d == 0) & P.id.map(seqs).isin(rl)
    P.loc[m, 'e'] = 0
    print(f'relabelled as human: {int(m.sum())}', flush=True)
T = P[P.d == 0]; r = (T.e == 1).sum() / (T.e == 0).sum()
pairs = T[T.e == 1].set_index('key')[['id']].join(P[P.d == 1].set_index('key')[['id']].rename(columns={'id': 'did'}), how='inner')
ent = P.set_index('id').e; mzs = P.set_index('id').mz
fr = {}
def load(ids):
    need = [i for i in ids if i not in fr]
    for k in range(0, len(need), 20000):
        ch = need[k:k + 20000]
        for i in ch: fr[i] = []
        for pid, mz in c.execute(f"select tpm.PRECURSOR_ID, tr.PRODUCT_MZ from TRANSITION tr join TRANSITION_PRECURSOR_MAPPING tpm on tpm.TRANSITION_ID=tr.ID where tpm.PRECURSOR_ID in ({','.join(map(str, ch))})"):
            fr[pid].append(mz)
    for i in need: fr[i] = np.array(fr[i])
def share(pid, rt, mz, drt, dmz, did):
    lo, hi = np.searchsorted(drt, rt - DT), np.searchsorted(drt, rt + DT, side='right')
    j = np.arange(lo, hi); j = j[(np.abs(dmz[j] - mz) <= W) & (did[j] != pid)]
    if len(j) == 0: return -1
    a = fr[pid]; best = 0
    for d in did[j]:
        b = fr[d]
        if len(a) and len(b): best = max(best, int((np.abs(a[:, None] - b[None, :]) <= b[None, :] * FT).any(axis=1).sum()))
    return best
cl = lambda s: np.where(s <= 1, '0-1', np.where(s <= 3, '2-3', '4+'))
for spec in sys.argv[5:]:
    arm, g = spec.split('=', 1)
    rows = []; acct = []
    for f in sorted(glob.glob(g)):
        t = pd.read_csv(f, sep='\t', usecols=['precursor_id', 'is_decoy', 'rt_apex', 'score', 'q_value'])
        t = t[t.score > -1e29].sort_values('score', ascending=False).drop_duplicates('precursor_id').set_index('precursor_id')
        t['e'] = ent.reindex(t.index).fillna(0).values; t['mz'] = mzs.reindex(t.index).values
        acc = t[(t.is_decoy == 0) & (t.q_value <= 0.01)]; thr = acc.score.min()
        don = acc[acc.e == 0].sort_values('rt_apex'); drt, dmz, did = don.rt_apex.values, don.mz.values, don.index.values
        x = pairs[pairs.id.isin(t.index) & pairs.did.isin(t.index)].copy()
        x['st'] = t.score.reindex(x.id).values; x['sd'] = t.score.reindex(x.did).values
        x['qb'] = np.fmin(t.q_value.reindex(x.id).values, t.q_value.reindex(x.did).values)
        x = x[x.qb <= 0.05]
        dmt = t[(t.is_decoy == 1) & (t.score >= thr)]
        load(set(did) | set(x.id) | set(x.did) | set(dmt.index) | set(acc.index))
        x['sE'] = [share(p, t.rt_apex[p], t.mz[p], drt, dmz, did) for p in x.id]
        x['sD'] = [share(p, t.rt_apex[p], t.mz[p], drt, dmz, did) for p in x.did]
        rows.append(x)
        ae = acc[acc.e == 1]
        acct.append(pd.DataFrame({'grp': ['E'] * len(ae) + ['D'] * len(dmt),
                                  's': [share(p, t.rt_apex[p], t.mz[p], drt, dmz, did) for p in ae.index] +
                                       [share(p, t.rt_apex[p], t.mz[p], drt, dmz, did) for p in dmt.index]}))
    x = pd.concat(rows); a = pd.concat(acct)
    x['win'] = x.st > x.sd
    print(f'== {arm}: pairs {len(x)} P {x.win.mean():.3f}  (W {W} Th, DT {DT} s, FTOL {FT * 1e6:.0f} ppm)', flush=True)
    for band, m in (('all q<=0.05', x.qb <= 0.05), ('q<=0.001', x.qb <= 0.001), ('0.001-0.01', (x.qb > 0.001) & (x.qb <= 0.01)), ('0.01-0.05', x.qb > 0.01)):
        y = x[m]
        neither = y[(y.sE <= 1) & (y.sD <= 1)]
        parts = [f'neither borrows P {neither.win.mean():.3f} (n {len(neither)})']
        for k in ('0-1', '2-3', '4+'):
            z = y[cl(y.sE.values) == k]
            parts.append(f'E share {k}: P {z.win.mean():.3f} (n {len(z)})')
        for k in ('0-1', '2-3', '4+'):
            z = y[cl(y.sD.values) == k]
            parts.append(f'D share {k}: P {z.win.mean():.3f} (n {len(z)})')
        print(f'  {band:11s} n {len(y)} P {y.win.mean():.3f} | ' + ' | '.join(parts), flush=True)
    a['c'] = cl(a.s.values)
    ct = a.groupby(['grp', 'c']).size().unstack(fill_value=0)
    nE, nD = (a.grp == 'E').sum(), (a.grp == 'D').sum()
    exp = r / (1 + r)
    line = f'  q<=0.01 accounting: E {nE}, counted decoys {nD}, expected E if exchangeable D r/(1+r) = {nD * exp:.0f}; by share class '
    line += ' | '.join(f'{k}: E {ct.loc["E", k] if k in ct.columns else 0} vs D*r/(1+r) {(ct.loc["D", k] if k in ct.columns else 0) * exp:.0f}' for k in ('0-1', '2-3', '4+'))
    print(line, flush=True)
