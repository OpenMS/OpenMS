#!/usr/bin/env python3
"""Step A, Astral: why a shadow rule is weak. For accepted entrapments (own q <= 0.01) whose library has an isobaric human
target (same charge, 10 ppm) sharing >= 2 fragments (10 ppm), where is that human partner's kept window and how good is it?
Reports |apex diff| distribution, partner q, partner score vs entrapment score. Same for counted decoys (score >= the lowest
accepted target) with such a partner, and for accepted human targets with such a partner. Diagnosis only.
Usage: isomer_partner.py <lib.pqp> <glob>"""
import sqlite3, sys, glob, numpy as np, pandas as pd
lib = sys.argv[1]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
P = pd.read_sql("""select p.ID id, p.DECOY d, p.CHARGE z, p.PRECURSOR_MZ mz, p.LIBRARY_RT lrt,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
fr = {}
for pid, mz in c.execute('select tpm.PRECURSOR_ID, tr.PRODUCT_MZ from TRANSITION tr join TRANSITION_PRECURSOR_MAPPING tpm on tpm.TRANSITION_ID=tr.ID'):
    fr.setdefault(pid, []).append(mz)
fr = {k: np.array(v) for k, v in fr.items()}
H = P[(P.d == 0) & (P.e == 0)].sort_values('mz'); hmz, hz, hid = H.mz.values, H.z.values, H.id.values
Pi = P.set_index('id')
pc = {}
def partners(pid):
    if pid in pc: return pc[pid]
    mz, z = Pi.mz[pid], Pi.z[pid]
    lo, hi = np.searchsorted(hmz, mz * (1 - 1e-5)), np.searchsorted(hmz, mz * (1 + 1e-5), side='right')
    out = []
    a = fr.get(pid, np.empty(0))
    for j in range(lo, hi):
        if hz[j] != z or hid[j] == pid: continue
        b = fr.get(hid[j], np.empty(0))
        k = int((np.abs(a[:, None] - b[None, :]) <= b[None, :] * 1e-5).any(axis=1).sum()) if len(a) and len(b) else 0
        if k >= 2: out.append((hid[j], k))
    pc[pid] = out
    return out
rows = []
for f in sorted(glob.glob(sys.argv[2])):
    t = pd.read_csv(f, sep='\t', usecols=['precursor_id', 'is_decoy', 'rt_apex', 'score', 'q_value'])
    t = t[t.score > -1e29].sort_values('score', ascending=False).drop_duplicates('precursor_id').set_index('precursor_id')
    t['e'] = Pi.e.reindex(t.index).fillna(0).values
    acc = t[(t.is_decoy == 0) & (t.q_value <= 0.01)]; thr = acc.score.min()
    grp = {'E_acc': acc[acc.e == 1].index, 'D_cnt': t[(t.is_decoy == 1) & (t.score >= thr)].index,
           'T_acc': acc[acc.e == 0].index}
    for g, ids in grp.items():
        for p in ids:
            ps = partners(p)
            if not ps: continue
            best = None
            for h, k in ps:
                if h not in t.index: continue
                dt = abs(t.rt_apex[h] - t.rt_apex[p])
                if best is None or dt < best[0]: best = (dt, t.q_value[h], t.score[h] - t.score[p], k)
            rows.append(dict(run=f, grp=g, p=p, has=best is not None, dt=best[0] if best else np.nan, hq=best[1] if best else np.nan,
                             dscore=best[2] if best else np.nan, k=best[3] if best else max(k for _, k in ps)))
x = pd.DataFrame(rows)
tot = {g: sum(1 for _ in []) for g in ('E_acc', 'D_cnt', 'T_acc')}
print(f'runs {len(glob.glob(sys.argv[2]))}; rows = group members with an isobaric human partner sharing >= 2 fragments')
for g, y in x.groupby('grp'):
    h = y[y.has]
    print(f'{g}: with partner {len(y)}; partner scored in run {len(h)}; |apex diff| <=3 s {np.mean(h.dt <= 3):.2f}, <=10 s {np.mean(h.dt <= 10):.2f}, '
          f'<=30 s {np.mean(h.dt <= 30):.2f}, median {h.dt.median():.1f} s; partner q<=0.01 {np.mean(h.hq <= 0.01):.2f}; '
          f'partner scores higher {np.mean(h.dscore > 0):.2f}; shared frags median {y.k.median():.0f}')
    co = h[h.dt <= 3]
    if len(co): print(f'    co-eluting (<=3 s): n {len(co)}, partner q<=0.01 {np.mean(co.hq <= 0.01):.2f}, partner scores higher {np.mean(co.dscore > 0):.2f}')
