"""Characterise the 'shufshare' targets (their shuffled decoy shares >=2 top-6 fragments with them; where O's excess sits,
score4_02.txt) against the rest: length, charge, distinct-residue fraction, C-terminal residue, share of entrapment."""
import sqlite3, numpy as np, pandas as pd
c = sqlite3.connect('file:merged.pqp?mode=ro', uri=True)
p = pd.read_sql('select p.ID pid, p.TRAML_ID tid, p.DECOY d, p.CHARGE z, pe.UNMODIFIED_SEQUENCE seq from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID', c)
acc = pd.read_sql('select m.PEPTIDE_ID pep, pr.PROTEIN_ACCESSION a from PEPTIDE_PROTEIN_MAPPING m join PROTEIN pr on pr.ID=m.PROTEIN_ID', c)
tr = pd.read_sql('select m.PRECURSOR_ID pid, t.PRODUCT_MZ f, t.LIBRARY_INTENSITY li, t.ID t from TRANSITION t join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID', c)
tr = tr.sort_values(['pid', 'li', 't'], ascending=[True, False, True]).groupby('pid').head(6)
F = tr.groupby('pid').f.apply(np.array).to_dict()
tp = pd.read_sql('select ID tpid, TRAML_ID tid from PRECURSOR', sqlite3.connect('file:targets.pqp?mode=ro', uri=True))
T = p[p.d == 0]; t2m = dict(zip(T.tid, T.pid)); tp2m = {a: t2m.get(b) for a, b in zip(tp.tpid, tp.tid)}
pepacc = pd.read_sql('select PRECURSOR_ID pid, PEPTIDE_ID pep from PRECURSOR_PEPTIDE_MAPPING', c).merge(acc, on='pep')
ent = pepacc.assign(e=pepacc.a.str.contains('_p_target')).groupby('pid').e.all()
Sd = p[p.tid.str.startswith('DECOYS_')]
bad = {}
for pid, tid in zip(Sd.pid, Sd.tid):
    t = tp2m.get(int(tid.split('_')[-1]))
    if t is None or pid not in F or t not in F: continue
    bad[t] = int((np.abs(F[pid][:, None] - F[t][None, :]) <= 34e-6 * F[pid][:, None]).any(1).sum()) >= 2
T = T[T.pid.isin(bad)].copy(); T['bad'] = T.pid.map(bad); T['ent'] = T.pid.map(ent).fillna(False)
T['L'] = T.seq.str.len(); T['distinct'] = T.seq.map(lambda s: len(set(s)) / len(s)); T['cterm'] = T.seq.str[-1]
T['maxrun'] = T.seq.map(lambda s: max(s.count(a) for a in set(s)) / len(s))
for b, g in T.groupby('bad'):
    print(f'shufshare={b}: n {len(g)} | entrapment {g.ent.mean():.3f} | length median {g.L.median():.0f} (<=8: {(g.L <= 8).mean():.3f}) | charge2 {(g.z == 2).mean():.3f} | distinct-residue frac {g.distinct.mean():.3f} | max single-residue frac {g.maxrun.mean():.3f} | Cterm K {(g.cterm == "K").mean():.3f} R {(g.cterm == "R").mean():.3f}')
