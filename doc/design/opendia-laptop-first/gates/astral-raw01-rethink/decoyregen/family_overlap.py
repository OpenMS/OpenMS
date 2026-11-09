"""Per decoy family: how many of a decoy's top-6 fragments fall within 34 ppm of its own target's top-6 fragments,
and the precursor m/z difference. A decoy at its target's m/z that shares fragments picks up the target's own signal."""
import sqlite3, numpy as np, pandas as pd
c = sqlite3.connect('file:merged.pqp?mode=ro', uri=True)
p = pd.read_sql('select ID pid, TRAML_ID tid, DECOY d, PRECURSOR_MZ mz from PRECURSOR', c)
tr = pd.read_sql('select m.PRECURSOR_ID pid, t.PRODUCT_MZ f, t.LIBRARY_INTENSITY li, t.ID t from TRANSITION t join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID', c)
tr = tr.sort_values(['pid', 'li', 't'], ascending=[True, False, True]).groupby('pid').head(6)
F = tr.groupby('pid').f.apply(np.array).to_dict()
tp = pd.read_sql('select ID tpid, TRAML_ID tid from PRECURSOR', sqlite3.connect('file:targets.pqp?mode=ro', uri=True))
t2m = dict(zip(p.tid[p.d == 0], p.pid[p.d == 0])); tp2m = {a: t2m.get(b) for a, b in zip(tp.tpid, tp.tid)}
mz = dict(zip(p.pid, p.mz))
rng = np.random.default_rng(0)
for tag in ('DECOY_', 'DECOYP_', 'DECOYS_'):
    q = p[p.tid.str.startswith(tag)].sample(20000, random_state=0)
    sh, dm = [], []
    for pid, tid in zip(q.pid, q.tid):
        t = tp2m.get(int(tid.split('_')[-1]))
        if t is None or pid not in F or t not in F: continue
        a, b = F[pid], F[t]
        sh.append(int((np.abs(a[:, None] - b[None, :]) <= 34e-6 * a[:, None]).any(1).sum())); dm.append(abs(mz[pid] - mz[t]))
    sh = np.array(sh); dm = np.array(dm)
    print(f'{tag:8s} n {len(sh)} shared top-6 fragments with own target: mean {sh.mean():.2f}, share>=3 {(sh >= 3).mean():.3f}, share>=4 {(sh >= 4).mean():.3f} | precursor |dmz| median {np.median(dm):.3f}, share <0.01 Th {(dm < 0.01).mean():.3f}')
