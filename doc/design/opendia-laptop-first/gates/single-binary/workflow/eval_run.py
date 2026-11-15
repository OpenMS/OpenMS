"""PLAN_SINGLE_BINARY: entrapment evaluation of an OpenDIAWorkflow <run>.precursors.tsv on its own population (labels read
here only). Pairs and folds come from the library by the make_population.py rules; metrics are evaluate_v2.py's (tdc,
n_at, entrapment ratio R). Usage: eval_run.py <library.pqp> <precursors.tsv> [...]"""
import sqlite3, sys, numpy as np, pandas as pd
SRC = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/pairlib-20261113/scripts/evaluate_v2.py'
LIB, FILES = sys.argv[1], sys.argv[2:]
sys.argv = [sys.argv[0]]
code = open(SRC).read().split('rng = np.random.default_rng(5)')[0]
code = code.replace("B + '/pairlib-20261113/lib/library.pqp", f"'{LIB}")
assert LIB in code
exec(code)
pre = pd.read_sql("""select p.ID precursor_id, p.TRAML_ID traml, p.DECOY decoy, pe.UNMODIFIED_SEQUENCE seq from PRECURSOR p
  join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID""", c)
target_seq = pre[pre.decoy == 0].set_index('traml').seq
pre['pair_key'] = np.where(pre.decoy == 1, pre.traml.str.replace('^DECOY_', '', regex=True), pre.traml)
pre['family'] = pre.pair_key.map(target_seq).fillna(pre.seq)
def fnv(s):
    h = 0xcbf29ce484222325
    for b in s.encode(): h = ((h ^ b) * 0x100000001b3) & 0xFFFFFFFFFFFFFFFF
    return h
pre['fold'] = pre.family.map(lambda s: fnv(s) % 3)
pre = pre.set_index('precursor_id')
for f in FILES:
    t = pd.read_csv(f, sep='\t', usecols=['candidate_row', 'precursor_id', 'is_decoy', 'score', 'q_value'])
    pop = pd.DataFrame({'sample_id': t.candidate_row.values, 'precursor_id': t.precursor_id.values, 'is_decoy': t.is_decoy.values,
                        'fold': pre.fold.loc[t.precursor_id].values, 'pair_key': pre.pair_key.loc[t.precursor_id].values})
    w, best = tdc(pop, t.score.values.astype(np.float64))
    tw = w[w.is_decoy == 0].assign(e=lambda d: d.precursor_id.map(ent).values)
    tb = best[best.is_decoy == 0].assign(e=lambda d: d.precursor_id.map(ent).values)
    n1 = int((tb.q <= 0.01).sum()); f1 = tb[tb.q <= 0.01].e.sum() * (1 + 1 / R) / max(n1, 1)
    # The tool's own q-values (independent check of the in-binary TDC).
    tt = t[(t.is_decoy == 0) & (t.q_value <= 0.01)]
    own = tt.precursor_id.map(ent).sum() * (1 + 1 / R) / max(len(tt), 1)
    print(f'{f.split("/")[-1]:48s} N@FDP_c 0.5/1/2%: {[n_at(tw.s.values, tw.e.values, x) for x in (0.005, 0.01, 0.02)]} | '
          f'q<=0.01 {n1} FDP_c {f1:.4f} | tool q<=0.01 {len(tt)} FDP_c {own:.4f}', flush=True)
