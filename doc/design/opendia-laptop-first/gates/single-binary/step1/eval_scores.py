"""PLAN_SINGLE_BINARY step 1 evaluation (labels read here only): N@FDP_c 0.5/1/2 % and q <= 0.01 N, FDP_c of score
arrays on the arm-E populations, with evaluate_v2.py's tdc / n_at / entrapment ratio. Usage: eval_scores.py <raw> <npy> [...]"""
import sys, numpy as np, pandas as pd
SRC = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/pairlib-20261113/scripts/evaluate_v2.py'
RAW, FILES = sys.argv[1], sys.argv[2:]
sys.argv = [sys.argv[0]]
G = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/gap-20261114'
# On codon the gap root is bound over the pairlib root, so its lib/library.pqp is the arm-E library (md5 e0461beb).
code = open(SRC).read().split('rng = np.random.default_rng(5)')[0]
code = code.replace("B + '/pairlib-20261113/lib/library.pqp", f"'{G}/armE_raw01/library.pqp")
assert 'armE_raw01/library.pqp' in code
exec(code)
pop = pd.read_csv(f'{G}/armE_{RAW}/population.tsv', sep='\t')
for f in FILES:
    w, best = tdc(pop, np.load(f).astype(np.float64))
    tw = w[w.is_decoy == 0].assign(e=lambda d: d.precursor_id.map(ent).values)
    tb = best[best.is_decoy == 0].assign(e=lambda d: d.precursor_id.map(ent).values)
    n1 = int((tb.q <= 0.01).sum()); f1 = tb[tb.q <= 0.01].e.sum() * (1 + 1 / R) / max(n1, 1)
    print(f'{RAW} {f.split("/")[-1]:28s} N@FDP_c 0.5/1/2%: {[n_at(tw.s.values, tw.e.values, t) for t in (0.005, 0.01, 0.02)]} | q<=0.01 {n1} FDP_c {f1:.4f}', flush=True)
