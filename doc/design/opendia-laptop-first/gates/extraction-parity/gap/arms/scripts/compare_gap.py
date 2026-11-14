"""PLAN_GAP comparison against the control (evaluation only; reads entrapment labels like evaluate_v2.py).
For each arm and raw: N@FDP_c 0.5/1/2 %, q <= 0.01 N and FDP_c, and the paired bootstrap 90 % interval of
Delta N@1 % and Delta N@2 % against the control arm (1,000 precursor resamples, seed 5, as evaluate_v2.py).
Usage: compare_gap.py <control> <arm> [...]   (run with the PAIRLIB scripts dir as cwd-independent import)."""
import sys, numpy as np, pandas as pd
SRC = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/pairlib-20261113/scripts/evaluate_v2.py'
CTRL, GAP_ARMS = sys.argv[1], sys.argv[2:]
sys.argv = [sys.argv[0]]
exec(open(SRC).read().split('rng = np.random.default_rng(5)')[0])  # defines its own ARMS; ours is GAP_ARMS (amendment 1)
ARMS = GAP_ARMS
rng = np.random.default_rng(5)
for raw in ('raw01', 'raw02'):
    pop = pd.read_csv(f'{O}/{raw}/population.tsv', sep='\t')
    res = {}
    for arm in [CTRL] + ARMS:
        w, best = tdc(pop, read(raw, arm, pop))
        tw = w[w.is_decoy == 0].assign(e=lambda d: d.precursor_id.map(ent).values)
        tb = best[best.is_decoy == 0].assign(e=lambda d: d.precursor_id.map(ent).values)
        n1 = int((tb.q <= 0.01).sum()); f1 = tb[tb.q <= 0.01].e.sum() * (1 + 1 / R) / max(n1, 1)
        res[arm] = tw.set_index('precursor_id')[['s', 'e']]
        line = f'{raw} {arm:10s} N@FDP_c 0.5/1/2%: {[n_at(tw.s.values, tw.e.values, t) for t in (0.005, 0.01, 0.02)]} | q<=0.01 {n1} FDP_c {f1:.4f}'
        if arm != CTRL:
            m = res[arm].join(res[CTRL].s.rename('b'), how='inner')
            a_, b_, e_ = m.s.values, m.b.values, m.e.values.astype(int); d1, d2 = [], []
            for _ in range(1000):
                i = rng.integers(0, len(m), len(m))
                d1.append(n_at(a_[i], e_[i], 0.01) - n_at(b_[i], e_[i], 0.01)); d2.append(n_at(a_[i], e_[i], 0.02) - n_at(b_[i], e_[i], 0.02))
            line += f' | vs {CTRL} 90%: @1% {np.percentile(d1, 5):+.0f}..{np.percentile(d1, 95):+.0f} @2% {np.percentile(d2, 5):+.0f}..{np.percentile(d2, 95):+.0f}'
        print(line, flush=True)
