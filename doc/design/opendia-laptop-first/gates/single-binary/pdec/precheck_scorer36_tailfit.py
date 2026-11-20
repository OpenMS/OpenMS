"""PRECHECK_SCORER36 control: tail-focused learners. Negatives = hardest decoy windows only
(decoy precursors whose best tool-score window is in the top 2 % of decoys, all their windows).
  C_tail: DIA-NN-labelled positives (ceiling, not a method)
  B_tail: label-free positives (tool q<=0.01 winners), 3 rounds
Diagnostic only. Usage: precheck_scorer36_tailfit.py"""
import numpy as np, pandas as pd, runpy, sys
if len(sys.argv) < 3: sys.argv = ['x', '/home/singjc/.jcode/scratch/dia2/gap7/D36f.tsv.gz', 'D36']
NEG_Q = float(sys.argv[3]) if len(sys.argv) > 3 else 0.98
# Reuse loaders/evaluator from precheck_scorer36 without its main body: exec the top half.
src = open('/home/singjc/.jcode/scratch/dia2/gap7/precheck_scorer36.py').read().split('kept = F.score.values')[0]
g = {}; exec(compile(src, 'precheck_scorer36_head', 'exec'), g)
F, pid, dec, gbm, evaluate, winners, d = g['F'], g['pid'], g['dec'], g['gbm'], g['evaluate'], g['winners'], g['d']
kept = F.score.values > -1e29
s_tool = np.where(kept, F.score.values, -1e30).astype(np.float64)
best = pd.Series(s_tool).groupby(pid).max()
dbest = best[best.index.isin(np.unique(pid[dec]))]
hard = dbest[dbest >= dbest.quantile(NEG_Q)].index
neg = dec & np.isin(pid, hard)
print(f'hard decoy precursors {len(hard)}, windows {neg.sum()}', flush=True)
evaluate(s_tool, 'A tool score')
drt = np.abs(F.rt_apex.values - (d.RT * 60).reindex(pid).values)
sC = gbm((~dec) & (drt <= 6), neg, seed=21)
evaluate(sC, 'C_tail DIA-NN-labelled vs hard decoys')
s = s_tool
for it in range(3):
    s = gbm(winners(s), neg, seed=30 + it)
    evaluate(s, f'B_tail label-free vs hard decoys round {it + 1}')
