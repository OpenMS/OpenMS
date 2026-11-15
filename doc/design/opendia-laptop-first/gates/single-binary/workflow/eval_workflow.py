"""PLAN_SINGLE_BINARY acceptance (labels read here only): OpenDIAWorkflow <run>.precursors.tsv (one row per precursor
winner, candidate_row = search feature row = population sample_id) to a score array in population order, then
eval_scores.py's metrics. Usage: eval_workflow.py <raw01|raw02> <precursors.tsv> [...]"""
import os, subprocess, sys, numpy as np, pandas as pd
G = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/gap-20261114'
RAW, FILES = sys.argv[1], sys.argv[2:]
pop = pd.read_csv(f'{G}/armE_{RAW}/population.tsv', sep='\t', usecols=['sample_id', 'precursor_id'])
row = pd.Series(np.arange(len(pop)), index=pop.sample_id.values)
out = []
for f in FILES:
    t = pd.read_csv(f, sep='\t', usecols=['candidate_row', 'precursor_id', 'score'])
    s = np.full(len(pop), -1e30, np.float32)
    r = row.loc[t.candidate_row.values].values
    assert (pop.precursor_id.values[r] == t.precursor_id.values).all(), 'candidate_row / precursor mismatch'
    s[r] = t.score.values
    p = os.path.join(os.environ.get('JCODE_SCRATCH_DIR', '/tmp'), os.path.basename(os.path.dirname(f)) + '.npy')
    np.save(p, s); out.append(p)
subprocess.run([sys.executable, f'{G}/scripts/eval_scores.py', RAW] + out, check=True)
