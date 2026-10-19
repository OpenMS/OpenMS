"""Percolator (-Y -U, seeds 1/42/137) on a ProSE PIN; target PSMs at q <= 0.01 and entrapment FDP (2 N_E / N, r = 1).
usage: perc_eval.py <pin> <workdir> [percolator]  -> prints JSON"""
import json, statistics, subprocess, sys
from pathlib import Path
pin, work = Path(sys.argv[1]), Path(sys.argv[2])
perc = sys.argv[3] if len(sys.argv) > 3 else '/home/user/bench/tp/Linux/x86_64/Percolator/percolator'
work.mkdir(parents=True, exist_ok=True)
n, ne = [], []
for seed in (1, 42, 137):
    t = work / f's{seed}.target.tsv'
    subprocess.run([perc, '-Y', '-U', '--seed', str(seed), '--results-psms', t, '--decoy-results-psms', work / f's{seed}.decoy.tsv', pin],
                   check=True, stdout=subprocess.DEVNULL, stderr=open(work / f's{seed}.log', 'w'))
    acc = [l.split('\t') for l in t.read_text().splitlines()[1:]]
    acc = [f for f in acc if float(f[2]) <= 0.01]
    n.append(len(acc)); ne.append(sum(all(p.startswith('ENTRAP_') for p in f[5:] if p) for f in acc))
print(json.dumps(dict(pin=str(pin), psms=n, mean_psms=round(statistics.mean(n), 1), entrap=ne,
                      fdp_combined_pct=round(200 * sum(ne) / sum(n), 2) if sum(n) else None)))
