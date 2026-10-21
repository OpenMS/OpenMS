"""Breakdown of wrong-site spike-in PSMs (q <= 0.01): true label residue -> reported label residue, per arm dir/seed."""
import collections, sys
sys.path.insert(0, __import__('os').path.dirname(__file__))
from mhc_site import parse, HEAVY
truth = {l.split('\t')[0]: (int(l.split('\t')[3]), l.split('\t')[2]) for l in open(sys.argv[1]).read().splitlines()[1:]}
for t in sys.argv[2:]:
    c = collections.Counter(); peps = collections.Counter()
    for line in open(t).read().splitlines()[1:]:
        f = line.split('\t')
        if float(f[2]) > 0.01: continue
        seq, mods = parse(f[4])
        if seq not in truth: continue
        h = [(p, a) for p, a, d in mods if a in HEAVY and abs(d - HEAVY[a]) < 0.01]
        if len(h) == 1 and h[0] != truth[seq]:
            c[f'{truth[seq][1]}->{h[0][1]}'] += 1; peps[seq] += 1
    print(t, dict(c), 'peptides:', dict(peps.most_common(8)))
