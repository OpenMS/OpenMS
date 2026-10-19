"""Rank-1 PSMs per spectrum (highest PIN 'score'): isotope_error distribution and target/decoy label, none vs int,
for all spectra and for the spectra whose rank-1 peptide changed."""
import gzip, sys, collections
from pathlib import Path
RUNS = Path('/home/user/bench/runs')
def top1(arm, d):
    best = {}
    with gzip.open(RUNS / arm / d / 'search.pin.tsv.gz', 'rt') as f:
        h = next(f).rstrip('\n').split('\t'); i_s, i_p, i_l, i_i = h.index('score'), h.index('Peptide'), h.index('Label'), h.index('isotope_error')
        for line in f:
            x = line.rstrip('\n').split('\t'); sc = float(x[i_s])
            if x[0] not in best or sc > best[x[0]][0]: best[x[0]] = (sc, x[i_p], x[i_l], x[i_i])
    return best
print('dataset\tarm\tsubset\tn\tiso0_pct\tdecoy_pct')
for d in sys.argv[1:]:
    a, b = top1('none', d), top1('int', d)
    ch = [k for k in a if k in b and a[k][1] != b[k][1]]
    for arm, t in (('none', a), ('int', b)):
        for name, keys in (('all', list(t)), ('changed', ch)):
            n = len(keys)
            print(f"{d}\t{arm}\t{name}\t{n}\t{100*sum(float(t[k][3])==0 for k in keys)/n:.1f}\t{100*sum(t[k][2]=='-1' for k in keys)/n:.1f}")
