"""andes 695acf7, --fragment-index-intensity-tiebreak off vs on: isotope error (andes convention: observed minus
theoretical, window -1..2) of the rank-1 PSM per spectrum (SpecId '<scan>_<scan>_<rank>_<i>' with rank 1) for all
spectra and for the spectra whose rank-1 peptide changed, and of the accepted PSMs (Percolator q <= 0.01, 3 seeds)."""
import collections, gzip, sys
from pathlib import Path
RUNS = Path('/home/user/bench/runs')


def rank1(arm, d):
    out, iso_of = {}, {}
    with gzip.open(RUNS / f'andes_{arm}_E' / d / 'native.pin.gz', 'rt') as f:
        h = next(f).rstrip('\n').split('\t'); ii, ip, il = h.index('isotope_error'), h.index('Peptide'), h.index('Label')
        for line in f:
            x = line.rstrip('\n').split('\t')
            iso_of[x[0]] = int(float(x[ii]))
            key, rank, _ = x[0].rsplit('_', 2)
            if rank == '1': out[key] = (x[ip], x[il], int(float(x[ii])))
    return out, iso_of


def dist(vals):
    c = collections.Counter(vals); n = len(vals)
    return n, ' '.join(f'{k:+d}:{100 * c[k] / n:.1f}' for k in sorted(c))


print('dataset\tarm\tsubset\tn\tisotope_error_pct\tdecoy_pct')
for d in sys.argv[1:]:
    (a, ia), (b, ib) = rank1('off', d), rank1('on', d)
    ch = [k for k in a if k in b and a[k][0] != b[k][0]]
    for arm, t, iso_of in (('off', a, ia), ('on', b, ib)):
        for name, keys in (('rank1_all', list(t)), ('rank1_changed', ch)):
            n, s = dist([t[k][2] for k in keys])
            print(f"{d}\t{arm}\t{name}\t{n}\t{s}\t{100 * sum(t[k][1] == '-1' for k in keys) / n:.1f}")
        acc = []
        for seed in (1, 42, 137):
            with gzip.open(RUNS / f'andes_{arm}_E' / d / f's{seed}' / 'target.tsv.gz', 'rt') as f:
                next(f)
                for line in f:
                    x = line.split('\t')
                    if float(x[2]) <= 0.01: acc.append(iso_of[x[0]])
        n, s = dist(acc)
        print(f"{d}\t{arm}\taccepted_3seeds\t{n}\t{s}\t")
