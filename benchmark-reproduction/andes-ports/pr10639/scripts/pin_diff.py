"""Rows and columns that differ between two PINs of the same search (row order compared as is).
usage: pin_diff.py <a.pin> <b.pin> -> one line: rows, differing rows, differing columns with counts, peptide/score/label differences per spectrum"""
import collections, sys
A = open(sys.argv[1]).read().splitlines(); B = open(sys.argv[2]).read().splitlines()
h = A[0].split('\t'); cols = collections.Counter(); rows = 0
for x, y in zip(A[1:], B[1:]):
    if x != y:
        rows += 1
        for i, (u, v) in enumerate(zip(x.split('\t'), y.split('\t'))):
            if u != v: cols[h[i]] += 1
def psms(L):
    d = collections.defaultdict(list); hh = L[0].split('\t')
    for r in L[1:]:
        f = r.split('\t'); d[f[0]].append((f[hh.index('Peptide')], f[hh.index('score')], f[hh.index('Label')]))
    return d
a, b = psms(A), psms(B)
spec = sum(1 for k in set(a) | set(b) if a.get(k) != b.get(k))
print(f'{sys.argv[1]} vs {sys.argv[2]}: rows {len(A)-1}/{len(B)-1}, differing rows {rows}, columns {dict(cols)}, spectra with other peptide/score/label {spec}')
