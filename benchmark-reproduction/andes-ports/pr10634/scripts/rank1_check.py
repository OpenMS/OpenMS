"""Compare the best PIN row (by 'score') per spectrum between two PINs: same peptide, label and features?"""
import sys
def best(p):
    b = {}
    with open(p) as f:
        h = f.readline().rstrip('\n').split('\t'); i_s, i_sc, i_p = h.index('ScanNr'), h.index('score'), h.index('Peptide')
        for l in f:
            c = l.rstrip('\n').split('\t')
            k = c[i_s]; s = float(c[i_sc])
            if k not in b or s > b[k][0]: b[k] = (s, c[i_p], c[1], tuple(c[3:i_p]))
    return b
a, b = best(sys.argv[1]), best(sys.argv[2])
same_pep = sum(a[k][1] == b[k][1] for k in a if k in b); same_all = sum(a[k] == b[k] for k in a if k in b)
print(f'spectra {len(a)}/{len(b)}  same rank-1 peptide {same_pep}  identical rank-1 row (label+features) {same_all}')
