"""How many accepted PSMs (q <= 0.01) are not the spectrum's search-engine rank-1 row, and how many of those are a
positional isomer of the rank-1 peptide (same residues, other site)?  usage: rank_switch.py input.pin target.tsv"""
import re, sys
def strip(p): return re.sub(r'\[[^\]]*\]|\([^)]*\)', '', p.split('.', 1)[1].rsplit('.', 1)[0])
best = {}
with open(sys.argv[1]) as f:
    h = f.readline().rstrip('\n').split('\t'); i_id, i_sc, i_p = h.index('SpecId'), h.index('score'), h.index('Peptide')
    for l in f:
        c = l.rstrip('\n').split('\t'); s = float(c[i_sc])
        if c[i_id] not in best or s > best[c[i_id]][0]: best[c[i_id]] = (s, c[i_p])
n = sw = iso = 0
with open(sys.argv[2]) as f:
    h = f.readline().rstrip('\n').split('\t')
    for l in f:
        c = l.rstrip('\n').split('\t')
        if float(c[2]) > 0.01: continue
        n += 1; r1 = best[c[0]][1]
        if c[4] != r1:
            sw += 1; iso += strip(c[4]) == strip(r1)
print(f'accepted {n}  not rank-1 {sw}  of which positional isomer of rank-1 {iso}')
