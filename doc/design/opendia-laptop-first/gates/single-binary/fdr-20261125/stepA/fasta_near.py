#!/usr/bin/env python3
"""Step A: are the cross-acquisition recurrent entrapment peptides near-identical to a human or contaminant sequence the
homolog-free rule missed? Diagnosis only. The rule excluded exact substrings of human proteins and Hamming <= 2 to
same-length human tryptic peptides (I=L). Here, against every same-length window of every protein in FASTA (human reviewed
plus contaminants plus iRT), I=L and also with deamidation-equivalent N=D, Q=E: minimum Hamming distance, and whether a
window differs only by isobaric or near-isobaric swaps (K/Q 0.036 Da, N/GG, Q/AG, W/DA ...: approximated by composition
equality). Compared with a random sample of non-recurrent library entrapment peptides of matched length.
Usage: fasta_near.py <fasta> <recurrent.tsv> <lib.pqp> <n random>"""
import sys, sqlite3, random, numpy as np
fa, rec, lib, NR = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
seqs = []; s = []
for l in open(fa):
    if l.startswith('>'):
        if s: seqs.append(''.join(s)); s = []
    else: s.append(l.strip())
if s: seqs.append(''.join(s))
norm1 = lambda p: p.replace('I', 'L')
norm2 = lambda p: p.replace('I', 'L').replace('N', 'D').replace('Q', 'E')
R = [l.strip() for l in open(rec)][1:]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
ent = sorted({r[0] for r in c.execute("""select distinct pe.UNMODIFIED_SEQUENCE from PEPTIDE pe
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID
  where instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 and pe.DECOY=0""")})
rng = random.Random(0)
bylen = {}
for p in ent:
    if p not in R: bylen.setdefault(len(p), []).append(p)
ctrl = []
for p in R:
    pool = bylen.get(len(p), [])
    ctrl += rng.sample(pool, min(NR, len(pool)))
def windows(L, norm):
    arr = []
    for q in seqs:
        q = norm(q)
        for i in range(len(q) - L + 1): arr.append(q[i:i + L])
    return np.frombuffer(''.join(arr).encode(), dtype=np.uint8).reshape(-1, L) if arr else None
cache = {}
def best(p, norm, key):
    L = len(p)
    if (L, key) not in cache: cache[(L, key)] = windows(L, norm)
    M = cache[(L, key)]
    q = np.frombuffer(norm(p).encode(), dtype=np.uint8)
    h = (M != q).sum(1)
    comp = np.sort(q)
    j = np.where(h <= 4)[0]
    same_comp = any((np.sort(M[k]) == comp).all() for k in j)
    return int(h.min()), same_comp
for lab, S in (('recurrent', R), ('control', ctrl)):
    v1 = [best(p, norm1, 1) for p in S]; v2 = [best(p, norm2, 2) for p in S]
    h1 = np.array([a for a, _ in v1]); h2 = np.array([a for a, _ in v2]); sc = np.array([b for _, b in v1])
    print(f'{lab}: n {len(S)}; min Hamming (I=L) to any FASTA window: <=2 {np.mean(h1 <= 2):.2f} <=3 {np.mean(h1 <= 3):.2f} '
          f'<=4 {np.mean(h1 <= 4):.2f} median {np.median(h1):.0f}; with N=D, Q=E <=2 {np.mean(h2 <= 2):.2f} <=3 {np.mean(h2 <= 3):.2f}; '
          f'same composition within Hamming 4 {np.mean(sc):.2f}', flush=True)
    if lab == 'recurrent':
        for p, (a, b), (a2, _) in zip(S, v1, v2): print(f'   {p}: H {a}, H(N=D,Q=E) {a2}, same-comp {b}')
