"""Per spectrum, whether two PINs of the same search differ only by the peptide tie order. A spectrum whose hits
changed is explained if every peptide in only one PIN's hits of that spectrum has a tie
partner, i.e. another peptidoform of the same protein with the same residue composition and number of Oxidation (M)
(another sequence, or the same sequence with the oxidation on another M)
(Trypsin/P, <= 2 missed cleavages, length 7-40, same terminal modifications)? Then the two have the same precursor
m/z and protein, and the order between them decided which one stayed in the candidate set), or is the
lowest-scoring hit of that spectrum in its PIN (moved over the report:top_hits cut by such a swap).
usage: tie_partner.py <fasta> <old.pin> <new.pin>"""
import collections, re, sys
fasta, a, b = sys.argv[1:4]
seqs = {}; acc = None
for l in open(fasta):
    if l.startswith('>'): acc = l[1:].split()[0]; seqs[acc] = []
    else: seqs[acc].append(l.strip())
seqs = {k: ''.join(v) for k, v in seqs.items()}
def peptides(s):
    cut = [0] + [i + 1 for i, c in enumerate(s) if c in 'KR' and i + 1 < len(s)] + [len(s)]
    for i in range(len(cut) - 1):
        for j in range(i + 1, min(i + 4, len(cut))):
            p = s[cut[i]:cut[j]]
            if 7 <= len(p) <= 40: yield p
    if s.startswith('M'):  # initiator methionine removed
        yield from (p for p in peptides(s[1:]) if s[1:].startswith(p))
def key(pep):
    core = pep.split('.', 1)[1].rsplit('.', 1)[0]
    ox = core.count('M[+15.99'); res = re.sub(r'n?\[.*?\]', '', core)
    return ''.join(sorted(res)), ox, core.startswith('n[')
def read(p):
    L = open(p).read().splitlines(); h = L[0].split('\t'); ip, ipr = h.index('Peptide'), h.index('Proteins')
    s = collections.defaultdict(list)
    for r in L[1:]:
        c = r.split('\t'); s[(c[0], c[h.index('ExpMass')])].append((c[ip], tuple(c[ipr:]), float(c[h.index('score')]), r))
    return s
A, B = read(a), read(b); n = collections.Counter(); bad = []
for k in set(A) | set(B):
    if [x[3] for x in A.get(k, [])] == [x[3] for x in B.get(k, [])]: continue
    n['changed'] += 1
    pa = {x[0]: x for x in A.get(k, [])}; pb = {x[0]: x for x in B.get(k, [])}
    if pa.keys() == pb.keys(): n['same peptides'] += 1; continue
    ok = True
    for src, other in ((pa, pb), (pb, pa)):
        low = min(x[2] for x in src.values())
        for pep in set(src) - set(other):
            kk = key(pep); core = re.sub(r'n?\[.*?\]', '', pep.split('.', 1)[1].rsplit('.', 1)[0])
            partner = any(''.join(sorted(p)) == kk[0] and p.count('M') >= kk[1] and (p != core or p.count('M') > kk[1])
                          for prot in src[pep][1] for p in set(peptides(seqs.get(prot, ''))))
            ok &= partner or src[pep][2] == low
    n['explained' if ok else 'unexplained'] += 1
    if not ok: bad.append(k[0])
print(b.split('/runs/')[-1], dict(n), bad[:5])
