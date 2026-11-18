"""Foreign-proteome entrapment FASTA (b1, label-free construction, inputs in inputs.sha256).
Targets: every unique unmodified target peptide of the HeLa experimental library, one entry each (as entlib v1).
Entrapments: Arabidopsis thaliana reviewed (UP000006548) fully tryptic peptides (cleave after K/R, not before P? no:
Trypsin/P as DIALibGen), 0-2 missed cleavages, length 7-44, that do not occur in any human reviewed protein (I=L) and are
not a HeLa library peptide (I=L). Sampled without replacement, matched to the HeLa peptides' length distribution
(per-length counts equal; lengths with too few candidates take from the nearest length), seeded."""
import sqlite3, random, re, collections
L=__import__('os').environ['HELA_PQP']
HUM=__import__('os').environ['HUM_FASTA']
def read(p):
    out=[]; s=[]
    for l in open(p):
        if l.startswith('>'):
            if s: out.append(''.join(s)); s=[]
        else: s.append(l.strip())
    if s: out.append(''.join(s))
    return out
def digest(seq, mc=2, lo=7, hi=44):
    sites=[i+1 for i,a in enumerate(seq) if a in 'KR']
    b=[0]+[x for x in sites if x<len(seq)]+[len(seq)]
    b=sorted(set(b)); peps=set()
    for i in range(len(b)-1):
        for j in range(i+1,min(i+2+mc,len(b))):
            p=seq[b[i]:b[j]]
            if lo<=len(p)<=hi and set(p)<=set('ACDEFGHIKLMNPQRSTVWY'): peps.add(p)
    return peps
c=sqlite3.connect(f'file:{L}?mode=ro',uri=True)
hela=sorted({r[0] for r in c.execute("select UNMODIFIED_SEQUENCE from PEPTIDE where DECOY=0")})
norm=lambda p:p.replace('I','L')
human=set()
for s in read(HUM): human|={norm(p) for p in digest(s)}
human|={norm(p) for p in hela}
# Also exclude any Arabidopsis peptide that is a substring of a human protein (I=L), to drop semi-specific human matches
humcat='|'.join(norm(s) for s in read(HUM))
cand=collections.defaultdict(set)
for s in read('ath.fasta'):
    for p in digest(s):
        if norm(p) not in human: cand[len(p)].add(p)
rng=random.Random(20261118)
need=collections.Counter(len(p) for p in hela); out=[]; used=set(); short=0
for ln in sorted(need):
    pool=sorted(cand[ln]-used); k=need[ln]
    if len(pool)<k:
        short+=k-len(pool)
        for d in range(1,10):
            for l2 in (ln-d,ln+d):
                pool+=sorted(cand.get(l2,set())-used-set(pool))
            if len(pool)>=k: break
    pick=rng.sample(pool,k); used|=set(pick); out+=pick
out=[p for p in out if norm(p) not in humcat]  # substring filter
print(f'HeLa targets {len(hela)}, Arabidopsis candidates {sum(map(len,cand.values()))}, picked {len(out)} (length-short {short}, after substring filter)')
with open('hela_ath.fasta','w') as f:
    for p in hela: f.write(f'>sp|{p}_target|{p}_target\n{p}\n')
    for p in out: f.write(f'>sp|{p}_p_target|{p}_ARATH_entrap\n{p}\n')
