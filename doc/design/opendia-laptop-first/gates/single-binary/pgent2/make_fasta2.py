"""PLAN_PGENT2 entrapment FASTA. As make_fasta.py (Arabidopsis UP000006548 reviewed, tryptic 0-2 MC, length 7-44, not a
human tryptic peptide or human-protein substring, I=L), plus the rule fixed by PRECHECK_FOREIGN: drop every candidate within
Hamming distance <= 2 (I=L) of any same-length human tryptic peptide (0-2 MC, 7-44) or HeLa library peptide.
Length-matched to the HeLa peptides, seed 20261119 (different from the precheck's 20261118)."""
import sqlite3, random, collections, numpy as np
src=open('make_fasta.py').read()
exec(src.split("c=sqlite3")[0])
L='/home/singjc/Documents/github/testing_auto_irt/diaPASEF_evosep_published/200623_hela_evosep_fractions_library_ptypic_decoy.pqp'
HUM='/home/singjc/Documents/github/testing_auto_irt/diaPASEF_evosep_published/uniprotkb_Human_AND_reviewed_true_AND_m_2026_02_26.fasta'
c=sqlite3.connect(f'file:{L}?mode=ro',uri=True)
hela=sorted({r[0] for r in c.execute("select UNMODIFIED_SEQUENCE from PEPTIDE where DECOY=0")})
norm=lambda p:p.replace('I','L')
hum=collections.defaultdict(set)
for s in read(HUM):
    for p in digest(s): hum[len(p)].add(norm(p))
for p in hela: hum[len(p)].add(norm(p))
A={ln:np.frombuffer(''.join(sorted(v)).encode(),dtype=np.uint8).reshape(-1,ln) for ln,v in hum.items()}
humcat='|'.join(norm(s) for s in read(HUM))
def far(p):
    M=A.get(len(p))
    if M is None: return True
    q=np.frombuffer(norm(p).encode(),dtype=np.uint8)
    return int((M!=q).sum(1).min())>2
cand=collections.defaultdict(list)
for s in read('ath.fasta'):
    for p in digest(s):
        if norm(p) not in hum[len(p)]: cand[len(p)].append(p)
cand={k:sorted(set(v)) for k,v in cand.items()}
rng=random.Random(20261119); need=collections.Counter(len(p) for p in hela); out=[]; used=set()
for ln in sorted(need):
    k=need[ln]; got=[]; pool=list(cand.get(ln,[])); rng.shuffle(pool)
    for p in pool:
        if len(got)==k: break
        if p in used or norm(p) in humcat or not far(p): continue
        got.append(p); used.add(p)
    d=1
    while len(got)<k and d<10:
        for l2 in (ln-d,ln+d):
            pool=list(cand.get(l2,[])); rng.shuffle(pool)
            for p in pool:
                if len(got)==k: break
                if p in used or norm(p) in humcat or not far(p): continue
                got.append(p); used.add(p)
        d+=1
    out+=got
print(f'HeLa targets {len(hela)}, homolog-free Arabidopsis entrapments {len(out)}')
with open('hela_ath2.fasta','w') as f:
    for p in hela: f.write(f'>sp|{p}_target|{p}_target\n{p}\n')
    for p in out: f.write(f'>sp|{p}_p_target|{p}_ARATH_entrap\n{p}\n')
