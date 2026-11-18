"""HeLa entrapment FASTA (FDRBench-style peptide pairs, label-free construction from the library only).
Each unique unmodified target peptide of the HeLa library is one entry `<pep>_target`, followed by `<shuf>_p_target`:
the residues shuffled with every K and R (cleavage sites) and the C-terminal residue kept in place, seeded per peptide,
rejecting shuffles equal to any target peptide or an earlier shuffle (50 tries; else the pair is dropped). The digest of
an entrapment therefore has the same cleavage pattern and length as its parent."""
import sqlite3, random, hashlib
L='/home/singjc/Documents/github/testing_auto_irt/diaPASEF_evosep_published/200623_hela_evosep_fractions_library_ptypic_decoy.pqp'
c=sqlite3.connect(f'file:{L}?mode=ro',uri=True)
peps=sorted({r[0] for r in c.execute("select pe.UNMODIFIED_SEQUENCE from PEPTIDE pe where pe.DECOY=0")})
tset=set(peps); used=set(); out=[]; dropped=0
for p in peps:
    rng=random.Random(int(hashlib.sha256(p.encode()).hexdigest()[:16],16))
    free=[i for i in range(len(p)-1) if p[i] not in 'KR']
    s=None
    for _ in range(50):
        vals=[p[i] for i in free]; rng.shuffle(vals); q=list(p)
        for i,v in zip(free,vals): q[i]=v
        cand=''.join(q)
        if cand not in tset and cand not in used: s=cand; break
    if s is None: dropped+=1; continue
    used.add(s); out.append((p,s))
with open('hela_entrap.fasta','w') as f:
    for p,s in out: f.write(f'>sp|{p}_target|{p}_target\n{p}\n>sp|{s}_p_target|{s}_p_target\n{s}\n')
print(f'target peptides {len(peps)}, pairs {len(out)}, dropped {dropped}')
