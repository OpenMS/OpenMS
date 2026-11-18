"""Pseudo-reversed decoy of every *digested* target peptide in a targets PQP (one entry per unique unmodified peptide),
for prediction with -generation:enzyme 'no cleavage'. Usage: make_decoy_fasta2.py targets.pqp out.fasta"""
import sys, sqlite3
pr=lambda p: p[:-1][::-1]+p[-1]
c=sqlite3.connect(f'file:{sys.argv[1]}?mode=ro',uri=True)
peps=sorted({r[0] for r in c.execute("select UNMODIFIED_SEQUENCE from PEPTIDE where DECOY=0")})
with open(sys.argv[2],'w') as f:
    for s in peps: f.write(f'>sp|DECOY_{s}|DECOY_{s}\n{pr(s)}\n')
print('target peptides', len(peps))
