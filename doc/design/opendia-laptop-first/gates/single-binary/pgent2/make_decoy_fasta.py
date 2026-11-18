"""Pseudo-reversed decoy sequences (all but the C-terminal residue reversed) of every unmodified peptide in a target FASTA
of single-peptide entries; header DECOY_<peptide>. Usage: make_decoy_fasta.py in.fasta out.fasta"""
import sys
pr=lambda p: p[:-1][::-1]+p[-1]
seqs=[l.strip() for l in open(sys.argv[1]) if not l.startswith('>')]
with open(sys.argv[2],'w') as f:
    for s in sorted(set(seqs)): f.write(f'>sp|DECOY_{s}|DECOY_{s}\n{pr(s)}\n')
