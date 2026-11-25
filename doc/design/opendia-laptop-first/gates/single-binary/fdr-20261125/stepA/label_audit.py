#!/usr/bin/env python3
"""Step A: entrapment label audit. The library build digests each entrapment FASTA entry again, so a foreign entry with
missed cleavages can yield sub-peptides that are HUMAN sequences. Such a precursor carries only entrapment protein
accessions, so evaluators count it as an entrapment although it is a human peptide that can be present.
Lists every library target peptide whose label is entrapment and whose sequence (I=L) is an exact substring of a protein in
the human FASTA used at build time. Output: seq per line (for the evaluator) and a summary. Diagnosis only.
Usage: label_audit.py <lib.pqp> <human.fasta> <out.tsv>"""
import sys, sqlite3, collections
lib, fa, out = sys.argv[1:4]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
ent = {}
for seq, e in c.execute("""select pe.UNMODIFIED_SEQUENCE, min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end)
  from PEPTIDE pe join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID
  where pe.DECOY=0 group by pe.ID"""):
    ent[seq] = min(e, ent.get(seq, 1))
E = sorted(s for s, e in ent.items() if e == 1)
print(f'library target peptides {len(ent)}, entrapment-labelled {len(E)}', flush=True)
prots = []; s = []
for l in open(fa):
    if l.startswith('>'):
        if s: prots.append(''.join(s)); s = []
    else: s.append(l.strip())
if s: prots.append(''.join(s))
prots = [p.replace('I', 'L') for p in prots]
bylen = collections.defaultdict(list)
for p in E: bylen[len(p)].append(p)
hit = []
for L in sorted(bylen):
    win = set()
    for p in prots:
        for i in range(len(p) - L + 1): win.add(p[i:i + L])
    h = [p for p in bylen[L] if p.replace('I', 'L') in win]
    hit += h
    print(f'  length {L}: entrapment peptides {len(bylen[L])}, human substrings {len(h)}', flush=True)
    del win
print(f'entrapment-labelled peptides that are exact human substrings (I=L): {len(hit)} of {len(E)} ({len(hit) / len(E):.4f})', flush=True)
open(out, 'w').write('seq\n' + '\n'.join(hit) + '\n')
