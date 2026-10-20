"""Check that every ENTRAP_ protein has the composition and tryptic cleavage sites of its target, so that every
entrapment peptide is mass-matched to a target peptide. Usage: check_entrapment.py db.fasta"""
import sys, re
from collections import Counter
seqs = {}; head = None; buf = []
for line in open(sys.argv[1]):
    line = line.strip()
    if line.startswith('>'):
        if head: seqs[head] = ''.join(buf)
        head, buf = line[1:].split()[0], []
    elif line: buf.append(line)
if head: seqs[head] = ''.join(buf)
def sites(s): return [m.start() for m in re.finditer(r'[KR](?!P)', s)]
n = bad_comp = bad_sites = identical_pep = total_pep = 0
for acc, s in seqs.items():
    if not acc.startswith('ENTRAP_'): continue
    t = seqs[acc[len('ENTRAP_'):]]; n += 1
    bad_comp += Counter(s) != Counter(t)
    if sites(s) != sites(t): bad_sites += 1; continue
    cut = [0] + [i + 1 for i in sites(t)] + [len(t)]
    for a, b in zip(cut, cut[1:]):
        if b - a >= 7:
            total_pep += 1; identical_pep += s[a:b] == t[a:b]
            assert Counter(s[a:b]) == Counter(t[a:b])
print(f"entrapment proteins {n}; composition differs {bad_comp}; cleavage sites differ {bad_sites}; "
      f"peptides >=7 aa {total_pep}, all composition-matched, identical to target {identical_pep}")
