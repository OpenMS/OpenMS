# T/E as andes docs/benchmarks/reproduce/build_databases.sh (bigbio/andes fa0195a) computes it, verbatim. Input: target + entrapment FASTA without decoys.
import sys, re
tgt = ent = 0
hdr = None; seq = []
def digest(s):
    sites = [0] + [i+1 for i in range(len(s)-1) if s[i] in "KR"] + [len(s)]
    sites = sorted(set(sites)); out = set()
    for i in range(len(sites)-1):
        for j in range(i+1, min(i+3, len(sites))):  # andes default: <= 1 missed cleavage
            p = s[sites[i]:sites[j]]
            if 6 <= len(p) <= 50: out.add(p)
    return out
T, E = set(), set()
for line in open(sys.argv[1]):
    if line.startswith(">"):
        if hdr is not None:
            (E if hdr.startswith(">ENTRAP_") else T).update(digest("".join(seq)))
        hdr = line; seq = []
    else: seq.append(line.strip())
if hdr is not None:
    (E if hdr.startswith(">ENTRAP_") else T).update(digest("".join(seq)))
E -= T  # a peptide also in a target protein is a target identification
r = len(T)/max(len(E),1)
print(f"    target peptides {len(T)}  entrapment-only peptides {len(E)}")
print(f"    T/E = {r:.3f}  ->  FDP = (hits/total) x {1+r:.3f}")
print(f"    NOTE: this factor is a property of YOUR build. Never assume 1:1 (factor 2).")
