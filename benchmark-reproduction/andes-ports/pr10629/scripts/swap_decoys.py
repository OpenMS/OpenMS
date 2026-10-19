"""Replace the decoys of a target + entrapment database by full-sequence reversals with prefix DECOY_ (the decoys of the
#10364 databases and of make_entrapment.py): targets and entrapment first, then their decoys, in input order.
usage: swap_decoys.py in.fasta out.fasta"""
import hashlib, sys
entries, h, seq = [], None, []
for line in open(sys.argv[1]):
    if line.startswith(">"):
        if h is not None: entries.append((h, "".join(seq)))
        h, seq = line[1:].rstrip("\n"), []
    else: seq.append(line.strip())
entries.append((h, "".join(seq)))
keep = [(h, s) for h, s in entries if not h.startswith("DECOY_")]
with open(sys.argv[2], "w") as f:
    for h, s in keep: f.write(f">{h}\n{s}\n")
    for h, s in keep: f.write(f">DECOY_{h.split()[0]}\n{s[::-1]}\n")
print(sys.argv[2], len(keep), hashlib.sha256(open(sys.argv[2], "rb").read()).hexdigest())
