"""Foreign-proteome entrapment (Wen et al. 2025): targets of a target-decoy FASTA + a foreign proteome as ENTRAP_ targets,
then full-reversal DECOY_ decoys of both (#10364 style).  make_foreign.py <target_decoy.fasta> <foreign.fasta> <out.fasta>"""
import hashlib, sys
def read(p):
    out, h, s = [], None, []
    for l in open(p):
        if l.startswith(">"):
            if h is not None: out.append((h, "".join(s)))
            h, s = l[1:].rstrip("\n"), []
        else: s.append(l.strip())
    out.append((h, "".join(s)))
    return out
targets = [(h, s) for h, s in read(sys.argv[1]) if not h.startswith("DECOY_")]
foreign = [("ENTRAP_" + h, s) for h, s in read(sys.argv[2])]
with open(sys.argv[3], "w") as f:
    for h, s in targets + foreign: f.write(f">{h}\n{s}\n")
    for h, s in targets + foreign: f.write(f">DECOY_{h.split()[0]}\n{s[::-1]}\n")
print(sys.argv[3], len(targets), len(foreign), hashlib.sha256(open(sys.argv[3], "rb").read()).hexdigest())
