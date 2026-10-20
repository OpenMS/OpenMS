"""Targets + DECOY_ full-sequence reversals, byte-identical to #10364's prepare_databases.py
(target: '>id description', one-line sequence; decoy: '>DECOY_id', reversed sequence)."""
import sys
def read(p):
    h, s = None, []
    for line in open(p):
        line = line.rstrip("\r\n")
        if line.startswith(">"):
            if h is not None: yield h, "".join(s)
            h, s = line[1:], []
        elif line: s.append(line.strip())
    if h is not None: yield h, "".join(s)
ents = []
for h, s in read(sys.argv[1]):
    ident, _, desc = h.partition(" ")
    ents.append((ident, desc, s))
with open(sys.argv[2], "w") as o:
    for i, d, s in ents: o.write(">" + i + " " + d + "\n" + s + "\n")
    for i, d, s in ents: o.write(">DECOY_" + i + "\n" + s[::-1] + "\n")
