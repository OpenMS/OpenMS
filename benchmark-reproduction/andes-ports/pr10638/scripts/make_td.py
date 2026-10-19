"""hye_td.fasta of the #10364 suite (prepare_databases.py without pyOpenMS): targets, then DECOY_ full reversals."""
import hashlib, sys
src, dst = sys.argv[1], sys.argv[2]
entries, head, seq = [], None, []
for line in open(src):
    line = line.rstrip("\n\r")
    if line.startswith(">"):
        if head is not None: entries.append((head, "".join(seq)))
        head, seq = line[1:], []
    elif line.strip(): seq.append(line.strip())
if head is not None: entries.append((head, "".join(seq)))
with open(dst, "w") as f:
    for h, s in entries:
        ident, _, desc = h.partition(" ")
        f.write(f">{ident} {desc}\n{s}\n")
    for h, s in entries: f.write(f">DECOY_{h.split(' ')[0]}\n{s[::-1]}\n")
print(len(entries), hashlib.sha256(open(dst, "rb").read()).hexdigest())
