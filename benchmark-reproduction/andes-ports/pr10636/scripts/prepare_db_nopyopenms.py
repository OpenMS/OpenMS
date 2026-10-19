"""Pure-Python twin of broad_benchmark/prepare_databases.py (FASTA part only); checks the frozen SHA256."""
import hashlib, json, sys
from pathlib import Path
R = Path(sys.argv[1])  # broad_benchmark dir
W = R.parent
def load(p):
    out, ident, desc, seq = [], None, "", []
    for line in open(p):
        line = line.rstrip("\n").rstrip("\r")
        if line.startswith(">"):
            if ident is not None: out.append((ident, desc, "".join(seq)))
            h = line[1:].strip(); parts = h.split(None, 1)
            ident = parts[0]; desc = parts[1] if len(parts) > 1 else ""; seq = []
        elif line.strip():
            seq.append(line.strip())
    if ident is not None: out.append((ident, desc, "".join(seq)))
    return out
src = R / "data/hye_targets.fasta"
allE = load(src)
expected = json.load(open(W / "databases.json"))
for db in ["hy", "hye", "human"]:
    e = allE
    if db == "hy": e = [x for x in e if x[0].endswith(("_HUMAN", "_YEAST")) or "Cont_" in x[0]]
    if db == "human": e = [x for x in e if x[0].endswith("_HUMAN") or "Cont_" in x[0]]
    out = R / "data" / (db + "_td.fasta")
    with out.open("w") as f:
        for i, d, s in e: f.write(">" + i + " " + d + "\n" + s + "\n")
        for i, d, s in e: f.write(">DECOY_" + i + "\n" + s[::-1] + "\n")
    h = hashlib.sha256(out.read_bytes()).hexdigest()
    print(db, len(e), h, "OK" if h == expected[db]["sha256"] else "MISMATCH " + expected[db]["sha256"])
