"""For accepted entrapment PSMs: sequence identity to the target peptide at the same position (the pair) and whether that
target peptide is itself accepted in the same file and seed. Positional pairing needs entrapment proteins with the target's
cleavage sites (all designs here except the whole-protein shuffle).  twin_check.py <db.fasta> <run_dir>..."""
import statistics, sys
db, runs = sys.argv[1], sys.argv[2:]
seqs, h = {}, None
for l in open(db):
    l = l.rstrip()
    if l.startswith(">"): h = l[1:].split()[0]; seqs[h] = []
    else: seqs[h].append(l)
seqs = {k: "".join(v) for k, v in seqs.items()}
pair = {}
for acc, s in seqs.items():
    if not acc.startswith("ENTRAP_"): continue
    t = seqs.get(acc[7:])
    if t is None or len(t) != len(s): continue
    st = [0] + [i + 1 for i in range(len(s) - 1) if s[i] in "KR"] + [len(s)]
    for i in range(len(st) - 1):
        for j in range(i + 1, min(i + 4, len(st))):
            if 7 <= st[j] - st[i] <= 40: pair.setdefault(s[st[i]:st[j]], t[st[i]:st[j]])
ident, twin, n = [], 0, 0
for r in runs:
    for sd in (1, 42, 137):
        rows = [l.split("\t") for l in open(f"{r}/s{sd}/target.tsv").read().splitlines()[1:]]
        acc = [f for f in rows if float(f[2]) <= 0.01]
        accepted_peps = {f[4].split(".")[1] for f in acc}
        for f in acc:
            prots = [p for p in f[5:] if p and not p.startswith("DECOY_")]
            if not prots or not all(p.startswith("ENTRAP_") for p in prots): continue
            pep = "".join(c for c in f[4].split(".")[1] if c.isalpha() and c.isupper())
            t = pair.get(pep)
            if t is None: continue
            n += 1
            ident.append(sum(a == b for a, b in zip(pep, t)) / len(pep))
            twin += any("".join(c for c in x if c.isalpha() and c.isupper()) == t for x in accepted_peps)
print(f"paired entrapment PSMs {n}  mean identity to pair {statistics.mean(ident):.3f}  identity>=0.6: {sum(x >= 0.6 for x in ident) / n:.1%}  pair target also accepted: {twin / n:.1%}")
