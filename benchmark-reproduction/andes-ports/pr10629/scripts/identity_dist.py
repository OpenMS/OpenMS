"""Identity of every entrapment peptide (Trypsin/P, 0-2 missed cleavages, 7-40) to the target peptide at the same position.
identity_dist.py <db.fasta>"""
import statistics, sys
seqs, h = {}, None
for l in open(sys.argv[1]):
    l = l.rstrip()
    if l.startswith(">"): h = l[1:].split()[0]; seqs[h] = []
    else: seqs[h].append(l)
seqs = {k: "".join(v) for k, v in seqs.items()}
ident = {}
for acc, s in seqs.items():
    if not acc.startswith("ENTRAP_"): continue
    t = seqs.get(acc[7:])
    if t is None or len(t) != len(s): continue
    st = [0] + [i + 1 for i in range(len(s) - 1) if s[i] in "KR"] + [len(s)]
    for i in range(len(st) - 1):
        for j in range(i + 1, min(i + 4, len(st))):
            if 7 <= st[j] - st[i] <= 40:
                e = s[st[i]:st[j]]
                ident.setdefault(e, sum(a == b for a, b in zip(e, t[st[i]:st[j]])) / len(e))
v = list(ident.values())
print(f"{sys.argv[1].split('/')[-1]}: entrapment peptides {len(v)}  mean identity {statistics.mean(v):.3f}  >=0.6: {sum(x >= 0.6 for x in v) / len(v):.1%}  ==1: {sum(x == 1 for x in v) / len(v):.2%}")
