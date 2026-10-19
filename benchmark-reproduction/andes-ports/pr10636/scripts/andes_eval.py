"""andes reference: per class (canonical / isoform, as andes assigns it: canonical if any protein is not an isoform),
target and decoy PSMs at Percolator q <= 0.01 (global; best PSM per scan), and the PSMs andes finally accepts after its
class FDR (<run>.q0p01.tsv). usage: andes_eval.py <run prefix>"""
import csv, re, sys
from collections import Counter
ISO = re.compile(r"^(?:DECOY_)?(?:sp|tr)\|[^|]+-\d+\|")
pre = sys.argv[1]
def cls(prots): return "canonical" if any(not ISO.match(p) for p in prots) else "isoform"
best = {}
for kind in ("target", "decoy"):
    with open(f"{pre}.percolator.{kind}.psms.txt") as f:
        next(f)
        for line in f:
            x = line.rstrip("\n").split("\t")
            scan = x[0].split("_")[0]
            s = float(x[1])
            if scan not in best or s > best[scan][0]:
                best[scan] = (s, float(x[2]), kind, cls(x[5:]))
c = Counter((k, c) for s, q, k, c in best.values() if q <= 0.01)
for cl in ("canonical", "isoform"):
    t, d = c[("target", cl)], c[("decoy", cl)]
    print(f"global q<=0.01 {cl}: targets {t} decoys {d} decoy/target {100*d/max(t,1):.2f}%")
rows = list(csv.DictReader(open(f"{pre}.q0p01.tsv"), delimiter="\t"))
pc = "Protein" if "Protein" in rows[0] else [k for k in rows[0] if "rotein" in k][0]
fin = Counter(cls(r[pc].split(";")) for r in rows)
print("andes final (class FDR):", dict(fin), "total", len(rows))
