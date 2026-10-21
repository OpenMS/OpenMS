"""Add benchmark class meta values to every PeptideHit of an idXML (for FalseDiscoveryRate -PSM_class <meta value>).
Streams the file line by line (whole runs do not fit in memory as a tree); expects one PeptideHit start tag per line,
as OpenMS writes idXML.

bench_len  peptide length
bench_mc   missed cleavages (Trypsin/P: internal K/R), "0" / "1" / "2"
bench_z    precursor charge, "2" / "3" / "4+" (1 stays "1")
bench_term "specific" (both termini tryptic in any protein; N-terminal Met clipping counts as specific) / "semi"
bench_iso  "canonical" if any protein of the hit is not a UniProt isoform (acc-N), else "isoform"
           (DECOY_ / ENTRAP_ prefixes stripped first, so decoys and entrapment carry their target's class)
usage: annotate.py in.idXML out.idXML
"""
import re, sys

ISO = re.compile(r"^(?:sp|tr)\|[^|]+-\d+\|")
ATTR = re.compile(r'(\w+)="([^"]*)"')
def strip(acc):
    for p in ("DECOY_", "ENTRAP_"):
        if acc.startswith(p): acc = acc[len(p):]
    return acc

acc = {}; n = 0
with open(sys.argv[1]) as fin, open(sys.argv[2], "w") as fout:
    for line in fin:
        s = line.lstrip()
        if s.startswith("<ProteinHit "):
            a = dict(ATTR.findall(s)); acc[a["id"]] = a["accession"]
        elif s.startswith("<PeptideHit "):
            a = dict(ATTR.findall(s))
            seq = re.sub(r"\([^)]*\)|\[[^\]]*\]", "", a["sequence"]).replace(".", "")
            mc = sum(1 for aa in seq[:-1] if aa in "KR")
            z = int(a["charge"])
            accs = [strip(acc[r]) for r in a.get("protein_refs", "").split()]
            iso = "canonical" if any(not ISO.match(x) for x in accs) or not accs else "isoform"
            before = a.get("aa_before", "").split(); after = a.get("aa_after", "").split()
            starts = a.get("start", "").split()
            ends = list(zip(before, after, starts + [""] * len(before))) or [("", "", "")]
            spec = any((b in "KR[-" or b == "" or (b == "M" and st == "1")) and (seq[-1] in "KR" or e in "]-")
                       for b, e, st in ends)
            ind = line[:len(line) - len(s)] + "\t"
            ups = "".join(f'{ind}<UserParam type="string" name="{k}" value="{v}"/>\n' for k, v in (
                ("bench_len", str(len(seq))), ("bench_term", "specific" if spec else "semi"),
                ("bench_mc", str(min(mc, 2))), ("bench_z", str(z) if z < 4 else "4+"), ("bench_iso", iso)))
            n += 1
            if s.rstrip().endswith("/>"):
                line = line.rstrip()[:-2].rstrip() + " >\n" + ups + line[:len(line) - len(s)] + "</PeptideHit>\n"
            else:
                line = line + ups
        fout.write(line)
print("annotated", n, "hits")
