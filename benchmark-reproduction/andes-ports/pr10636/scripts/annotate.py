"""Add benchmark class meta values to every PeptideHit of an idXML (for FalseDiscoveryRate -PSM_class <meta value>).

bench_len  peptide length
bench_mc   missed cleavages (Trypsin/P: internal K/R), "0" / "1" / "2"
bench_z    precursor charge, "2" / "3" / "4+"
bench_term "specific" (both termini tryptic in any protein) / "semi"
bench_iso  "canonical" if any protein of the hit is not a UniProt isoform (acc-N), else "isoform"
           (DECOY_ / ENTRAP_ prefixes stripped first, so decoys and entrapment carry their target's class)
usage: annotate.py in.idXML out.idXML
"""
import re, sys
from lxml import etree as ET

ISO = re.compile(r"^(?:sp|tr)\|[^|]+-\d+\|")
def strip(acc):
    for p in ("DECOY_", "ENTRAP_"):
        if acc.startswith(p): acc = acc[len(p):]
    return acc

tree = ET.parse(sys.argv[1], ET.XMLParser(huge_tree=True, remove_blank_text=False))
root = tree.getroot()
ns = root.nsmap.get(None)
q = (lambda t: "{%s}%s" % (ns, t)) if ns else (lambda t: t)
acc = {ph.get("id"): ph.get("accession") for ph in root.iter(q("ProteinHit"))}
n = 0
for hit in root.iter(q("PeptideHit")):
    seq = re.sub(r"\([^)]*\)|\[[^\]]*\]|^\.|\.$", "", hit.get("sequence"))
    seq = seq.replace(".", "")
    mc = sum(1 for aa in seq[:-1] if aa in "KR")
    z = int(hit.get("charge"))
    accs = [strip(acc[r]) for r in (hit.get("protein_refs") or "").split()]
    iso = "canonical" if any(not ISO.match(a) for a in accs) or not accs else "isoform"
    before = (hit.get("aa_before") or "").split(); after = (hit.get("aa_after") or "").split()
    starts = (hit.get("start") or "").split()
    ends = list(zip(before, after, starts + [""] * len(before))) or [("", "", "")]
    # N-terminal Met clipping (start 1 after M) counts as specific
    spec = any((b in "KR[-" or b == "" or (b == "M" and st == "1")) and (seq[-1] in "KR" or a in "]-") for b, a, st in ends)
    for name, val in (("bench_len", str(len(seq))), ("bench_term", "specific" if spec else "semi"), ("bench_mc", str(min(mc, 2))), ("bench_z", str(z) if z < 4 else "4+"), ("bench_iso", iso)):
        up = ET.SubElement(hit, q("UserParam"), type="string", name=name, value=val)
    n += 1
tree.write(sys.argv[2], xml_declaration=True, encoding="UTF-8")
print("annotated", n, "hits")
