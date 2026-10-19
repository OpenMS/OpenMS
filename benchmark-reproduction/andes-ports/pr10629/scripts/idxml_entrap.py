"""Print sequence and entrapment flag (all proteins ENTRAP_) of accepted (q <= 0.01) hits of a q-value idXML."""
import sys, xml.etree.ElementTree as ET
for path in sys.argv[1:]:
    root = ET.parse(path).getroot(); acc = {}
    for ph in root.iter("ProteinHit"): acc[ph.get("id")] = ph.get("accession")
    for pid in root.iter("PeptideIdentification"):
        h = pid.find("PeptideHit")
        if h is None or float(h.get("score")) > 0.01: continue
        refs = h.get("protein_refs", "").split()
        print(path, h.get("sequence"), all(acc[r].startswith("ENTRAP_") for r in refs) and bool(refs))
