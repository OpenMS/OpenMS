"""Independent count (plain XML parsing, no OpenMS) of accepted and entrapment PSMs in ProSE's in-process Percolator
idXML (main score = q-value), compared with EntrapmentEvaluator -te_ratio 1.  check_inproc.py <bin_dir> <idXML>..."""
from pathlib import Path
import os, subprocess, sys, tempfile, xml.etree.ElementTree as ET
bin_dir = sys.argv[1]; ok = True
with tempfile.NamedTemporaryFile(suffix=".tsv") as t:
    subprocess.run([f"{bin_dir}/EntrapmentEvaluator", "-in"] + sys.argv[2:] + ["-te_ratio", "1", "-out", t.name, "-no_progress"], check=True, capture_output=True,
                   env=dict(os.environ, LD_LIBRARY_PATH=str(Path(bin_dir).parent / "lib"), OPENMS_DATA_PATH=str(Path(bin_dir).parent / "share/OpenMS")))
    ee = {l.split("\t")[0]: (int(l.split("\t")[3]), int(l.split("\t")[4])) for l in open(t.name).read().splitlines()[1:] if l.split("\t")[1] == "all"}
for path in sys.argv[2:]:
    root = ET.parse(path).getroot()
    acc = {ph.get("id"): ph.get("accession") for ph in root.iter("ProteinHit")}
    n = ne = 0
    for pid in root.iter("PeptideIdentification"):
        assert "q-value" in pid.get("score_type")
        hits = pid.findall("PeptideHit")
        if not hits: continue
        h = min(hits, key=lambda x: float(x.get("score")))
        td = next((u.get("value") for u in h.findall("UserParam") if u.get("name") == "target_decoy"), "target")
        if td == "decoy" or float(h.get("score")) > 0.01: continue
        n += 1; refs = h.get("protein_refs", "").split()
        ne += bool(refs) and all(acc[r].startswith("ENTRAP_") for r in refs)
    key = path if len(sys.argv) > 3 else path.split("/")[-1]
    e = ee.get(key) or ee.get(path)
    ok &= (n, ne) == e
    print(path, "xml", (n, ne), "EntrapmentEvaluator", e)
print("identical:", ok)
