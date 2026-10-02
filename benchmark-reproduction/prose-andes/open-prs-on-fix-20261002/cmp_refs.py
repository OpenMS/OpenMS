"""Compare a build's replayed TOPP outputs with an OpenMS checkout's ProSE idXML references.

Ignores db paths, dates and version strings. usage: cmp_refs.py <build> <checkout>
"""
import difflib, re, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent
build, checkout = sys.argv[1], Path(sys.argv[2])
def norm(s):
    s = re.sub(r' db="[^"]*"', '', s)
    s = re.sub(r'date="[^"]*"', '', s)
    s = re.sub(r'search_engine_version="[^"]*"', '', s)
    s = re.sub(r'name="spectra_data" value="[^"]*"', 'name="spectra_data"', s)
    s = re.sub(r'<\?xml-stylesheet[^>]*>\n', '', s)
    s = re.sub(r'<UserParam type="string" name="(db_path|database)"[^>]*/>\n', '', s)
    return s
bad = 0
for ref_path in sorted((checkout / "src/tests/topp").glob("ProSE_*_out.idXML")):
    test = ref_path.name[:-len("_out.idXML")]
    out = ROOT / "topp-emulation" / build / test / "out.idXML"
    if not out.exists():
        print(test, "no replay"); continue
    a, b = norm(ref_path.read_text()).splitlines(), norm(out.read_text()).splitlines()
    d = [l for l in difflib.unified_diff(a, b, "ref", "replay", lineterm="", n=0) if not l.startswith(("---", "+++", "@@"))]
    print(test, "MATCH" if not d else f"{len(d)} differing lines")
    if d:
        bad += 1
        print("\n".join(d[:20]))
sys.exit(1 if bad else 0)
