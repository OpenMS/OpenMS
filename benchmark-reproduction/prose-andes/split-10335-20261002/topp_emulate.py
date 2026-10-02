"""Replay the ProSE TOPP test searches through the harness driver and diff two builds.

The wheel's libOpenMS has no TOPPBase, so the ProSE TOPP tool cannot be linked here. The
TOPP tests' Search:* parameters (INI plus command-line overrides) are fed to
ProSEAlgorithm::search() of each harness build instead; the resulting idXML files show the
exact search-parameter and hit differences a source change causes in the TOPP outputs.

usage: topp_emulate.py <build_a> <build_b> --openms <checkout providing src/tests/topp>
"""

import argparse
import difflib
import json
import os
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BUNDLE = ROOT.parent / "prose-andes-reproduction"
ENV = dict(
    os.environ,
    OMP_NUM_THREADS="1",
    OPENMS_DATA_PATH=str(BUNDLE / "nightly/pyopenms/share/OpenMS"),
    PYTHONPATH=str(BUNDLE / "nightly") + os.pathsep + str(BUNDLE / "deps"),
)

# test: (ini or None, mzML, fasta, Search:* overrides)
TESTS = {
    "ProSE_1": ("ProSE_1.ini", "SimpleSearchEngine_1.mzML", "SimpleSearchEngine_1.fasta", {}),
    "ProSE_2": ("ProSE_2.ini", "SimpleSearchEngine_1.mzML", "SimpleSearchEngine_1.fasta", {}),
    "ProSE_4": ("ProSE_4.ini", "SimpleSearchEngine_1_shifted_7ppm.mzML", "SimpleSearchEngine_1.fasta", {}),
    "ProSE_5_chunked": ("ProSE_1.ini", "SimpleSearchEngine_1.mzML", "SimpleSearchEngine_1.fasta", {"database:chunk_size": 5}),
    "ProSE_6_fdr": ("ProSE_6_fdr.ini", "THIRDPARTY/DatabaseSuitability_in_spec.mzML",
                    "THIRDPARTY/DatabaseSuitability_database.fasta", {}),
    "ProSE_6_fdr_chunked": ("ProSE_6_fdr.ini", "THIRDPARTY/DatabaseSuitability_in_spec.mzML",
                            "THIRDPARTY/DatabaseSuitability_database.fasta", {"database:chunk_size": 20}),
    "ProSE_8_auto": ("ProSE_6_fdr.ini", "THIRDPARTY/DatabaseSuitability_in_spec.mzML",
                     "THIRDPARTY/DatabaseSuitability_database.fasta", {"decoys": "auto"}),
    "ProSE_10": (None, "ProSE_10_duplicate_peaks.mzML", "ProSE_10.fasta",
                 {"fragment:mass_tolerance": 10.0, "fragment:mass_tolerance_unit": "ppm"}),
}

READ_INI = """
import json, sys
import pyopenms as p
param = p.Param()
p.ParamXMLFile().load(sys.argv[1], param)
out = {}
for key in param.keys():
    k = key.decode() if isinstance(key, bytes) else key
    if not k.startswith("ProSE:1:Search:"):
        continue
    v = param.getValue(key)
    if isinstance(v, bytes):
        v = v.decode()
    if isinstance(v, bool):  # INI items of type bool; ProSE takes them as "true"/"false" strings
        v = "true" if v else "false"
    if isinstance(v, list):
        v = [x.decode() if isinstance(x, bytes) else x for x in v]
    out[k[len("ProSE:1:Search:"):]] = v
print(json.dumps(out))
"""


def params_for(topp, ini, overrides):
    params = {}
    if ini:
        params = json.loads(subprocess.check_output(["python3", "-c", READ_INI, str(topp / ini)], env=ENV, text=True))
    params.update(overrides)
    return params


def run(build, name, topp, out_root):
    ini, mzml, fasta, overrides = TESTS[name]
    out = out_root / build / name
    out.mkdir(parents=True, exist_ok=True)
    params = params_for(topp, ini, overrides)
    # Keys unknown to an older build only trigger a warning; keep the replay identical.
    (out / "params.json").write_text(json.dumps(params, indent=2))
    with (out / "search.log").open("w") as log:
        rc = subprocess.run([str(ROOT / "builds" / build / "prose"), str(topp / mzml), str(topp / fasta),
                             str(out / "params.json"), str(out / "hits.tsv"), str(out / "out.idXML")],
                            env=ENV, stdout=log, stderr=log).returncode
    return rc, out


def normalized(path):
    text = path.read_text().splitlines()
    return [re.sub(r'date="[^"]*"', 'date=""', re.sub(r'db="[^"]*"', 'db=""', line)) for line in text]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--openms", type=Path, required=True)
    ap.add_argument("--tests", nargs="*", default=list(TESTS))
    ap.add_argument("--out", type=Path, default=ROOT / "topp-emulation")
    args = ap.parse_args()
    topp = args.openms / "src/tests/topp"
    for name in args.tests:
        ra, oa = run(args.a, name, topp, args.out)
        rb, ob = run(args.b, name, topp, args.out)
        print(f"== {name}: rc {args.a}={ra} {args.b}={rb}")
        if ra or rb:
            continue
        diff = list(difflib.unified_diff(normalized(oa / "out.idXML"), normalized(ob / "out.idXML"),
                                         f"{args.a}/{name}", f"{args.b}/{name}", n=1, lineterm=""))
        print("\n".join(diff) if diff else "   identical idXML (dates/db path ignored)")


if __name__ == "__main__":
    main()
