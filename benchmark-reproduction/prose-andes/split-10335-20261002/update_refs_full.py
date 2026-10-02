"""Apply the change between two builds' replayed ProSE TOPP outputs to the TOPP references.

The replay of <base> must equal the reference line for line, except for lines the harness cannot
reproduce (search_engine_version, spectra_data path, db, date). The line-level diff from the <base>
replay to the <build> replay is then applied to the reference, and every changed line must contain
one of the allowed substrings (or be the extra_features line).
usage: update_refs_full.py <base build> <new build> <openms checkout to edit> <allowed substring>...
"""
import re, sys, difflib
from pathlib import Path
ROOT = Path(__file__).resolve().parent
base, build, checkout, allowed = sys.argv[1], sys.argv[2], Path(sys.argv[3]), sys.argv[4:]
harness_only = re.compile(r'search_engine_version=|name="spectra_data"| db="| date="')
for test in ["ProSE_1", "ProSE_2", "ProSE_4", "ProSE_10"]:
    ref_path = checkout / "src/tests/topp" / f"{test}_out.idXML"
    ref = ref_path.read_text().split("\n")
    nodate = lambda s: re.sub(r' date="[^"]*"', ' date=""', s)  # replay run times differ
    old = nodate((ROOT / "topp-emulation" / base / test / "out.idXML").read_text()).split("\n")
    new = nodate((ROOT / "topp-emulation" / build / test / "out.idXML").read_text()).split("\n")
    assert len(ref) == len(old), test
    mismatch = [i for i, (r, o) in enumerate(zip(ref, old)) if r != o]
    assert all(harness_only.search(ref[i]) for i in mismatch), (test, [ref[i] for i in mismatch if not harness_only.search(ref[i])])
    out, changed = [], 0
    for tag, i1, i2, j1, j2 in difflib.SequenceMatcher(None, old, new, autojunk=False).get_opcodes():
        if tag == "equal":
            out.extend(ref[i1:i2])
            continue
        for l in old[i1:i2] + new[j1:j2]:
            assert any(a in l for a in allowed) or "extra_features" in l, (test, tag, l)
            assert not harness_only.search(l), (test, l)
        out.extend(new[j1:j2])
        changed += (i2 - i1) + (j2 - j1)
    ref_path.write_text("\n".join(out))
    print(f"{test}: {changed} changed lines")
