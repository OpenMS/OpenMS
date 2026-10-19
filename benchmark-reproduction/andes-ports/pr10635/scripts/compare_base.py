#!/usr/bin/env python3
"""Option-off isolation check: PIN, Percolator results and idXML of the PR build (annotate:isotope_mimic=false) against
the merge-base build, same inputs and parameters. idXML is compared without the run date and search_engine_version
(which names the git branch).

usage: compare_base.py <runs dir> <out tsv>"""
import hashlib, re, sys
from pathlib import Path

runs, out = Path(sys.argv[1]), Path(sys.argv[2])
DATE = re.compile(rb'search_engine_version="[^"]*"|date="[^"]*"|creation_date="[^"]*"|<\?xml-stylesheet[^>]*>')

def h(p, norm=False):
    b = Path(p).read_bytes()
    if norm:
        b = DATE.sub(b"", b)
    return hashlib.sha256(b).hexdigest()

lines = ["dataset\tpin_sha256_base\tpin_identical\tidxml_identical_normalized\tpercolator_s1_s42_s137_identical"]
for d in sorted((runs / "base" / "plain" / "off").iterdir()):
    pr = runs / "pr" / "plain" / "off" / d.name
    pin = h(d / "search.pin") == h(pr / "search.pin")
    idx = h(d / "search.idXML", True) == h(pr / "search.idXML", True)
    perc = all(h(d / f"s{s}" / "target.tsv") == h(pr / f"s{s}" / "target.tsv") for s in ("1", "42", "137"))
    lines.append(f"{d.name}\t{h(d / 'search.pin')}\t{pin}\t{idx}\t{perc}")
out.write_text("\n".join(lines) + "\n")
print(out.read_text())
