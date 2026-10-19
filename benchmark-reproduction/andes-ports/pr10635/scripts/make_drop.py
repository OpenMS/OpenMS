#!/usr/bin/env python3
"""'drop' arm (the andes #102 behaviour): remove the PSMs the isotope_mimic=true PIN flags from the
isotope_mimic=false PIN, then rescore with Percolator (seeds 1/42/137). The spectrum's next candidate competes instead.

usage: make_drop.py <runs/pr/<db>> <dataset> <percolator>"""
import csv, subprocess, sys
from pathlib import Path

base, ds, perc = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
with open(base / "on" / ds / "search.pin") as f:
    r = csv.reader(f, delimiter="\t")
    h = next(r)
    i_s, i_p, i_f = h.index("SpecId"), h.index("Peptide"), h.index("isotope_mimic")
    flagged = {(row[i_s], row[i_p]) for row in r if row[i_f] == "1"}
out = base / "drop" / ds
out.mkdir(parents=True, exist_ok=True)
kept = dropped = 0
with open(base / "off" / ds / "search.pin") as f, open(out / "search.pin", "w") as g:
    h = f.readline()
    g.write(h)
    cols = h.rstrip("\n").split("\t")
    i_s, i_p = cols.index("SpecId"), cols.index("Peptide")
    for line in f:
        row = line.rstrip("\n").split("\t")
        if (row[i_s], row[i_p]) in flagged:
            dropped += 1
        else:
            kept += 1
            g.write(line)
(out / "drop.txt").write_text(f"kept {kept} dropped {dropped}\n")
for seed in ["1", "42", "137"]:
    d = out / f"s{seed}"
    d.mkdir(exist_ok=True)
    with open(d / "percolator.log", "w") as log:
        subprocess.run([perc, "-Y", "-U", "--seed", seed, "--results-psms", str(d / "target.tsv"),
                        "--decoy-results-psms", str(d / "decoy.tsv"), str(out / "search.pin")],
                       check=True, stdout=log, stderr=subprocess.STDOUT)
print(ds, "kept", kept, "dropped", dropped)
