"""Write the Sage-only TMTpro spectra (no ProSE candidate) of one file to a small mzML, with Sage's peptides."""
import json, sys
from compare import accepted, candidates, SAGE, PROSE, B
import pyopenms as p
d = sys.argv[1]; arm = "base"
sa, pa = accepted(SAGE / d), accepted(PROSE / d / arm)
pc = candidates(PROSE / d / arm / "native.pin", "score")
want = {s: sa[s] for s in sa if s not in pa and s not in pc}
exp = p.MSExperiment(); p.MzMLFile().load(str(B / f"prose-andes-reproduction/broad_benchmark/data/{d}.mzML"), exp)
keep = p.MSExperiment()
for s in exp:
    if s.getMSLevel() == 2 and int(s.getNativeID().split("scan=")[-1]) in want: keep.addSpectrum(s)
p.MzMLFile().store(f"probe/{d}_nocand.mzML", keep)
json.dump(want, open(f"probe/{d}_nocand.json", "w"))
print(keep.size(), "spectra written")
