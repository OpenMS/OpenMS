"""Share of matchable b/y ions that ProSE's deisotoping removes, per file.

For up to N spectra Sage accepts (majority of seeds), every b/y ion of Sage's peptide (charge 1, and 2 for precursor
charge >= 3) that has a raw peak within 20 ppm counts as matchable. After Deisotoper::deisotopeAndSingleCharge with
ProSE's arguments (and variants), an ion is kept if a peak remains within 20 ppm of its m/z."""
import json, random, sys, collections
from compare import accepted, SAGE, MANIFEST, B
import pyopenms as p
import csv
N = int(sys.argv[1]); files = sys.argv[2:]
VARIANTS = {
    "prose (start_check=2)": dict(use_decreasing_model=True, start_intensity_check=2),
    "start_check=1": dict(use_decreasing_model=True, start_intensity_check=1),
    "no decreasing model": dict(use_decreasing_model=False, start_intensity_check=2),
}
def ions(seq, z):
    out = []
    for i in range(1, seq.size()):
        for kind, frag, rt in (("b", seq.getPrefix(i), p.Residue.ResidueType.BIon), ("y", seq.getSuffix(i), p.Residue.ResidueType.YIon)):
            for c in (1, 2):
                if c == 1 or z >= 3: out.append(frag.getMonoWeight(rt, c) / c)
    return out
def present(mzs, m):
    import bisect
    k = bisect.bisect_left(mzs, m * (1 - 20e-6))
    return k < len(mzs) and mzs[k] <= m * (1 + 20e-6)
for d in files:
    row = next(r for r in MANIFEST if r["id"] == d)
    import os
    if os.environ.get("SOURCE") == "prose":
        import gzip
        from compare import PROSE
        sa, peps = {}, {}
        votes = collections.Counter()
        for sd in (1, 42, 137):
            with open(PROSE / d / "base" / f"s{sd}" / "target.tsv") as f:
                for r in csv.DictReader(f, delimiter="\t"):
                    if float(r["q-value"]) <= 0.01:
                        scan = int(r["PSMId"].rsplit("_", 2)[-2]); votes[scan] += 1
                        peps[scan] = r["peptide"].split(".", 1)[1].rsplit(".", 1)[0]
        sa = {k: 1 for k, v in votes.items() if v >= 2}
    else:
        sa = accepted(SAGE / d)
        peps = {}
        with open(SAGE / d / "input.pin") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if r["rank"] == "1": peps[int(r["ScanNr"])] = r["Peptide"][2:-2] if r["Peptide"].startswith("-.") else r["Peptide"]
    scans = sorted(sa); random.Random(1).shuffle(scans); scans = set(scans[:N])
    exp = p.MSExperiment(); p.MzMLFile().load(str(B / f"prose-andes-reproduction/broad_benchmark/data/{d}.mzML"), exp)
    tot = collections.Counter()
    for spec in exp:
        if spec.getMSLevel() != 2: continue
        scan = int(spec.getNativeID().split("scan=")[-1])
        if scan not in scans: continue
        s = p.MSSpectrum(spec); s.sortByPosition()
        z = spec.getPrecursors()[0].getCharge()
        theo = ions(p.AASequence.fromString(peps[scan]), z)
        raw = [pk.getMZ() for pk in s]
        here = [m for m in theo if present(raw, m)]
        tot["matchable"] += len(here)
        for name, kw in VARIANTS.items():
            dd = p.MSSpectrum(s)
            p.Deisotoper.deisotopeAndSingleCharge(dd, 20.0, True, 1, 3, False, 3, 10, False, False, False, kw["use_decreasing_model"], kw["start_intensity_check"])
            kept = [pk.getMZ() for pk in dd]
            tot[name] += sum(1 for m in here if not present(kept, m))
    print(f"{d:22s} {row['label']:7s} spectra {len(scans):4d} matchable ions {tot['matchable']:6d}  removed by: " +
          "  ".join(f"{k} {100*tot[k]/tot['matchable']:5.1f}%" for k in VARIANTS), flush=True)
