import collections, csv, json, sys, statistics
from compare import accepted, candidates, SAGE, PROSE, B
d = sys.argv[1]; arm = sys.argv[2] if len(sys.argv) > 2 else "base"
sa, pa = accepted(SAGE / d), accepted(PROSE / d / arm)
pc = candidates(PROSE / d / arm / "native.pin", "score")
prov = {s["scan"]: s for s in json.loads((B / f"prose-andes-reproduction/broad_benchmark/provenance/{d}.json").read_text())["scans"]}
sage_rows = {}
with open(SAGE / d / "input.pin") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["rank"] == "1": sage_rows[int(r["ScanNr"])] = r
miss = [s for s in sa if s not in pa and s not in pc]
print(d, "Sage-accepted spectra without any ProSE candidate:", len(miss))
mp = collections.Counter(int(float(sage_rows[s]["matched_peaks"])) for s in miss)
print("Sage matched_peaks:", sorted(mp.items()))
both = [s for s in sa if pa.get(s) == sa[s]]
print("matched_peaks among shared:", statistics.median(int(float(sage_rows[s]["matched_peaks"])) for s in both))
for s in miss[:12]:
    r = sage_rows[s]
    print(s, "z", prov[s]["charge"], "mz", round(prov[s]["precursor_mz"], 3), "peptide", r["Peptide"][:60], "matched", r["matched_peaks"], "longest_b/y", r["longest_b"], r["longest_y"])
