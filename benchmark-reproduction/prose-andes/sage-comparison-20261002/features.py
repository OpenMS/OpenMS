"""Median ProSE rank-1 features: spectra both engines accept, Sage-only spectra where ProSE's rank 1 is Sage's peptide
(rank1_not_accepted), and rank-1 decoys. Also Sage's own features for the same groups."""
import collections, csv, gzip, json, statistics, sys
from compare import accepted, candidates, classify, SAGE, PROSE, MANIFEST, group_of, sequence
arm = sys.argv[1]; groups = sys.argv[2:]
PF = ["score", "hyperscore_zscore", "delta_score", "ln_num_candidates", "matched_prefix_ions", "matched_suffix_ions",
      "longest_peptide_ion_sequence", "matched_ion_current_fraction", "complementary_ions_fraction", "absdm", "isotope_error", "peplen"]
SF = ["ln(hyperscore)", "ln(-poisson)", "ln(delta_next)", "matched_peaks", "longest_b", "longest_y", "ln(matched_intensity_pct)", "ln(precursor_ppm)"]
for g in groups:
    vals = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in MANIFEST:
        d = r["id"]
        if group_of(d) != g: continue
        sa, pa = accepted(SAGE / d), accepted(PROSE / d / arm)
        pc = candidates(PROSE / d / arm / "native.pin", "score")
        cls = {}
        for s, q in sa.items():
            c = classify(s, q, pa, pc)
            cls[s] = "both" if c == "same" else c
        rank1 = {}
        with gzip.open(PROSE / d / arm / "native.pin.gz", "rt") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                s = int(row["ScanNr"])
                if s not in rank1 or float(row["score"]) > float(rank1[s]["score"]): rank1[s] = row
        for s, row in rank1.items():
            k = "decoy" if row["Label"] == "-1" else cls.get(s)
            if k in ("both", "rank1_not_accepted", "decoy"):
                for f_ in PF: vals[k][f_].append(float(row[f_]))
        with open(SAGE / d / "input.pin") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                if row["rank"] != "1": continue
                s = int(row["ScanNr"]); k = "decoy" if row["Label"] == "-1" else cls.get(s)
                if k in ("both", "rank1_not_accepted", "decoy"):
                    for f_ in SF: vals[k]["sage:" + f_].append(float(row[f_]))
    print(f"\n== {g} ({arm}); n = " + ", ".join(f"{k} {len(vals[k]['score'])}" for k in ("both", "rank1_not_accepted", "decoy")))
    print(f"{'feature':34s}{'both':>10s}{'rank1_not_acc':>15s}{'decoy rank1':>13s}")
    for f_ in PF + ["sage:" + x for x in SF]:
        print(f"{f_:34s}" + "".join(f"{statistics.median(vals[k][f_]) if vals[k][f_] else float('nan'):>{w}.3f}" for k, w in (("both", 10), ("rank1_not_accepted", 15), ("decoy", 13))))
