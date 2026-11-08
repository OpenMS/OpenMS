"""Raw 01 design run: pick target and decoy CNN winners in matched z bands (target/decoy labels only, no entrapment labels)."""
import numpy as np, pandas as pd
c = pd.read_csv("../pipeline/crossfit_it0/crossfit.tsv", sep="\t")
w = c[c.precursor_best == 1].copy()
for f, g in w[w.is_decoy == 1].groupby("fold"):
    med = g.logit.median(); mad = 1.4826 * (g.logit - med).abs().median()
    m = w.fold == f; w.loc[m, "z"] = (w.loc[m, "logit"] - med) / mad
bands = [(10, 99), (8, 10), (6, 8), (4, 6), (0, 2)]
rows = []
for lo, hi in bands:
    for d in (0, 1):
        g = w[(w.is_decoy == d) & (w.z >= lo) & (w.z < hi)]
        s = g.sample(min(len(g), 60), random_state=11)
        for _, r in s.iterrows():
            rows.append((int(r.precursor_id), d, f"{lo}-{hi}", r.z, r.sample_id))
s = pd.DataFrame(rows, columns=["precursor_id", "is_decoy", "band", "z", "sample_id"])
print(s.groupby(["band", "is_decoy"]).size().unstack().to_string())
s.to_csv("selection.tsv", sep="\t", index=False)
open("precursors.txt", "w").write("# raw01 design run, matched z bands\n" + "\n".join(map(str, s.precursor_id)) + "\n")
