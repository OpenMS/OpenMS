"""Post-mortem on spent raw 02: score distributions of human targets, entrapment targets and decoys (not a gate)."""
import sys
sys.path.insert(0, "../astral-entrapment-raw02-20261108")
import numpy as np, pandas as pd, evaluate as ev
lab = ev.library_labels("../astral-entrapment-20261107/lib/library.pqp")
c = pd.read_csv("../astral-entrapment-raw02-20261108/crossfit/crossfit.tsv", sep="\t")
w = c[c.precursor_best == 1].copy()
# pooled standardisation, same as -pooled_tdc
for f, g in w[w.is_decoy == 1].groupby("fold"):
    med = g.logit.median(); mad = 1.4826 * (g.logit - med).abs().median()
    m = w.fold == f; w.loc[m, "z"] = (w.loc[m, "logit"] - med) / mad
w["cls"] = np.where(w.is_decoy == 1, "decoy", np.where(w.precursor_id.map(lab) == True, "entrap", "human"))
print(w.cls.value_counts().to_string())
print("\nz quantiles per class")
print(w.groupby("cls").z.quantile([.5, .9, .99, .999]).unstack().round(2).to_string())
print("\nsurvival counts N(z>t) per class (all precursor winners, before pair competition)")
rows = []
for t in [2, 4, 6, 8, 10, 12, 15, 20]:
    s = w[w.z > t].cls.value_counts()
    rows.append((t, s.get("human", 0), s.get("entrap", 0), s.get("decoy", 0)))
print(pd.DataFrame(rows, columns=["t", "human", "entrap", "decoy"]).to_string(index=False))
w.to_parquet("winners.parquet")
