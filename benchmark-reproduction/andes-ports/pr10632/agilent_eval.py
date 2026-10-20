"""Agilent 6550 set (PXD069532, no ground truth): per-precursor quantities from an MSstats CSV.
Reports precursors quantified in all 3 runs of both conditions, the median CV per condition (linear intensities,
as ProteoBench), and agreement of log2(NPS/OCBS) with FragPipe-IonQuant MBR (combined_ion.tsv) on unmodified
precursors present in both: Pearson r and median |difference|."""
import sys, numpy as np, pandas as pd
FP = ["NPS_4", "NPS_5", "NPS_6", "OCBS_10", "OCBS_11", "OCBS_12"]
fp = pd.read_csv(sys.argv[1], sep="\t", low_memory=False)
fp = fp[fp["Assigned Modifications"].isna()]
fp["key"] = fp["Peptide Sequence"] + "/" + fp["Charge"].astype(str)
fpi = fp.groupby("key")[[s + " Intensity" for s in FP]].sum().replace(0, np.nan)
fp_fc = np.log2(fpi.iloc[:, :3].mean(axis=1) / fpi.iloc[:, 3:].mean(axis=1)).dropna()
print("file\tprecursors_6of6\tCV_NPS\tCV_OCBS\tCV_median\tn_vs_fragpipe\tpearson_log2FC\tmedian_abs_dlog2FC")
for p in sys.argv[2:]:
    d = pd.read_csv(p)
    d["key"] = d["PeptideSequence"] + "/" + d["PrecursorCharge"].astype(str)
    m = d.groupby(["key", "Reference"])["Intensity"].sum().unstack()
    m = m[[c for c in sorted(m.columns)]].replace(0, np.nan)
    full = m.dropna()
    a, b = full.iloc[:, :3], full.iloc[:, 3:]
    cva = (a.std(axis=1, ddof=1) / a.mean(axis=1)); cvb = (b.std(axis=1, ddof=1) / b.mean(axis=1))
    fc = np.log2(a.mean(axis=1) / b.mean(axis=1))
    unmod = fc[~fc.index.str.contains(r"\(")]
    j = pd.concat([unmod.rename("o"), fp_fc.rename("f")], axis=1, join="inner")
    print(f"{p}\t{len(full)}\t{cva.median():.4f}\t{cvb.median():.4f}\t{pd.concat([cva, cvb]).median():.4f}\t"
          f"{len(j)}\t{j['o'].corr(j['f']):.4f}\t{(j['o'] - j['f']).abs().median():.4f}", flush=True)
