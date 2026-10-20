"""Per-run counts of features and elution model fits (model_status) in ProteomicsLFQ -feat_dir checkpoints."""
import sys, glob, collections
import pyarrow.parquet as pq
print("run\tfeatures\tvalid_fit\tinvalid_fit\tno_fit_meta")
tot = collections.Counter()
for p in sorted(glob.glob(sys.argv[1] + "/*.featureParquet")):
    t = pq.read_table(p, columns=["metavalues", "depth"])
    # top-level features only (depth 0); other rows of the checkpoint are not features
    mv = [m for m, d in zip(t.column("metavalues").to_pylist(), t.column("depth").to_pylist()) if d == 0]
    c = collections.Counter()
    for m in mv:
        st = next((x["value"] for x in (m or []) if x["name"] == "model_status"), None)
        c["features"] += 1
        c["no_fit_meta" if st is None else ("valid_fit" if st.startswith("0") else "invalid_fit")] += 1
    tot.update(c)
    print(p.split("/")[-1].split("_fg")[0] + "\t" + "\t".join(str(c[k]) for k in ("features", "valid_fit", "invalid_fit", "no_fit_meta")))
print("total\t" + "\t".join(str(tot[k]) for k in ("features", "valid_fit", "invalid_fit", "no_fit_meta")))
