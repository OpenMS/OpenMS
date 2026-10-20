"""Compare the top-level features (rt, mz, charge, intensity, quality) of two -feat_dir checkpoints of one run."""
import sys, pyarrow.parquet as pq
cols = ["rt", "mz", "charge", "intensity", "quality", "depth"]
def load(p):
    t = pq.read_table(p, columns=cols).to_pandas()
    return t[t.depth == 0].drop(columns="depth").sort_values(["mz", "rt", "charge", "intensity"]).reset_index(drop=True)
a, b = load(sys.argv[1]), load(sys.argv[2])
if len(a) == len(b) and (a.values == b.values).all():
    print(f"identical ({len(a)} features)")
else:
    m = a.merge(b, on=["rt", "mz", "charge", "intensity"], how="outer", indicator=True)
    print(f"DIFFER: {len(a)} vs {len(b)} features, {(m._merge == 'both').sum()} identical in rt/mz/charge/intensity")
