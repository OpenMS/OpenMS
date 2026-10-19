"""Compare EntrapmentEvaluator with the frozen harness logic (frozen_eval.py) on the same Percolator results.

  compare.py <bin_dir> <velos|shuffled> <te_ratio> <fasta|-> <run_dir>...

For every run_dir/s{1,42,137}: frozen_eval.py on target.tsv, EntrapmentEvaluator (-te_ratio) on target.idXML.
Prints per-file rows with both N, N_E, FDP and whether they agree, then the pooled FDP of all inputs."""
import os, subprocess, sys, tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
bin_dir, mode, te, fasta, runs = Path(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4], [Path(p) for p in sys.argv[5:]]
tsvs = [r / f"s{s}" / "target.tsv" for r in runs for s in (1, 42, 137)]
ids = [t.with_suffix(".idXML") for t in tsvs]
fr = subprocess.run([sys.executable, "-I", HERE / "frozen_eval.py", mode] + ([fasta] if mode == "velos" else []) + tsvs,
                    check=True, capture_output=True,
                   env=dict(os.environ, LD_LIBRARY_PATH=str(Path(bin_dir).parent / "lib"), OPENMS_DATA_PATH=str(Path(bin_dir).parent / "share/OpenMS")), text=True).stdout
frozen = {}
for line in fr.splitlines():
    if line.startswith("#"): print(line); continue
    f = line.split("\t")
    if f[0] == "file": continue
    frozen[f[0]] = (int(f[1]), int(f[2]), float(f[-1]))
with tempfile.NamedTemporaryFile(suffix=".tsv") as out:
    subprocess.run([bin_dir / "EntrapmentEvaluator", "-in"] + ids + ["-te_ratio", te, "-out", out.name, "-no_progress"],
                   check=True, capture_output=True,
                   env=dict(os.environ, LD_LIBRARY_PATH=str(Path(bin_dir).parent / "lib"), OPENMS_DATA_PATH=str(Path(bin_dir).parent / "share/OpenMS")))
    ee = {}
    for line in open(out.name).read().splitlines()[1:]:
        f = line.split("\t")
        if f[1] == "all": ee[f[0]] = (int(f[3]), int(f[4]), float(f[5]))
print("file\tfrozen_N\tfrozen_NE\tfrozen_FDP\tEE_N\tEE_NE\tEE_FDP\tsame")
ok = True
for t, i in zip(tsvs, ids):
    a, b = frozen[str(t)], ee[str(i) if len(ids) > 1 else i.name]
    same = a[0] == b[0] and a[1] == b[1] and abs(a[2] - b[2]) < 5e-5  # EntrapmentEvaluator prints 4 significant digits
    ok &= same
    print(f"{t.parent.parent.name}/{t.parent.name}\t{a[0]}\t{a[1]}\t{a[2]:.5f}\t{b[0]}\t{b[1]}\t{b[2]:.5f}\t{same}")
n, ne = sum(v[0] for v in frozen.values()), sum(v[1] for v in frozen.values())
tot = ee.get("all files")
print(f"pooled\tfrozen N {n} N_E {ne}\tEE {tot}\tall_rows_identical {ok}")
