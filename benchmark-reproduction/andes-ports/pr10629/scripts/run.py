"""Search one dataset with one database and rescore it, for the PR #10629 entrapment benchmark.

  run.py <bin_dir> <dataset> <database.fasta> <out_dir>

ProSE replays the frozen `native_dedup:<dataset>:default` parameters of OpenMS#10364 (jobs.json); keys the build does not
define are skipped and logged. ProSE runs with its default in-process Percolator (-Search:rescore true) and writes
  inproc.idXML  in-process Percolator q-values (one seed, the default), best hit per spectrum
  native.pin    ProSE's own scores, all reported hits (pre-rescoring), the input of external Percolator
External Percolator (OpenMS/THIRDPARTY Linux/x86_64, -Y -U, seeds 1/42/137) on input.pin (native.pin with SpecId
<dataset>_<scan>_<rank>) writes s<seed>/target.tsv, which
IDFileConverter (-score_type qvalue) turns into s<seed>/target.idXML for EntrapmentEvaluator."""
import json, os, subprocess, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORK = HERE.parent / "work"
PERCOLATOR = HERE.parent / "thirdparty/Linux/x86_64/Percolator/percolator"
JOBS = {j["id"]: j for j in json.load(open(HERE.parent / "frozen/z/prose-andes-reproduction/jobs.json"))}
MANIFEST = {r["id"]: r for r in json.load(open(WORK / "broad_benchmark/manifest.json"))}

bin_dir, dataset, fasta, out = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3]), Path(sys.argv[4])
out.mkdir(parents=True, exist_ok=True)
env = dict(os.environ, OMP_NUM_THREADS="4", OPENMS_DATA_PATH=str(bin_dir.parent / "share/OpenMS"),
           LD_LIBRARY_PATH=str(bin_dir.parent / "lib"))
log = open(out / "run.log", "a")


def run(cmd, **kw):
    print("$", " ".join(map(str, cmd)), file=log, flush=True)
    t = time.monotonic()
    subprocess.run(["/usr/bin/time", "-v"] + [str(c) for c in cmd], check=True, stdout=log, stderr=log, env=env, **kw)
    return time.monotonic() - t


# ProSE parameter names of this build
ini = out / "prose_defaults.ini"
if not ini.exists():
    subprocess.run([bin_dir / "ProSE", "-write_ini", ini], check=True, stdout=subprocess.DEVNULL, env=env)
known = ini.read_text()
params = JOBS[f"native_dedup:{dataset}:default"]["parameters"]
args, skipped = [], []
for k, v in params.items():
    if f'name="{k.split(":")[-1]}"' not in known or (":" in k and f'NODE name="{k.split(":")[0]}"' not in known):
        skipped.append(k); continue
    args.append(f"-Search:{k}")
    args += [str(x) for x in v] if isinstance(v, list) else [str(v)]
print("skipped (not defined by this build):", skipped, file=log, flush=True)
mzml = WORK / "broad_benchmark/data" / f"{dataset}.mzML"
if not (out / "native.pin").exists():
    secs = run([bin_dir / "ProSE", "-in", mzml, "-database", fasta, "-out_idxml", out / "inproc.idXML",
                "-out_pin", out / "native.pin", "-threads", "4", "-no_summary"] + args)
    print(f"prose_seconds {secs:.1f}", file=log, flush=True)
# SpecId as the frozen native_dedup/export.py writes it (<dataset>_<scan>_<rank>), which the frozen evaluation scripts
# parse; every other column is ProSE's own (no further PIN normalisation, proteins stay tab-separated)
if not (out / "input.pin").exists():
    rank = {}
    with open(out / "native.pin") as f, open(out / "input.pin.part", "w") as g:
        g.write(f.readline())
        for line in f:
            spec, rest = line.split("\t", 1)
            scan = rest.split("\t", 2)[1]
            rank[scan] = rank.get(scan, -1) + 1
            g.write(f"{dataset}_{scan}_{rank[scan]}\t{rest}")
    (out / "input.pin.part").rename(out / "input.pin")
for seed in (1, 42, 137):
    s = out / f"s{seed}"
    s.mkdir(exist_ok=True)
    if not (s / "target.tsv").exists():
        run([PERCOLATOR, "-Y", "-U", "--seed", seed, "--results-psms", s / "target.tsv",
             "--decoy-results-psms", s / "decoy.tsv", out / "input.pin"])
    if not (s / "target.idXML").exists():
        (s / "target.psms").write_bytes((s / "target.tsv").read_bytes())  # IDFileConverter detects the format by extension
        run([bin_dir / "IDFileConverter", "-in", s / "target.psms", "-score_type", "qvalue",
             "-out", s / "target.idXML"])
print("done", file=log, flush=True)
