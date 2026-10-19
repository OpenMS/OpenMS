#!/usr/bin/env python3
"""Run one ProSE search (TOPP tool) for the PR #10635 benchmark and rescore its PIN with external Percolator.

Parameters are the frozen #10364 'native_dedup:<dataset>:default' parameters (jobs.json of the
prose-andes-reproduction bundle) plus the overrides of the arm. Parameters the build does not know are skipped and logged.

usage: run_search.py --bin <openms bin dir> --jobs <jobs.json> --dataset <id> --mzml <file> --db <fasta>
                     --out <dir> [--set key=value ...] [--percolator <exe>] [--seeds 1,42,137]
"""
import argparse, json, os, subprocess, sys, time
import xml.etree.ElementTree as ET
from pathlib import Path


def ini_params(prose):
    """Names of all Search:* parameters of this ProSE build."""
    out = subprocess.run([prose, "-write_ini", "-"], check=True, capture_output=True, text=True).stdout
    root = ET.fromstring(out)
    names = set()

    def walk(node, prefix):
        for c in node:
            if c.tag == "NODE":
                walk(c, prefix + [c.get("name")])
            elif c.tag in ("ITEM", "ITEMLIST"):
                names.add(":".join(prefix + [c.get("name")]))
    walk(root, [])
    # strip "ProSE:1:" prefix
    return {n.split(":", 2)[2] for n in names if n.count(":") >= 2}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--jobs", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--mzml", required=True)
    ap.add_argument("--db", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--set", action="append", default=[], help="key=value or key=v1|v2 for lists")
    ap.add_argument("--percolator")
    ap.add_argument("--seeds", default="1,42,137")
    ap.add_argument("--threads", default="4")
    a = ap.parse_args()

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    prose = str(Path(a.bin) / "ProSE")
    job = next(j for j in json.load(open(a.jobs)) if j["id"] == f"native_dedup:{a.dataset}:default")
    params = dict(job["parameters"])
    for s in a.set:
        k, v = s.split("=", 1)
        params[k] = v.split("|") if "|" in v else v
    known = ini_params(prose)
    cmd = [prose, "-in", a.mzml, "-database", a.db, "-out_idxml", str(out / "search.idXML"),
           "-out_pin", str(out / "search.pin"), "-threads", a.threads, "-no_summary"]
    skipped = []
    for k, v in params.items():
        key = "Search:" + k
        if key not in known:
            skipped.append(k)
            continue
        cmd.append("-" + key)
        cmd += [str(x) for x in v] if isinstance(v, list) else [str(v)]
    (out / "search.command.json").write_text(json.dumps({"command": cmd, "skipped_parameters": skipped}, indent=1))
    t0 = time.monotonic()
    with open(out / "search.log", "w") as log:
        subprocess.run(["/usr/bin/time", "-v"] + cmd, check=True, stdout=log, stderr=subprocess.STDOUT,
                       env=dict(os.environ, OMP_NUM_THREADS=a.threads))
    (out / "search.seconds").write_text(f"{time.monotonic() - t0:.2f}\n")

    if a.percolator:
        for seed in a.seeds.split(","):
            d = out / f"s{seed}"
            d.mkdir(exist_ok=True)
            with open(d / "percolator.log", "w") as log:
                subprocess.run([a.percolator, "-Y", "-U", "--seed", seed, "--results-psms", str(d / "target.tsv"),
                                "--decoy-results-psms", str(d / "decoy.tsv"), str(out / "search.pin")],
                               check=True, stdout=log, stderr=subprocess.STDOUT)


if __name__ == "__main__":
    main()
