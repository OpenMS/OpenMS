"""Feature ablation for #9975 without new searches: drop PeptDeep features from a finished arm's PIN and rescore.

For each dataset and arm, writes results/<dataset>/<arm>-<ablation>/ with the reduced input.pin, Percolator
3.09.0 (-Y -U, seeds 1/42/137) results as in run.py, and a summary.json (seed_psms, mean_psms) that
evaluate_split.compare()/counts() can read. Ablations:
  nort  without rt_abs_error
  noms2 without ms2_cosine, ms2_spectral_angle, ms2_pearson, ms2_frac_pred_found
"""
import argparse
import concurrent.futures
import csv
import json
import statistics
from pathlib import Path

import run

ROOT = Path(__file__).resolve().parent
ABLATIONS = {"nort": ["rt_abs_error"], "noms2": ["ms2_cosine", "ms2_spectral_angle", "ms2_pearson", "ms2_frac_pred_found"]}


def ablate(dataset, arm, ablation):
    src = ROOT / "results" / dataset / arm
    out = ROOT / "results" / dataset / f"{arm}-{ablation}"
    if (out / "summary.json").exists():
        return json.loads((out / "summary.json").read_text())
    out.mkdir(parents=True, exist_ok=True)
    saved = json.loads((src / "summary.json").read_text())
    assert run.sha(src / "input.pin") == saved["input_checks"]["sha256"]
    drop = ABLATIONS[ablation]
    with (src / "input.pin").open() as f, (out / "input.pin").open("w", newline="") as g:
        rows = csv.reader(f, delimiter="\t")
        header = next(rows)
        assert all(d in header for d in drop), (dataset, arm, drop)
        keep = [i for i, h in enumerate(header) if h not in drop]
        w = csv.writer(g, delimiter="\t", lineterminator="\n")
        w.writerow([header[i] for i in keep])
        n = 0
        for r in rows:
            # Proteins is the open-ended last column; keep everything from it on
            w.writerow([r[i] for i in keep[:-1]] + r[keep[-1]:])
            n += 1
    check = dict(sha256=run.sha(out / "input.pin"), candidate_rows=n)
    seeds = [run.rescore(out, s, check) for s in (1, 42, 137)]
    result = dict(dataset=dataset, arm=f"{arm}-{ablation}", source_arm=arm, dropped=drop, input_checks=check,
                  seed_psms=[s["accepted"] for s in seeds], mean_psms=statistics.mean(s["accepted"] for s in seeds),
                  native=saved["native"], native_tsv_sha256=saved["native_tsv_sha256"], native_pin_sha256=saved["native_pin_sha256"])
    (out / "summary.json").write_text(json.dumps(result, indent=2))
    (out / "input.pin").unlink()
    print(json.dumps(dict(dataset=dataset, arm=result["arm"], seed_psms=result["seed_psms"])), flush=True)
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arms", nargs="+", required=True)
    ap.add_argument("--ablations", nargs="+", default=list(ABLATIONS))
    ap.add_argument("--workers", type=int, default=2)
    args = ap.parse_args()
    jobs = [(d.parent.name, d.name, a) for arm in args.arms for d in sorted((ROOT / "results").glob(f"*/{arm}"))
            if (d / "summary.json").exists() for a in args.ablations]
    print("Planned:", len(jobs), flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        for f in concurrent.futures.as_completed([pool.submit(ablate, *j) for j in jobs]):
            f.result()
    print("Done", flush=True)


if __name__ == "__main__":
    main()
