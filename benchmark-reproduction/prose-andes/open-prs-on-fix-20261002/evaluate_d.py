"""Do the open PRs' gains hold on develop with the #10391 deisotoping fix?

Each candidate is compared twice against its own baseline on the same files:
  old: the PR arm against develop f5ea2d04 (arm "base"), as benchmarked on 2026-10-01;
  new: the PR merged with develop 3a47278 against develop 3a47278 (arm "d_base").
Per file: mean Percolator PSMs over seeds 1/42/137, the difference and z = diff / SE with
SE = sqrt(var_arm / 3 + var_base / 3). A file is up (z >= 2), down (z <= -2) or flat. Group means are
means over the group's files of the seed means; the group change is the ratio of group means.
Files an arm does not search fall back to the listed arms (#10379's default does not calibrate
low-resolution files, so there it equals its baseline; checked on two files).
"""
import csv
import json
import math
import statistics
from pathlib import Path

from evaluate import GROUPS, MANIFEST, group_of

ROOT = Path(__file__).resolve().parent
# candidate: (old arms, old base, new arms, new base)
CANDIDATES = {
    "#10378 priors on (opt-in)": (["f78x_priors"], "base", ["d_priors"], "d_base"),
    "#10379 default (calibration auto)": (["f79_default", "base"], "base", ["d_79", "d_base"], "d_base"),
    "#10379 + mass_accuracy (opt-in)": (["f79_mass"], "base", ["d_79_mass"], "d_base"),
    "#10335 merged with develop": (["pr10335"], "base", ["d_35"], "d_base"),
}


def load(d, a):
    p = ROOT / "results" / d / a / "summary.json"
    return json.loads(p.read_text()) if p.exists() else None


def compare(d, arms, base):
    b = load(d, base)
    arm = next((a for a in arms if load(d, a)), None)
    if b is None or arm is None:
        return None
    s = load(d, arm)
    va, vb = statistics.variance(s["seed_psms"]), statistics.variance(b["seed_psms"])
    se = math.sqrt(va / 3 + vb / 3) or 1.0
    diff = s["mean_psms"] - b["mean_psms"]
    z = diff / se
    return dict(arm=arm, base=b["mean_psms"], mean=s["mean_psms"], delta_pct=100 * diff / b["mean_psms"], z=z,
                verdict="up" if z >= 2 else "down" if z <= -2 else "flat",
                native_delta=s["native"]["psms"] - b["native"]["psms"])


def main():
    rows = []
    for cand, (old_arms, old_base, new_arms, new_base) in CANDIDATES.items():
        for r in MANIFEST:
            d = r["id"]
            for when, arms, base in (("old", old_arms, old_base), ("new", new_arms, new_base)):
                c = compare(d, arms, base)
                if c:
                    rows.append(dict(candidate=cand, when=when, dataset=d, group=group_of(d), **c))
    with (ROOT / "holds_per_file.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        for x in rows:
            w.writerow({k: (round(v, 2) if isinstance(v, float) else v) for k, v in x.items()})
    out = []
    for cand in CANDIDATES:
        out.append(f"\n## {cand}")
        for when in ("old", "new"):
            sub = [x for x in rows if x["candidate"] == cand and x["when"] == when]
            if not sub:
                out.append(f"  {when}: no results")
                continue
            n = {v: sum(x["verdict"] == v for x in sub) for v in ("up", "flat", "down")}
            out.append(f"  {when}: files up {n['up']} / flat {n['flat']} / down {n['down']} of {len(sub)}")
        out.append(f"  {'group':22s} {'n':>2s} {'old base':>9s} {'old arm':>9s} {'old Δ':>8s}   {'new base':>9s} {'new arm':>9s} {'new Δ':>8s}  new per file (z)")
        for _, g in GROUPS:
            o = [x for x in rows if x["candidate"] == cand and x["when"] == "old" and x["group"] == g]
            nw = [x for x in rows if x["candidate"] == cand and x["when"] == "new" and x["group"] == g]
            if not o and not nw:
                continue
            def agg(sub):
                if not sub:
                    return float("nan"), float("nan"), float("nan")
                b = sum(x["base"] for x in sub) / len(sub)
                a = sum(x["mean"] for x in sub) / len(sub)
                return b, a, 100 * (a / b - 1)
            ob, oa, od = agg(o)
            nb, na, nd = agg(nw)
            per = " ".join(f"{x['delta_pct']:+.1f}({x['z']:+.1f})" for x in nw)
            out.append(f"  {g:22s} {len(nw or o):2d} {ob:9.1f} {oa:9.1f} {od:+7.2f}%   {nb:9.1f} {na:9.1f} {nd:+7.2f}%  {per}")
    text = "\n".join(out)
    (ROOT / "holds_summary.txt").write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
