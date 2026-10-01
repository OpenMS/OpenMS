"""Per-file robustness of each candidate change against develop.

For every file: mean Percolator PSMs over seeds 1/42/137 for the arm and for develop, the difference, and a
z-score against the seed-to-seed spread (SE = sqrt(var_arm / 3 + var_dev / 3), sample variances over the
three seeds). A file counts as up (z >= 2), down (z <= -2) or flat. Also the native TDC difference, which
does not depend on Percolator.
"""
import csv
import json
import math
import statistics
import sys
from pathlib import Path

from evaluate import GROUPS, MANIFEST, group_of

ROOT = Path(__file__).resolve().parent
# candidate: list of (arm, files) — the first arm that has a summary for a file is used
CANDIDATES = {
    "#10377 (cutoff)": ["pr68"],
    "#10378 head 241f9a7 (priors on)": ["pr66_priors"],
    "#10378 b02b538 (priors on)": ["p_zauto", "f78_priors", "pr66_priors"],  # Da files: p_zauto (identical PIN); ppm files: unchanged code path
    "#10378 523214b (priors on, cross-fitted)": ["f78x_priors"],
    "#10379 head 16ff034 (default)": ["pr65_default"],
    "#10379 eb74991 (default)": ["f79_default", "h_frag", "base"],
    "#10379 eb74991 + mass_accuracy": ["f79_mass", "m_shift_frag"],
    "#10335 merged": ["pr10335"],
}


def load(d, a):
    p = ROOT / "results" / d / a / "summary.json"
    return json.loads(p.read_text()) if p.exists() else None


def main():
    rows = []
    for cand, arms in CANDIDATES.items():
        for r in MANIFEST:
            d = r["id"]
            s = next((load(d, a) for a in arms if load(d, a)), None)
            # #10379's default does not calibrate low-resolution files: those equal develop
            if s is None:
                continue
            b = load(d, "base")
            va, vb = statistics.variance(s["seed_psms"]), statistics.variance(b["seed_psms"])
            se = math.sqrt(va / 3 + vb / 3) or 1.0
            diff = s["mean_psms"] - b["mean_psms"]
            z = diff / se
            rows.append(dict(candidate=cand, dataset=d, group=group_of(d), arm=next(a for a in arms if load(d, a)),
                             develop=round(b["mean_psms"], 1), mean=round(s["mean_psms"], 1), delta_pct=round(100 * diff / b["mean_psms"], 2),
                             z=round(z, 1), verdict="up" if z >= 2 else "down" if z <= -2 else "flat",
                             native_delta=s["native"]["psms"] - b["native"]["psms"]))
    with (ROOT / "robustness.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
    for cand in CANDIDATES:
        sub = [x for x in rows if x["candidate"] == cand]
        print(f"\n{cand}: files up {sum(x['verdict']=='up' for x in sub)}, flat {sum(x['verdict']=='flat' for x in sub)}, "
              f"down {sum(x['verdict']=='down' for x in sub)} of {len(sub)}")
        for _, g in GROUPS:
            gs = [x for x in sub if x["group"] == g]
            if gs:
                print(f"  {g:22s} " + "  ".join(f"{x['delta_pct']:+6.2f}% (z {x['z']:+5.1f}, native {x['native_delta']:+d})" for x in gs))


if __name__ == "__main__":
    main()
