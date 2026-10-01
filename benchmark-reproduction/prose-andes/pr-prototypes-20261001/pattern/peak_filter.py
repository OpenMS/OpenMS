"""Peak-filter arms (peaks:window_top 40 / 100 instead of 20) against develop on the 14 high-resolution files."""
import csv, json, statistics
from pathlib import Path
from evaluate import MANIFEST, group_of
ROOT = Path(__file__).resolve().parent.parent
L = lambda d, a: json.loads((ROOT / "results" / d / a / "summary.json").read_text())


def z(a, b):
    s = ((statistics.variance(a["seed_psms"]) + statistics.variance(b["seed_psms"])) / 3) ** 0.5
    return (a["mean_psms"] - b["mean_psms"]) / s if s else 0.0


loss = {r["dataset"]: r for r in csv.DictReader(open(ROOT / "pattern" / "filter_loss.tsv"), delimiter="\t")}
rows = []
for r in MANIFEST:
    d = r["id"]
    if r["fragment_unit"] != "ppm": continue
    b = L(d, "base")
    row = dict(dataset=d, group=group_of(d), frag_tic_lost_top20_pct=loss[d]["tic_lost_top20_pct"],
               develop=round(b["mean_psms"], 1), develop_native=b["native"]["psms"])
    for a in ["wt40", "wt100"]:
        s = L(d, a)
        row.update({f"{a}": round(s["mean_psms"], 1), f"{a}_seeds": " ".join(map(str, s["seed_psms"])),
                    f"{a}_pct": round(100 * (s["mean_psms"] / b["mean_psms"] - 1), 2), f"{a}_z": round(z(s, b), 1),
                    f"{a}_native_pct": round(100 * (s["native"]["psms"] / b["native"]["psms"] - 1), 1)})
    rows.append(row)
with open(ROOT / "pattern" / "peak_filter.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
groups = {}
for x in rows: groups.setdefault(x["group"], []).append(x)
for g, xs in groups.items():
    print(g, *(f"{a} {statistics.mean(x[a + '_pct'] for x in xs):+.2f}" for a in ["wt40", "wt100"]),
          f"native100 {statistics.mean(x['wt100_native_pct'] for x in xs):+.1f}")
