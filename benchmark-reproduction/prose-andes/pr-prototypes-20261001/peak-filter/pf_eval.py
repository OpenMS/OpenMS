"""Dense-spectrum peak quota (peakfilter_bench) against develop (base): per file and per group, plus identity checks."""
import csv, json, statistics, sys
from pathlib import Path
from evaluate import MANIFEST, group_of
ROOT = Path(__file__).resolve().parent
ARMS = sys.argv[1:] or ["pf", "pf_l10", "pf_l30", "wt100"]


def load(d, a):
    p = ROOT / "results" / d / a / "summary.json"
    return json.loads(p.read_text()) if p.exists() else None


def dense_share(d, a):
    log = ROOT / "results" / d / a / "search.log"
    for line in log.read_text().splitlines() if log.exists() else []:
        if "MS2 spectra are dense" in line:
            n, _, total = line.split("] ", 1)[1].split()[:3]
            return round(100 * int(n) / int(total), 1)
    return 0.0


rows, checks = [], []
for r in MANIFEST:
    d = r["id"]; b = load(d, "base")
    for a in ["pf_off", "pf"]:
        s = load(d, a)
        if s and (a == "pf_off" or r["fragment_unit"] == "Da"):
            same = s["native_pin_sha256"] == b["native_pin_sha256"] and s["native_tsv_sha256"] == b["native_tsv_sha256"]
            checks.append(dict(dataset=d, arm=a, identical_to_develop=same))
    for a in ARMS:
        s = load(d, a)
        if not s: continue
        se = ((statistics.variance(s["seed_psms"]) + statistics.variance(b["seed_psms"])) / 3) ** 0.5
        z = (s["mean_psms"] - b["mean_psms"]) / se if se else 0.0
        rows.append(dict(arm=a, dataset=d, group=group_of(d), dense_pct=dense_share(d, a) if a.startswith("pf") else "",
                         develop=round(b["mean_psms"], 1), mean=round(s["mean_psms"], 1),
                         delta_pct=round(100 * (s["mean_psms"] / b["mean_psms"] - 1), 2), z=round(z, 1),
                         verdict="up" if z >= 2 else "down" if z <= -2 else "flat",
                         native_delta_pct=round(100 * (s["native"]["psms"] / b["native"]["psms"] - 1), 1)))
with open(ROOT / "pf_per_file.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
for c in checks: print("identity", c)
for a in ARMS:
    sub = [x for x in rows if x["arm"] == a]
    if not sub: continue
    print(f"\n== {a}: up {sum(x['verdict']=='up' for x in sub)} / flat {sum(x['verdict']=='flat' for x in sub)} / down {sum(x['verdict']=='down' for x in sub)} of {len(sub)}")
    for g in dict.fromkeys(x["group"] for x in sub):
        gs = [x for x in sub if x["group"] == g]
        print(f"  {g:22s} dense {statistics.mean(x['dense_pct'] or 0 for x in gs):5.1f}%  "
              f"{statistics.mean(x['delta_pct'] for x in gs):+6.2f}%  native {statistics.mean(x['native_delta_pct'] for x in gs):+5.1f}%  "
              + " ".join(f"{x['delta_pct']:+.1f}({x['z']:+.1f})" for x in gs))
