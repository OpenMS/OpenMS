"""Collect the per-arm evaluation files into one table per database (results/summary_<db>.tsv)."""
import csv, glob, os, sys

B = "/home/user/bench"
ARMS = ["base_off", "pr_off", "q0.01", "q0.005", "q0.001", "pr_off_fdr0.01", "q0.001_fdr0.01", "pr_off_noprotfdr", "q0.001_noprotfdr"]
for db in ("hye", "entrap"):
    rows = []
    for arm in ARMS:
        d = f"{B}/lfq/{db}/{arm}"
        if not os.path.exists(f"{d}/eval_transfers.tsv"):
            continue
        r = {"db": db, "arm": arm}
        t = list(csv.DictReader(open(f"{d}/eval_transfers.tsv"), delimiter="\t"))[0]
        r.update({k: v for k, v in t.items() if k != "label"})
        for kind in ("msstats", "rtlocal.msstats"):
            f = f"{d}/eval_pb_{kind}.tsv"
            if os.path.exists(f) and os.path.getsize(f):
                p = list(csv.DictReader(open(f), delimiter="\t"))[0]
                tag = "none" if kind == "msstats" else "rtlocal"
                for k in ("precursors", "median_abs_eps", "median_abs_eps_ECOLI", "CV_median"):
                    r[f"pb_{tag}_{k}"] = p[k]
        log = open(f"{d}/lfq.log").read()
        for line in log.splitlines():
            if "PIP-ECHO donor gate" in line:
                r["gate_line"] = line.strip()
            if line.startswith("PIP-ECHO: added"):
                r["pipecho_added"] = line.strip()
        rows.append(r)
    if not rows:
        continue
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(f"{B}/results/summary_{db}.tsv", "w") as f:
        w = csv.DictWriter(f, keys, delimiter="\t", restval="", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print(f"results/summary_{db}.tsv", len(rows), "arms")
