#!/usr/bin/env python3
"""Second PR head (precursor:isotope_mimic=remove) against the arms of evaluate.py.

Arms compared, per db (plain = suite databases, entrap = entrapment databases):
  off (first head; identical to the merge-base), on (first head: flag + Percolator feature), drop (simulated andes:
  flagged PIN rows removed before Percolator), remove (second head, precursor:isotope_mimic=remove), annotate (second
  head, flag without Percolator feature; plain counts must equal off).
Also: second head with the option off against the merge-base (byte identity), and the precursor residual of accepted
deamidated PSMs per arm.

usage: evaluate_remove.py <runs dir> <out dir>"""
import hashlib, math, re, statistics, sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from evaluate import DATASETS, DEAM, SEEDS, group, is_entrap, read_perc, read_pin, zcmp  # noqa: E402

runs, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
ARMS = [("off", "pr"), ("on", "pr"), ("drop", "pr"), ("annotate", "pr2"), ("remove", "pr2")]
GAP = 1.0033548378 - 0.984016


def arm_dir(db, arm, build, ds):
    return runs / build / db / arm / ds


per_file = ["db\tdataset\tgroup\tarm\tpsms_s1/s42/s137\tmean\tdelta_vs_off\tz\tverdict\tentrapment_psms_3seeds"]
group_rows = defaultdict(lambda: defaultdict(list))
resid = defaultdict(lambda: defaultdict(float))
for db in ["plain", "entrap"]:
    for ds in DATASETS:
        res = {}
        for arm, build in ARMS:
            d = arm_dir(db, arm, build, ds)
            if not (d / "s137" / "target.tsv").exists():
                continue
            pin = read_pin(d / "search.pin")
            n, ne = [], 0
            for s in SEEDS:
                acc = read_perc(d / f"s{s}" / "target.tsv")
                n.append(len(acc))
                ne += sum(is_entrap(p) for _, _, p in acc)
                for sid, pep, _ in acc:
                    if DEAM not in pep:
                        continue
                    e = pin[(sid, pep)]["err"] * 1000
                    k = pep.count(DEAM)
                    r = resid[(db, group(ds), arm)]
                    r["deam"] += 1 / 3
                    r["near0"] += (abs(e) <= 5) / 3
                    r["mispick"] += (-26 <= e <= -12) / 3
                    r["mimic"] += (abs(e + GAP * 1000 * k) < abs(e)) / 3
            res[arm] = (n, ne)
        if "off" not in res:
            continue
        for arm, (n, ne) in res.items():
            diff, z, verdict = zcmp(n, res["off"][0])
            per_file.append("\t".join(map(str, [db, ds, group(ds), arm, "/".join(map(str, n)),
                                                 f"{statistics.mean(n):.1f}", f"{diff:+.1f}", f"{z:.1f}", verdict, ne])))
            g = group_rows[(db, group(ds), arm)]
            g["mean"].append(statistics.mean(n))
            g["verdict"].append(verdict)
            g["n"].append(sum(n))
            g["ne"].append(ne)
(out / "remove_per_file.tsv").write_text("\n".join(per_file) + "\n")

groups = ["Velos CID", "HF-X HCD", "Astral HCD", "Lumos HCD LFQ"]
lines = ["db\tgroup\tarm\tmean_psms\tdelta_vs_off_pct\tup/flat/down\tfdp_combined_pct\test_true_psms\t"
         "accepted_deamidated\tdeam_mispick_pct\tdeam_mimic_rule_pct\tdeam_near0_pct"]
for db in ["plain", "entrap"]:
    for gname in groups + ["All"]:
        offm = None
        for arm, _ in ARMS:
            keys = [(db, gname, arm)] if gname != "All" else [(db, g2, arm) for g2 in groups]
            gs = [group_rows[k] for k in keys if k in group_rows]
            if not gs or sum(len(g["mean"]) for g in gs) != (3 if gname != "All" else 12):
                continue
            means = [m for g in gs for m in g["mean"]]
            mean = statistics.mean(means) if gname != "All" else sum(means)
            if arm == "off":
                offm = mean
            verd = [v for g in gs for v in g["verdict"]]
            n = sum(sum(g["n"]) for g in gs)
            ne = sum(sum(g["ne"]) for g in gs)
            fdp = 100 * 2 * ne / n if db == "entrap" else float("nan")
            r = defaultdict(float)
            for k in keys:
                for kk, v in resid[k].items():
                    r[kk] += v
            dn = r["deam"] or 1
            lines.append("\t".join(map(str, [
                db, gname, arm, f"{mean:.1f}", f"{100 * (mean / offm - 1):+.2f}" if offm else "",
                f"{verd.count('up')}/{verd.count('flat')}/{verd.count('down')}", f"{fdp:.2f}",
                f"{mean * (1 - fdp / 100):.1f}" if db == "entrap" else "nan",
                f"{r['deam']:.1f}", f"{100 * r['mispick'] / dn:.1f}", f"{100 * r['mimic'] / dn:.1f}",
                f"{100 * r['near0'] / dn:.1f}"])))
(out / "remove_per_group.tsv").write_text("\n".join(lines) + "\n")

# second head with the option off vs merge-base
NORM = re.compile(rb'search_engine_version="[^"]*"|date="[^"]*"|creation_date="[^"]*"')


def h(p, norm=False):
    b = Path(p).read_bytes()
    return hashlib.sha256(NORM.sub(b"", b) if norm else b).hexdigest()


iso = ["dataset\tpin_identical\tidxml_identical_normalized\tpercolator_s1_s42_s137_identical"]
for ds in DATASETS:
    a, b = runs / "base" / "plain" / "off" / ds, runs / "pr2" / "plain" / "off" / ds
    if not (b / "s137" / "target.tsv").exists():
        continue
    iso.append(f"{ds}\t{h(a / 'search.pin') == h(b / 'search.pin')}\t"
               f"{h(a / 'search.idXML', True) == h(b / 'search.idXML', True)}\t"
               f"{all(h(a / f's{s}' / 'target.tsv') == h(b / f's{s}' / 'target.tsv') for s in SEEDS)}")
(out / "isolation_base_vs_pr2_off.tsv").write_text("\n".join(iso) + "\n")
print((out / "remove_per_group.tsv").read_text())
print((out / "isolation_base_vs_pr2_off.tsv").read_text())
