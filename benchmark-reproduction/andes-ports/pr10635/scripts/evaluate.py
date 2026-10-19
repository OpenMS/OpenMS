#!/usr/bin/env python3
"""Evaluate the PR #10635 (annotate:isotope_mimic) benchmark.

Layout: <runs>/<build>/<db>/<arm>/<dataset>/{search.pin, s<seed>/target.tsv, s<seed>/decoy.tsv}
  arms: off (Deamidated N/Q variable, isotope_mimic=false), on (isotope_mimic=true), nodeam (frozen #10364 mods),
        drop (off PIN without the rows the on-PIN flags, rescored again: the andes #102 behaviour)
Derived arm: post = off accepted PSMs minus flagged ones (filter after rescoring).

Per file and seed: target PSMs at q <= 0.01; entrapment PSMs (every protein ENTRAP_) with combined FDP
2 N_E / N (r = 1, shuffled mass-matched entrapment); class splits (flagged / deamidated unflagged / no deamidation)
with PSMs, entrapment FDP and decoys at the same q-value cut (D/T); precursor residual of accepted deamidated PSMs.
z = diff / sqrt(var_arm/3 + var_base/3) over seeds; up z >= 2, down z <= -2.

usage: evaluate.py <runs dir> <out dir>
"""
import csv, json, math, statistics, sys
from collections import defaultdict
from pathlib import Path

SEEDS = ["1", "42", "137"]
GROUPS = {"velos": "Velos CID", "hfx": "HF-X HCD", "astral": "Astral HCD", "lumos_lfq": "Lumos HCD LFQ"}
DATASETS = ["velos_125_R1", "velos_5000_R2", "velos_25000_R3", "hfx_A2", "hfx_B1", "hfx_B3",
            "astral_A2", "astral_B1", "astral_B3", "lumos_lfq_5192", "lumos_lfq_5194", "lumos_lfq_5199"]
DEAM = "[+0.984016]"


def group(ds):
    return next(v for k, v in GROUPS.items() if ds.startswith(k))


def read_pin(p):
    """{(SpecId, Peptide): row dict} with charge, residual (peptide - spectrum, Da), flag, isotope error."""
    rows = {}
    with open(p) as f:
        r = csv.reader(f, delimiter="\t")
        h = next(r)
        ix = {n: i for i, n in enumerate(h)}
        zc = [(int(n[6:]), ix[n]) for n in h if n.startswith("charge") and n[6:].isdigit()]
        for row in r:
            z = next(c for c, i in zc if row[i] == "1")
            rows[(row[ix["SpecId"]], row[ix["Peptide"]])] = dict(
                label=int(row[ix["Label"]]), z=z, err=-float(row[ix["dm"]]) * z,
                flag=int(row[ix["isotope_mimic"]]) if "isotope_mimic" in ix else None,
                iso=int(float(row[ix["isotope_error"]])), proteins=row[ix["Proteins"]], score=float(row[ix["score"]]),
                raw=[v for n, v in zip(h, row) if n != "isotope_mimic"])
    return rows


def read_perc(p, qmax=0.01):
    out = []
    with open(p) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if float(row["q-value"]) <= qmax:
                # Percolator writes further proteins as extra tab-separated fields
                out.append((row["PSMId"], row["peptide"], ";".join([row["proteinIds"]] + (row.get(None) or []))))
    return out


def is_entrap(proteins):
    acc = [a for a in proteins.replace("\t", ";").split(";") if a]
    return bool(acc) and all(a.startswith("ENTRAP_") for a in acc)


def cls(pep, flag):
    if flag:
        return "flagged"
    return "deam" if DEAM in pep else "plain"


def stats(accepted, decoys, flags):
    """accepted: list of (specid, pep, proteins) targets at q<=0.01; decoys: same for decoys."""
    d = defaultdict(lambda: dict(n=0, ne=0, nd=0))
    for s, pep, prot in accepted:
        c = cls(pep, flags.get((s, pep), 0))
        for k in (c, "all"):
            d[k]["n"] += 1
            d[k]["ne"] += is_entrap(prot)
    for s, pep, prot in decoys:
        c = cls(pep, flags.get((s, pep), 0))
        for k in (c, "all"):
            d[k]["nd"] += 1
    return {k: dict(v) for k, v in d.items()}


def zcmp(a, b):
    ma, mb = statistics.mean(a), statistics.mean(b)
    se = math.sqrt(statistics.variance(a) / 3 + statistics.variance(b) / 3) if len(a) > 1 else 0
    z = (ma - mb) / se if se else (0.0 if ma == mb else math.copysign(99, ma - mb))
    return ma - mb, z, "up" if z >= 2 else "down" if z <= -2 else "flat"


def main():
    runs, out = Path(sys.argv[1]) / "pr", Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    per_seed, checks, residual, competitor = [], [], [], []
    for db in ["plain", "entrap"]:
        for ds in DATASETS:
            base = runs / db
            if not (base / "on" / ds / "search.pin").exists():
                continue
            pin_on = read_pin(base / "on" / ds / "search.pin")
            pin_off = read_pin(base / "off" / ds / "search.pin")
            flags = {k: v["flag"] for k, v in pin_on.items()}
            same = pin_on.keys() == pin_off.keys() and all(pin_on[k]["raw"] == pin_off[k]["raw"] for k in pin_on)
            checks.append(dict(db=db, dataset=ds, check="off PIN == on PIN without isotope_mimic column", ok=same))
            # recompute the rule from the PIN residual: shift = 0.019339 Da per deamidation
            gap = 1.0033548378 - 0.984016
            mism = sum(1 for (s, pep), v in pin_on.items()
                       if v["flag"] != int(pep.count(DEAM) > 0 and abs(v["err"] + gap * pep.count(DEAM)) < abs(v["err"])))
            checks.append(dict(db=db, dataset=ds, check=f"flag == rule recomputed from PIN dm (mismatches {mism})",
                               ok=mism == 0))
            nrow = len(pin_on)
            nflag = sum(flags.values())
            ndeam = sum(DEAM in k[1] for k in pin_on)
            checks.append(dict(db=db, dataset=ds, check=f"PIN rows {nrow}, deamidated {ndeam}, flagged {nflag}", ok=True))
            for arm in ["nodeam", "off", "on", "drop"]:
                for seed in SEEDS:
                    d = base / arm / ds / f"s{seed}"
                    if not (d / "target.tsv").exists():
                        continue
                    acc = read_perc(d / "target.tsv")
                    dec = read_perc(d / "decoy.tsv")
                    st = stats(acc, dec, flags if arm != "nodeam" else {})
                    per_seed.append(dict(db=db, dataset=ds, arm=arm, seed=seed, stats=st))
                    if arm == "off":
                        accp = [a for a in acc if not flags.get((a[0], a[1]), 0)]
                        decp = [a for a in dec if not flags.get((a[0], a[1]), 0)]
                        per_seed.append(dict(db=db, dataset=ds, arm="post", seed=seed, stats=stats(accp, decp, flags)))
                        # precursor residual of accepted deamidated PSMs (andes #101 metric)
                        b = defaultdict(int)
                        for s, pep, _ in acc:
                            if DEAM not in pep:
                                continue
                            e = pin_off[(s, pep)]["err"] * 1000
                            b["deam"] += 1
                            b["near0"] += abs(e) <= 5
                            b["mispick"] += -26 <= e <= -12
                            b["flagged"] += flags[(s, pep)]
                            b[f"iso{pin_off[(s, pep)]['iso']}"] += 1
                        residual.append(dict(db=db, dataset=ds, seed=seed, **b))
                        if seed == "1":
                            # does the unmodified form compete for the flagged spectrum?
                            by_spec = defaultdict(list)
                            for (s, pep), v in pin_off.items():
                                by_spec[s].append((pep, v))
                            c = defaultdict(int)
                            for s, pep, _ in acc:
                                if not flags[(s, pep)]:
                                    continue
                                c["flagged_accepted"] += 1
                                plain = pep.replace(DEAM, "")
                                alt = [v for p2, v in by_spec[s] if p2 == plain]
                                if alt:
                                    c["unmodified_in_candidates"] += 1
                                    c["unmodified_scores_higher"] += alt[0]["score"] > pin_off[(s, pep)]["score"]
                                if is_entrap(pin_off[(s, pep)]["proteins"]):
                                    c["entrapment"] += 1
                            competitor.append(dict(db=db, dataset=ds, **c))
    json.dump(dict(per_seed=per_seed, checks=checks, residual=residual, competitor=competitor),
              open(out / "raw.json", "w"), indent=0)

    # per-file table
    def get(db, ds, arm, cl="all"):
        v = [r["stats"].get(cl, dict(n=0, ne=0, nd=0)) for r in per_seed
             if r["db"] == db and r["dataset"] == ds and r["arm"] == arm]
        return v

    lines = ["db\tdataset\tgroup\tarm\tpsms_s1/s42/s137\tmean\tdelta_vs_off\tz\tverdict\tfdp_combined_pct\t"
             "flagged_n\tflagged_fdp_pct\tflagged_decoy_pct\tdeam_unflagged_n\tdeam_unflagged_fdp_pct\t"
             "deam_unflagged_decoy_pct\tplain_n\tplain_fdp_pct\tplain_decoy_pct"]
    group_rows = defaultdict(lambda: defaultdict(list))
    for db in ["plain", "entrap"]:
        for ds in DATASETS:
            off = [x["n"] for x in get(db, ds, "off")]
            if len(off) < 3:
                continue
            for arm in ["nodeam", "off", "on", "post", "drop"]:
                v = get(db, ds, arm)
                if len(v) < 3:
                    continue
                n = [x["n"] for x in v]
                diff, z, verdict = zcmp(n, off)
                fdp = 100 * 2 * sum(x["ne"] for x in v) / sum(n) if db == "entrap" else float("nan")

                def cstat(cl):
                    c = get(db, ds, arm, cl)
                    cn = sum(x["n"] for x in c)
                    return (cn / 3, 100 * 2 * sum(x["ne"] for x in c) / cn if cn and db == "entrap" else float("nan"),
                            100 * sum(x["nd"] for x in c) / cn if cn else float("nan"))
                f, du, pl = cstat("flagged"), cstat("deam"), cstat("plain")
                lines.append("\t".join(map(str, [db, ds, group(ds), arm, "/".join(map(str, n)),
                                                 f"{statistics.mean(n):.1f}", f"{diff:+.1f}", f"{z:.1f}", verdict,
                                                 f"{fdp:.2f}", f"{f[0]:.1f}", f"{f[1]:.2f}", f"{f[2]:.2f}",
                                                 f"{du[0]:.1f}", f"{du[1]:.2f}", f"{du[2]:.2f}",
                                                 f"{pl[0]:.1f}", f"{pl[1]:.2f}", f"{pl[2]:.2f}"])))
                g = group_rows[(db, group(ds), arm)]
                g["mean"].append(statistics.mean(n))
                g["verdict"].append(verdict)
                g["n"].append(sum(n))
                g["ne"].append(sum(x["ne"] for x in v))
                for cl, t in (("flagged", f), ("deam", du), ("plain", pl)):
                    cc = get(db, ds, arm, cl)
                    g[cl + "_n"].append(sum(x["n"] for x in cc))
                    g[cl + "_ne"].append(sum(x["ne"] for x in cc))
                    g[cl + "_nd"].append(sum(x["nd"] for x in cc))
    (out / "per_file.tsv").write_text("\n".join(lines) + "\n")

    glines = ["db\tgroup\tarm\tn_files\tmean_psms\tdelta_vs_off_pct\tup/flat/down\tfdp_combined_pct\test_true_psms\t"
              "flagged_psms\tflagged_fdp_pct\tflagged_decoy_pct\tdeam_unflagged_psms\tdeam_unflagged_fdp_pct\t"
              "deam_unflagged_decoy_pct\tplain_psms\tplain_fdp_pct\tplain_decoy_pct"]
    for db in ["plain", "entrap"]:
        for gname in list(GROUPS.values()) + ["All"]:
            offm = None
            for arm in ["off", "nodeam", "on", "post", "drop"]:
                keys = [(db, gname, arm)] if gname != "All" else [(db, g2, arm) for g2 in GROUPS.values()]
                gs = [group_rows[k] for k in keys if k in group_rows]
                if not gs:
                    continue
                means = [m for g in gs for m in g["mean"]]
                mean = statistics.mean(means) if gname != "All" else sum(means)
                if arm == "off":
                    offm = mean
                verd = [v for g in gs for v in g["verdict"]]
                tot = lambda key: sum(sum(g[key]) for g in gs)
                fdp = 100 * 2 * tot("ne") / tot("n") if db == "entrap" else float("nan")
                row = [db, gname, arm, len(means), f"{mean:.1f}",
                       f"{100 * (mean / offm - 1):+.2f}" if offm else "",
                       f"{verd.count('up')}/{verd.count('flat')}/{verd.count('down')}", f"{fdp:.2f}",
                       f"{mean * (1 - fdp / 100):.1f}" if db == "entrap" else "nan"]
                for cl in ("flagged", "deam", "plain"):
                    n = tot(cl + "_n")
                    row += [f"{n / 3 / (len(means) if gname != 'All' else 1):.1f}",
                            f"{100 * 2 * tot(cl + '_ne') / n:.2f}" if n and db == "entrap" else "nan",
                            f"{100 * tot(cl + '_nd') / n:.2f}" if n else "nan"]
                glines.append("\t".join(map(str, row)))
    (out / "per_group.tsv").write_text("\n".join(glines) + "\n")

    rl = ["db\tdataset\tseed\taccepted_deamidated\tnear0_|err|<=5mDa\tmispick_-26..-12mDa\tflagged\tiso_errors"]
    for r in residual:
        isos = ",".join(f"{k}:{v}" for k, v in sorted(r.items()) if k.startswith("iso"))
        rl.append("\t".join(map(str, [r["db"], r["dataset"], r["seed"], r.get("deam", 0), r.get("near0", 0),
                                       r.get("mispick", 0), r.get("flagged", 0), isos])))
    (out / "residual_off.tsv").write_text("\n".join(rl) + "\n")
    cl_ = ["db\tdataset\tflagged_accepted_s1\tunmodified_form_in_candidates\tunmodified_scores_higher\tentrapment"]
    for r in competitor:
        cl_.append("\t".join(map(str, [r["db"], r["dataset"], r.get("flagged_accepted", 0),
                                        r.get("unmodified_in_candidates", 0), r.get("unmodified_scores_higher", 0),
                                        r.get("entrapment", 0)])))
    (out / "competitor_off.tsv").write_text("\n".join(cl_) + "\n")
    (out / "checks.tsv").write_text("db\tdataset\tcheck\tok\n" + "\n".join(
        f"{c['db']}\t{c['dataset']}\t{c['check']}\t{c['ok']}" for c in checks) + "\n")
    print((out / "per_group.tsv").read_text())


if __name__ == "__main__":
    main()
