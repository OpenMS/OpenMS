"""#9975 (PeptDeep rescoring features) rebased onto develop f6c680f, 2026-10-02.

Comparisons (per file mean Percolator PSMs over seeds 1/42/137, z = diff / SE with
SE = sqrt(var_arm / 3 + var_base / 3); up z >= 2, down z <= -2, else flat):
  pd_on   (#9975, peptdeep:enable=true, defaults: instrument QE, NCE auto) against re_auto (develop f6c680f).
  pd_inst (the same with peptdeep's own instrument group per instrument) against pd_on and re_auto.
Identity: pd_off (#9975 with its defaults, PeptDeep off) must equal re_auto byte for byte.
Entrapment: natural UPS1/yeast on Velos (r = 3.38), shuffled paired-mass (r = 1) on HF-X, Astral, TMTpro.
Cost: wall time, CPU time and peak RSS of the search (Percolator excluded).
"""
import csv
import json
import re
import statistics
from pathlib import Path

import evaluate_split as es
from evaluate import GROUPS, MANIFEST, group_of

ROOT = Path(__file__).resolve().parent
COMPARISONS = [
    ("#9975 PeptDeep features (defaults: instrument QE, NCE auto)", "pd_on", "re_auto"),
    ("#9975 + calibration on each spectrum's best target hit", "pd_cal", "re_auto"),
    ("calibration on each spectrum's best target hit, against #9975 as rebased", "pd_cal", "pd_on"),
    ("#9975 without rt_abs_error (Percolator rerun on the same PIN)", "pd_on-nort", "re_auto"),
    ("#9975 without the four MS2 features (Percolator rerun on the same PIN)", "pd_on-noms2", "re_auto"),
    ("calibration fix, without rt_abs_error", "pd_cal-nort", "re_auto"),
    ("calibration fix, without the four MS2 features", "pd_cal-noms2", "re_auto"),
    ("#9975 with peptdeep's instrument group per instrument, against QE", "pd_inst", "pd_on"),
    ("#9975 with peptdeep's instrument group per instrument, against develop", "pd_inst", "re_auto"),
]
IDENTITY = [("pd_off", "re_auto", "#9975 with PeptDeep off (its default) equals develop f6c680f")]
VELOS_ARMS = ["re_auto", "pd_on", "pd_on-nort", "pd_on-noms2", "pd_cal", "pd_cal-nort", "pd_cal-noms2", "pd_inst"]
SHUFFLED = {"HF-X HCD": ["hfx_A2", "hfx_B1", "hfx_B3"], "Astral HCD": ["astral_A2", "astral_B1", "astral_B3"],
            "Exploris 480 TMTpro": ["eclipse_tmtpro_10855", "eclipse_tmtpro_10858", "eclipse_tmtpro_10863"]}


def yields():
    rows, out = [], []
    for name, arm, base in COMPARISONS:
        sub = []
        for r in MANIFEST:
            c = es.compare(r["id"], arm, base)
            if c:
                sub.append(dict(comparison=name, arm=arm, baseline=base, dataset=r["id"], group=group_of(r["id"]), **c))
        rows += sub
        if not sub:
            out.append(f"\n## {name} ({arm} vs {base}): no results")
            continue
        n = {v: sum(x["verdict"] == v for x in sub) for v in ("up", "flat", "down")}
        total_b, total_a = sum(x["base"] for x in sub), sum(x["mean"] for x in sub)
        out.append(f"\n## {name} ({arm} vs {base}): files up {n['up']} / flat {n['flat']} / down {n['down']} of {len(sub)};"
                   f" sum {total_b:.1f} -> {total_a:.1f} ({100 * (total_a / total_b - 1):+.2f}%)")
        for _, g in GROUPS:
            gs = [x for x in sub if x["group"] == g]
            if not gs:
                continue
            b = sum(x["base"] for x in gs) / len(gs)
            a = sum(x["mean"] for x in gs) / len(gs)
            per = " ".join(f"{x['delta_pct']:+.1f}({x['z']:+.1f})" for x in gs)
            out.append(f"  {g:22s} {len(gs):2d} {b:9.1f} {a:9.1f} {100 * (a / b - 1):+7.2f}%  {per}")
    with (ROOT / "pd_per_file.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        for x in rows:
            w.writerow({k: (round(v, 2) if isinstance(v, float) else v) for k, v in x.items()})
    return out


def identity():
    out = ["\n## Identity checks (native TSV and PIN SHA-256)"]
    for arm, ref, what in IDENTITY:
        files = [r["id"] for r in MANIFEST if es.load(r["id"], arm) and es.load(r["id"], ref)]
        same = [d for d in files if es.compare(d, arm, ref)["identical"]]
        out.append(f"  {'PASS' if files and len(same) == len(files) else 'FAIL'} {arm} == {ref} on {len(same)}/{len(files)} files ({what}): {' '.join(files)}")
    return out


def velos():
    """evaluate_split.velos() writes split_entrapment_velos.tsv; keep that file and move this table aside."""
    kept = ROOT / "split_entrapment_velos.tsv"
    previous = kept.read_bytes() if kept.exists() else None
    es.VELOS_ARMS = VELOS_ARMS
    text = es.velos()
    kept.rename(ROOT / "pd_entrapment_velos.tsv")
    if previous is not None:
        kept.write_bytes(previous)
    return text


def shuffled():
    out = ["\n## Shuffled paired-mass entrapment (r = 1; combined FDP = 2 N_E / N, estimated true PSMs N - 2 N_E;"
           " seed means summed over files)"]
    for g, files in SHUFFLED.items():
        for arm in ("a35_raw_local_E", "re_auto_E", "pd_on_E", "pd_on_E-nort", "pd_on_E-noms2",
                    "pd_cal_E", "pd_cal_E-nort", "pd_cal_E-noms2"):
            c = [es.counts(d, arm) for d in files]
            if any(x is None for x in c):
                out.append(f"  {g:20s} {arm:16s} missing")
                continue
            n, e = sum(x[0] for x in c), sum(x[1] for x in c)
            out.append(f"  {g:20s} {arm:16s} N {n:8.1f}  N_E {e:5.1f}  FDP {200 * e / n:.2f}%  true {n - 2 * e:8.1f}")
    return out


def cost():
    out = ["\n## Search cost per file (wall s, CPU s, peak RSS MiB; median over files, OMP_NUM_THREADS=4, 2 searches in parallel)"]
    rows = []
    for r in MANIFEST:
        for arm in ("re_auto", "pd_on", "pd_cal", "pd_inst"):
            s = es.load(r["id"], arm)
            if s:
                rows.append(dict(dataset=r["id"], group=group_of(r["id"]), arm=arm, wall_s=round(s["search"]["seconds"], 1),
                                 cpu_s=round(s["search"]["user_seconds"] + s["search"]["system_seconds"], 1),
                                 rss_mib=round(s["search"]["max_rss_kib"] / 1024)))
    for arm in ("re_auto", "pd_on", "pd_cal", "pd_inst"):
        sub = [x for x in rows if x["arm"] == arm]
        if sub:
            out.append(f"  {arm:8s} n={len(sub):2d}  wall {statistics.median(x['wall_s'] for x in sub):7.1f}"
                       f"  CPU {statistics.median(x['cpu_s'] for x in sub):7.1f}  RSS {statistics.median(x['rss_mib'] for x in sub):6.0f}")
    paired = {}
    for x in rows:
        paired.setdefault(x["dataset"], {})[x["arm"]] = x
    for arm in ("pd_on", "pd_cal"):
        ratios = [(d, v[arm]["cpu_s"] / v["re_auto"]["cpu_s"], v[arm]["wall_s"] / v["re_auto"]["wall_s"],
                   v[arm]["rss_mib"] - v["re_auto"]["rss_mib"]) for d, v in paired.items() if arm in v and "re_auto" in v]
        if ratios:
            out.append(f"  {arm} / re_auto: CPU x{statistics.median(r[1] for r in ratios):.2f} (range {min(r[1] for r in ratios):.2f}-{max(r[1] for r in ratios):.2f}),"
                       f" wall x{statistics.median(r[2] for r in ratios):.2f}, RSS {statistics.median(r[3] for r in ratios):+.0f} MiB (median)")
    out.append("  pd_inst ran while Percolator-only ablations and diagnostics used the same 4 cores; its cost is not comparable.")
    with (ROOT / "pd_cost.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["dataset", "group", "arm", "wall_s", "cpu_s", "rss_mib"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return out


def calibration():
    """NCE chosen and retention-time calibration, from each pd_on / pd_inst search log."""
    out = ["\n## PeptDeep calibration per file (from search.log)"]
    for arm in ("pd_on", "pd_cal", "pd_inst"):
        for r in MANIFEST:
            log = ROOT / "results" / r["id"] / arm / "search.log"
            if not log.exists():
                continue
            lines = [l.split("] ", 1)[1] for l in log.read_text().splitlines() if l.startswith("[PeptDeepRescoring]")]
            out.append(f"  {arm:8s} {r['id']:22s} " + " | ".join(re.sub(r"run '[^']*': ", "", l) for l in lines))
    return out


def main():
    text = "\n".join(yields() + identity() + velos() + shuffled() + cost() + calibration()).lstrip("\n")
    (ROOT / "pd_summary.txt").write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
