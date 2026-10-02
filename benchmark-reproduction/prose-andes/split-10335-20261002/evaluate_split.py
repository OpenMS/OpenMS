"""#10335 split into #10394 (deduplication) and #10397 (multiply charged fragments), 2026-10-02.

Comparisons (arm against baseline on the same files; per file mean Percolator PSMs over seeds 1/42/137,
z = diff / SE with SE = sqrt(var_arm / 3 + var_base / 3); up z >= 2, down z <= -2, else flat):
  #10394:  d_dedup (develop 3a47278 + #10394, = develop 9517361) against d_base (develop 3a47278).
  #10335:  d_35 (#10335 136dde3 merged with 3a47278) against d_dedup.
  #10397:  d_fc (low-resolution files, hfx_A2, astral_A2) and d_fc_hr (the other high-resolution files) against d_dedup.
  Ablation of #10335 on the six low-resolution files (its build, scoring:method x scoring:fragment_charges)
  and its opt-in options, each against d_35 (#10335's defaults) and against d_dedup.
Entrapment: natural UPS1/yeast on Velos (entrapment_velos.py; combined FDP with r = 3.38, summed over
3 files x 3 seeds), shuffled paired-mass (r = 1) on HF-X, Astral and TMTpro (evaluate_entrap_d.py).
"""
import csv
import json
import math
import statistics
from pathlib import Path

import entrapment_velos as ev
from evaluate import GROUPS, MANIFEST, group_of

ROOT = Path(__file__).resolve().parent
LOW_RES = [r["id"] for r in MANIFEST if r["fragment_unit"] == "Da"]
COMPARISONS = [  # name, arm, baseline
    ("#10394 deduplication", "d_dedup", "d_base"),
    ("#10335 on develop with #10394", "d_35", "d_dedup"),
    ("#10397 multiply charged fragments", "d_fc+d_fc_hr", "d_dedup"),
    ("#10335 on develop with #10397 (2f30b40)", "d_35", "d_fc"),
    ("#10335 build: HyperScore, single charge", "a35_hs_single", "d_dedup"),
    ("#10335 build: HyperScore, multiple charges", "a35_hs_multi", "d_dedup"),
    ("#10335 build: calibrated, single charge", "a35_cal_single", "d_dedup"),
    ("#10335 build: calibrated, multiple charges (= #10335 default)", "d_35", "d_dedup"),
    ("#10335 opt-in: fragment:query_spectrum=raw", "a35_raw", "d_35"),
    ("#10335 opt-in: raw + annotate:local_fragment_evidence", "a35_raw_local", "d_35"),
    ("#10335 opt-in: annotate:local_fragment_evidence", "a35_local", "d_35"),
    ("#10335 opt-in: local evidence, high-resolution files", "a35_local_hr", "d_35"),
    ("port: raw retrieval, low-resolution files", "re_r", "d_fc"),
    ("port: raw + local evidence", "re_rl", "d_fc+d_fc_hr"),
    ("PR defaults (query_spectrum=auto, local evidence on)", "re_auto", "d_fc+d_fc_hr"),
]
# arms whose outputs must equal another arm byte for byte (native TSV and PIN)
IDENTITY = [("d_dedup_off", "d_base", "peptide:deduplicate=false restores develop 3a47278"),
            ("d_fc+d_fc_hr", "d_dedup", "#10397 leaves deisotoped (high-resolution) spectra unchanged"),
            ("d_fc", "a35_hs_multi", "#10397 equals #10335 with HyperScore and multiple charges"),
            ("d_35", "d_dedup", "#10335 equals develop with #10394 on deisotoped spectra"),
            ("re_off", "d_fc+d_fc_hr", "the port with both options off equals develop"),
            ("re_rl", "a35_raw_local", "the port reproduces #10335's raw + local arm on high-resolution files"),
            ("re_auto", "re_rl", "the PR defaults equal raw + local on high-resolution files")]
VELOS_ARMS = ["d_base", "d_dedup", "d_35", "d_fc", "a35_hs_single", "a35_hs_multi", "a35_cal_single", "a35_local", "re_r", "re_rl", "re_auto"]


def load(d, a):
    p = ROOT / "results" / d / a / "summary.json"
    return json.loads(p.read_text()) if p.exists() else None


def compare(d, arm, base):
    arm = next((a for a in arm.split("+") if load(d, a)), arm)
    base = next((a for a in base.split("+") if load(d, a)), base)
    s, b = load(d, arm), load(d, base)
    if s is None or b is None:
        return None
    va, vb = statistics.variance(s["seed_psms"]), statistics.variance(b["seed_psms"])
    se = math.sqrt(va / 3 + vb / 3)
    diff = s["mean_psms"] - b["mean_psms"]
    z = diff / se if se else (0.0 if diff == 0 else math.copysign(math.inf, diff))
    return dict(base=b["mean_psms"], mean=s["mean_psms"], delta_pct=100 * diff / b["mean_psms"], z=z,
                verdict="up" if z >= 2 else "down" if z <= -2 else "flat",
                native_delta=s["native"]["psms"] - b["native"]["psms"],
                identical=s["native_tsv_sha256"] == b["native_tsv_sha256"] and s["native_pin_sha256"] == b["native_pin_sha256"])


def yields():
    rows = []
    for name, arm, base in COMPARISONS:
        for r in MANIFEST:
            c = compare(r["id"], arm, base)
            if c:
                rows.append(dict(comparison=name, arm=arm, baseline=base, dataset=r["id"], group=group_of(r["id"]), **c))
    with (ROOT / "split_per_file.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        for x in rows:
            w.writerow({k: (round(v, 2) if isinstance(v, float) else v) for k, v in x.items()})
    out = []
    for name, arm, base in COMPARISONS:
        sub = [x for x in rows if x["comparison"] == name and x["arm"] == arm]
        if not sub:
            out.append(f"\n## {name} ({arm} vs {base}): no results")
            continue
        n = {v: sum(x["verdict"] == v for x in sub) for v in ("up", "flat", "down")}
        out.append(f"\n## {name} ({arm} vs {base}): files up {n['up']} / flat {n['flat']} / down {n['down']} of {len(sub)}")
        for _, g in GROUPS:
            gs = [x for x in sub if x["group"] == g]
            if not gs:
                continue
            b = sum(x["base"] for x in gs) / len(gs)
            a = sum(x["mean"] for x in gs) / len(gs)
            per = " ".join(f"{x['delta_pct']:+.1f}({x['z']:+.1f})" for x in gs)
            ident = " identical" if all(x["identical"] for x in gs) else ""
            out.append(f"  {g:22s} {len(gs):2d} {b:9.1f} {a:9.1f} {100 * (a / b - 1):+7.2f}%  {per}{ident}")
    return out


def identity():
    out = ["\n## Identity checks (native TSV and PIN SHA-256)"]
    for arm, ref, what in IDENTITY:
        files = [r["id"] for r in MANIFEST if any(load(r["id"], a) for a in arm.split("+")) and any(load(r["id"], a) for a in ref.split("+"))]
        if ref == "d_dedup" or "high-resolution" in what:
            files = [d for d in files if d not in LOW_RES]
        same = [d for d in files if compare(d, arm, ref)["identical"]]
        out.append(f"  {'PASS' if files and len(same) == len(files) else 'FAIL'} {arm} == {ref} on {len(same)}/{len(files)} files ({what}): {' '.join(files)}")
    return out


def velos():
    entries = ev.read_fasta()
    db, _ = ev.classes(entries)
    sample_peps, ent_peps = set(), set()
    for h, s in entries:
        s2 = s.replace("I", "L")
        acc = h.split("|")[1] if "|" in h else ""
        if "_HUMAN" in h.split()[0] and "Cont_" not in h and acc not in ev.UPS1:
            ent_peps |= ev.digest(s2)
        else:
            sample_peps |= ev.digest(s2)
    r = len(ent_peps - sample_peps) / len(sample_peps)
    cache = {}

    def category(pep):
        if pep not in cache:
            cache[pep] = "sample" if pep in db["sample"] else "ups1" if pep in db["ups1"] else "entrapment" if pep in db["human"] else "unknown"
        return cache[pep]

    out = [f"\n## Velos natural entrapment (r = {r:.3f}; combined FDP = N_E (1 + 1/r) / N over 3 files x 3 seeds)",
           f"  {'arm':16s} {'N':>6s} {'N_E':>4s} {'FDP':>7s}"]
    rows = []
    for arm in VELOS_ARMS:
        n = e = 0
        for d in ev.FILES:
            for seed in (1, 42, 137):
                p = ROOT / "results" / d / arm / f"s{seed}" / "target.tsv"
                if not p.exists():
                    break
                with p.open() as f:
                    for row in csv.DictReader(f, delimiter="\t"):
                        if float(row["q-value"]) <= 0.01:
                            n += 1
                            e += category(ev.stripped(row["peptide"])) == "entrapment"
            else:
                continue
            n = 0
            break
        if n:
            fdp = 100 * e * (1 + 1 / r) / n
            rows.append(dict(arm=arm, psms=n, entrapment=e, fdp_combined=round(fdp, 3)))
            out.append(f"  {arm:16s} {n:6d} {e:4d} {fdp:6.3f}%")
        else:
            out.append(f"  {arm:16s} incomplete")
    with (ROOT / "split_entrapment_velos.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["arm", "psms", "entrapment", "fdp_combined"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return out


def counts(d, arm):
    """Seed means of accepted target PSMs and of entrapment hits (all proteins ENTRAP_*), as evaluate_entrap_d.py."""
    n, e = [], []
    for s in (1, 42, 137):
        p = ROOT / "results" / d / arm / f"s{s}" / "target.tsv"
        if not p.exists():
            return None
        acc = [l.split("\t") for l in p.read_text().splitlines()[1:]]
        acc = [f for f in acc if float(f[2]) <= 0.01]
        n.append(len(acc))
        e.append(sum(all(x.startswith("ENTRAP_") for x in f[5:] if x) for f in acc))
    return statistics.mean(n), statistics.mean(e)


def shuffled():
    EG = {"HF-X HCD": ["hfx_A2", "hfx_B1", "hfx_B3"], "Astral HCD": ["astral_A2", "astral_B1", "astral_B3"],
          "Exploris 480 TMTpro": ["eclipse_tmtpro_10855", "eclipse_tmtpro_10858", "eclipse_tmtpro_10863"]}
    out =["\n## Shuffled paired-mass entrapment (r = 1; combined FDP = 2 N_E / N, seed means summed over files)"]
    for g, files in EG.items():
        for name, arm in (("develop 3a47278", lambda d: "tds1_E" if d.startswith("eclipse") else "hds1_E"),
                          ("#10394 (= 2f30b40)", lambda d: "d_dedup_E"),
                          ("+ local evidence", lambda d: "a35_local_E"),
                          ("+ raw + local", lambda d: "a35_raw_local_E")):
            c = [counts(d, arm(d)) for d in files]
            if any(x is None for x in c):
                out.append(f"  {g:20s} {name:16s} missing")
                continue
            n, e = sum(x[0] for x in c), sum(x[1] for x in c)
            out.append(f"  {g:20s} {name:16s} N {n:8.1f}  N_E {e:5.1f}  FDP {200 * e / n:.2f}%")
    return out


def main():
    text = "\n".join(yields() + identity() + velos() + shuffled()).lstrip("\n")
    (ROOT / "split_summary.txt").write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
