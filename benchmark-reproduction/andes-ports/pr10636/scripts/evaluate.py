"""Aggregate the class-FDR benchmark.

Input: <root>/<setting>/<dataset>/s<seed>/<arm>.tsv.gz from extract.py (arm: global, or the class definition used
with -PSM_class). Accepted = best hit per spectrum, not decoy, q <= 0.01. Entrapment FDP (combined, r = 1) = 2 N_E / N.
Per file: mean over seeds 1/42/137; per group: N and N_E are seed means summed over the group's files (as in
open-prs-on-fix-20261002/evaluate_entrap_d.py). A file is up/down when |z| >= 2, z = delta / sqrt(var_a/3 + var_b/3).
usage: evaluate.py <root> <out_dir>
"""
import csv, gzip, json, math, statistics, sys
from collections import defaultdict
from pathlib import Path

SEEDS = (1, 42, 137)
GROUPS = {"velos": "Velos CID", "hfx": "HF-X HCD", "astral": "Astral HCD", "lumos_lfq": "Lumos HCD LFQ",
          "lumos_tmt": "Lumos CID TMT", "eclipse": "Exploris 480 TMTpro", "tims": "timsTOF HT", "liver": "timsTOF Pro HLA-I"}
COL = {"mod": "modclass", "bench_mc": "mc", "bench_z": "z", "bench_iso": "iso", "bench_term": "term", "bench_len": "len"}
ALPHA = 0.01

def group(d):
    for k, v in GROUPS.items():
        if d.startswith(k): return v
    return d

def load(p):
    with gzip.open(p, "rt") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def classify(r, cdef):
    v = r[COL[cdef]]
    if cdef == "mod" and v != "unmodified":
        return v
    return v

def stats(rows, cdef):
    """per class: accepted N, entrapment N_E, accepted decoys D, and all-q targets/decoys"""
    out = defaultdict(lambda: [0, 0, 0])
    for r in rows:
        if float(r["q"]) > ALPHA: continue
        c = classify(r, cdef)
        for key in (c, "ALL"):
            if r["decoy"] == "1": out[key][2] += 1
            else:
                out[key][0] += 1
                out[key][1] += int(r["entrap"])
    return out

def main():
    root, outdir = Path(sys.argv[1]), Path(sys.argv[2]); outdir.mkdir(parents=True, exist_ok=True)
    # data[(setting, dataset, cdef)][arm][seed] = {class: [N, NE, D]}
    data = defaultdict(lambda: defaultdict(dict))
    for setting in sorted(p for p in root.iterdir() if p.is_dir()):
        for ds in sorted(p for p in setting.iterdir() if p.is_dir()):
            for s in SEEDS:
                sd = ds / f"s{s}"
                g = sd / "global.tsv.gz"
                if not g.exists(): continue
                grows = load(g)
                for cfile in sorted(sd.glob("*.tsv.gz")):
                    cdef = cfile.name[:-7]
                    if cdef == "global": continue
                    data[(setting.name, ds.name, cdef)]["global"][s] = stats(grows, cdef)
                    data[(setting.name, ds.name, cdef)]["class"][s] = stats(load(cfile), cdef)
    per_file = []
    for (setting, ds, cdef), arms in sorted(data.items()):
        if len(arms["global"]) != 3 or len(arms["class"]) != 3:
            print("incomplete", setting, ds, cdef, file=sys.stderr); continue
        classes = sorted(set().union(*[set(arms[a][s]) for a in arms for s in SEEDS]))
        for c in classes:
            row = dict(setting=setting, dataset=ds, group=group(ds), class_def=cdef, cls=c)
            for a in ("global", "class"):
                for i, name in enumerate(("N", "NE", "D")):
                    vals = [arms[a][s].get(c, [0, 0, 0])[i] for s in SEEDS]
                    row[f"{a}_{name}"] = statistics.mean(vals)
                    if name == "N": row[f"{a}_N_var"] = statistics.variance(vals)
            se = math.sqrt(row["global_N_var"] / 3 + row["class_N_var"] / 3)
            d = row["class_N"] - row["global_N"]
            row["z"] = (d / se) if se > 0 else (0.0 if d == 0 else math.copysign(math.inf, d))
            per_file.append(row)
    keys = list(per_file[0].keys()) if per_file else []
    with open(outdir / "per_file.tsv", "w") as f:
        w = csv.DictWriter(f, keys, delimiter="\t", lineterminator="\n"); w.writeheader()
        for r in per_file: w.writerow({k: (f"{v:.3f}" if isinstance(v, float) else v) for k, v in r.items()})
    # per group and suite
    agg = defaultdict(lambda: defaultdict(float)); files = defaultdict(lambda: [0, 0, 0])
    for r in per_file:
        for gname in (r["group"], "ALL files"):
            k = (r["setting"], r["class_def"], gname, r["cls"])
            for a in ("global", "class"):
                for n in ("N", "NE", "D"): agg[k][f"{a}_{n}"] += r[f"{a}_{n}"]
            if r["cls"] == "ALL" and gname != "ALL files" or (gname == "ALL files" and r["cls"] == "ALL"):
                pass
            if r["cls"] == "ALL":
                fz = files[k]
                fz[0 if r["z"] >= 2 else (2 if r["z"] <= -2 else 1)] += 1
    with open(outdir / "per_group.tsv", "w") as f:
        f.write("setting\tclass_def\tgroup\tclass\tglobal_N\tglobal_NE\tglobal_FDP_pct\tglobal_decoyFDR_pct\tclass_N\tclass_NE\tclass_FDP_pct\tclass_decoyFDR_pct\tdelta_N_pct\tfiles_up_flat_down\n")
        for k in sorted(agg):
            v = agg[k]
            fdp = lambda a: 200 * v[f"{a}_NE"] / v[f"{a}_N"] if v[f"{a}_N"] else float("nan")
            dfdr = lambda a: 100 * v[f"{a}_D"] / v[f"{a}_N"] if v[f"{a}_N"] else float("nan")
            dn = 100 * (v["class_N"] - v["global_N"]) / v["global_N"] if v["global_N"] else float("nan")
            fz = files.get(k)
            f.write("\t".join(map(str, [*k, f"{v['global_N']:.1f}", f"{v['global_NE']:.1f}", f"{fdp('global'):.2f}",
                                        f"{dfdr('global'):.2f}", f"{v['class_N']:.1f}", f"{v['class_NE']:.1f}",
                                        f"{fdp('class'):.2f}", f"{dfdr('class'):.2f}", f"{dn:+.2f}",
                                        "/".join(map(str, fz)) if fz else ""])) + "\n")
    print((outdir / "per_group.tsv").read_text())

if __name__ == "__main__":
    main()
