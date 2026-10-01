"""Effects and data properties by analyzer class (ion trap / Orbitrap / Astral / timsTOF) and by group."""
import csv, json, statistics
from pathlib import Path
from evaluate import MANIFEST, group_of
ROOT = Path(__file__).resolve().parent.parent
BENCH = json.loads((ROOT.parent / "prose-andes-reproduction/broad_benchmark/summary.json").read_text())["yield_results"]
ANDES = {(r["dataset"], r["arm"]): r["mean_psms"] for r in BENCH if r["arm"].startswith("andes")}
tsv = lambda name: {r["dataset"]: r for r in csv.DictReader(open(ROOT / "pattern" / name), delimiter="\t")}
COV, SPEC, TAIL = tsv("covariates.tsv"), tsv("spectra.tsv"), tsv("fragment_tails.tsv")
CLASS = {"Velos CID": "ion trap CID", "Lumos CID TMT": "ion trap CID", "HF-X HCD": "Orbitrap HCD",
         "Lumos HCD LFQ": "Orbitrap HCD", "Exploris 480 TMTpro": "Orbitrap HCD", "Astral HCD": "Astral (TOF)",
         "timsTOF HT": "timsTOF (QTOF)"}


def psms(d, arms):
    for a in arms:
        p = ROOT / "results" / d / a / "summary.json"
        if p.exists(): return json.loads(p.read_text())["mean_psms"]


rows = []
for r in MANIFEST:
    d = r["id"]; g = group_of(d); base = psms(d, ["base"])
    eff = lambda arms: None if psms(d, arms) is None else 100 * (psms(d, arms) / base - 1)
    t = TAIL.get(d)
    rows.append(dict(
        dataset=d, group=g, cls=CLASS[g], per100=float(SPEC[d]["per100"]), top20_tic=float(SPEC[d]["top20_tic_pct"]),
        id_rate=float(COV[d]["id_rate"]), marginal=100 * float(COV[d]["marginal"]), z3=float(COV[d]["z3plus"]),
        sd_weak=float(t["sd_weak"]) if t else None,
        cut_ratio=float(t["decoy_beyond_w_pct"]) / float(t["target_beyond_w_pct"]) if t and float(t["target_beyond_w_pct"]) else None,
        priors=eff(["f78x_priors"]), pr10335=eff(["pr10335"]),
        f79=eff(["f79_default", "h_frag"]) if r["fragment_unit"] == "ppm" else None,
        f79_mass=eff(["f79_mass", "m_shift_frag"]) if r["fragment_unit"] == "ppm" else None,
        wt40=eff(["wt40"]) if r["fragment_unit"] == "ppm" else None, wt100=eff(["wt100"]) if r["fragment_unit"] == "ppm" else None,
        andes_gap=100 * (ANDES[(d, "andes")] / base - 1),
        andes_highres_gap=100 * (ANDES[(d, "andes_highres")] / base - 1) if (d, "andes_highres") in ANDES else None))

KEYS = ["per100", "top20_tic", "id_rate", "marginal", "z3", "sd_weak", "cut_ratio", "priors", "pr10335", "f79", "f79_mass",
        "wt40", "wt100", "andes_gap", "andes_highres_gap"]


def agg(sub, by):
    out = []
    for k in dict.fromkeys(x[by] for x in sub):
        grp = [x for x in sub if x[by] == k]
        row = {by: k, "n": len(grp)}
        for key in KEYS:
            v = [x[key] for x in grp if x[key] is not None]
            row[key] = round(statistics.mean(v), 2) if v else ""
        out.append(row)
    return out


with open(ROOT / "pattern" / "classes.tsv", "w", newline="") as f:
    for by in ["cls", "group"]:
        a = agg(rows, by)
        w = csv.DictWriter(f, fieldnames=list(a[0]), delimiter="\t"); w.writeheader(); w.writerows(a); f.write("\n")
        print("\t".join(a[0])); [print("\t".join(str(v) for v in x.values())) for x in a]; print()
