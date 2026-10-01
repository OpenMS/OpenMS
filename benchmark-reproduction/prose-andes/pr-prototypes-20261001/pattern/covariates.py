"""Per-file covariates against per-file effects (Spearman rank correlations across files)."""
import csv, json, statistics
from pathlib import Path
from evaluate import MANIFEST, group_of
ROOT = Path(__file__).resolve().parent.parent
PROV = ROOT.parent / "prose-andes-reproduction/broad_benchmark/provenance"
L = lambda d, a: json.loads((ROOT / "results" / d / a / "summary.json").read_text())


def qcount(d, arm, q):
    n = []
    for s in (1, 42, 137):
        with open(ROOT / "results" / d / arm / f"s{s}" / "target.tsv") as f:
            n.append(sum(float(r["q-value"]) <= q for r in csv.DictReader(f, delimiter="\t")))
    return statistics.mean(n)


def rank(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i]); r = [0] * len(xs)
    for k, i in enumerate(order): r[i] = k
    return r


def spearman(a, b):
    ra, rb = rank(a), rank(b); n = len(a)
    ma, mb = sum(ra) / n, sum(rb) / n
    cov = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    return cov / (sum((x - ma) ** 2 for x in ra) * sum((y - mb) ** 2 for y in rb)) ** 0.5


rows = []
for r in MANIFEST:
    d = r["id"]; b = L(d, "base")
    prov = json.loads((PROV / f"{d}.json").read_text())
    ch = prov["source_charges"]; tot = sum(ch.values())
    eff = lambda a: 100 * (L(d, a)["mean_psms"] / b["mean_psms"] - 1)
    rows.append(dict(dataset=d, group=group_of(d), unit=r["fragment_unit"], id_rate=round(100 * b["mean_psms"] / 8000, 1),
                     marginal=round(qcount(d, "base", 0.05) / qcount(d, "base", 0.01) - 1, 3),
                     z3plus=round(100 * sum(v for k, v in ch.items() if int(k) >= 3) / tot, 1),
                     priors=round(eff("f78x_priors"), 2), pr10335=round(eff("pr10335"), 2),
                     frag_window=round(eff("h_frag"), 2) if r["fragment_unit"] == "ppm" else "",
                     mass_scorer=round(eff("pr65s_mass"), 2) if r["fragment_unit"] == "ppm" else ""))
with open(ROOT / "pattern" / "covariates.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
for x in rows: print("\t".join(str(v) for v in x.values()))
print()
for eff in ["priors", "pr10335", "frag_window", "mass_scorer"]:
    sub = [x for x in rows if x[eff] != ""]
    for cov in ["id_rate", "marginal", "z3plus"]:
        print(f"{eff:12s} vs {cov:9s} rho={spearman([x[cov] for x in sub], [x[eff] for x in sub]):+.2f} (n={len(sub)})")
