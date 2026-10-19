"""Evaluate one arm directory (search.pin.tsv + s<seed>/target.tsv from run_search.sh).

Per seed: target PSMs at q <= 0.01 (entrapment PSMs count as targets), split by chimeric_rank (from the PIN; rank 1
when the column is absent), unique peptides (modified sequence, flanks removed), entrapment PSMs (all proteins
ENTRAP_), combined entrapment FDP 2*N_E/N (r = 1), and PSMs at an entrapment FDP of 1% (largest q-ordered prefix with
2*N_E/N <= 0.01, the measure andes reports). Prints one JSON line per seed and the seed means.
usage: evaluate.py <arm dir> [<arm dir> ...]"""
import csv, json, sys
from pathlib import Path

SEEDS = (1, 42, 137)


def pin_ranks(pin):
    ranks = {}
    with open(pin) as f:
        header = f.readline().rstrip("\n").split("\t")
        i_id = header.index("SpecId") if "SpecId" in header else header.index("PSMId")
        i_r = header.index("chimeric_rank") if "chimeric_rank" in header else None
        for line in f:
            row = line.rstrip("\n").split("\t")
            if row[0].lower() == "defaultdirection": continue
            ranks[row[i_id]] = int(float(row[i_r])) if i_r is not None else 1
    return ranks


def strip(pep):
    return pep[2:-2] if len(pep) > 4 and pep[1] == "." and pep[-2] == "." else pep


def seed_stats(arm, seed, ranks):
    rows = []
    with open(arm / f"s{seed}" / "target.tsv") as f:
        header = f.readline().rstrip("\n").split("\t")
        iq, ip, iprot = header.index("q-value"), header.index("peptide"), header.index("proteinIds")
        for line in f:
            r = line.rstrip("\n").split("\t")
            prots = r[iprot:]
            rows.append((float(r[iq]), r[0], strip(r[ip]), all(p.startswith("ENTRAP_") for p in prots if p)))
    rows.sort(key=lambda x: x[0])  # Percolator's order: by q, then score
    acc = [x for x in rows if x[0] <= 0.01]
    out = {"seed": seed}
    for label, sel in (("all", lambda x: True), ("r1", lambda x: ranks.get(x[1], 1) == 1), ("r2", lambda x: ranks.get(x[1], 1) == 2)):
        a = [x for x in acc if sel(x)]
        ne = sum(x[3] for x in a)
        out[label] = {"psms": len(a), "peptides": len({x[2] for x in a}), "entrap": ne,
                      "fdp": round(2 * ne / len(a), 5) if a else None}
    # rank-2 PSMs at q <= 0.01: same peptide as an accepted rank-1 PSM of the same spectrum; peptides no rank-1 PSM has
    r1 = [x for x in acc if ranks.get(x[1], 1) == 1]
    r1_scan_pep = {(x[1].split("#")[0], x[2]) for x in r1}
    r1_peps = {x[2] for x in r1}
    r2 = [x for x in acc if ranks.get(x[1], 1) == 2]
    out["r2_same_as_r1_in_scan"] = sum(1 for x in r2 if (x[1].split("#")[0], x[2]) in r1_scan_pep)
    out["r2_new_peptides"] = len({x[2] for x in r2} - r1_peps)
    # PSMs at an entrapment FDP of 1% (over the q-ordered list)
    best, ne = 0, 0
    for k, x in enumerate(rows, 1):
        ne += x[3]
        if 2 * ne / k <= 0.01: best = k
    out["psms_at_fdp1"] = best
    out["peptides_at_fdp1"] = len({x[2] for x in rows[:best]})
    out["r2_at_fdp1"] = sum(1 for x in rows[:best] if ranks.get(x[1], 1) == 2)
    return out


def main():
    for d in sys.argv[1:]:
        arm = Path(d)
        ranks = pin_ranks(arm / "input.pin")
        per = [seed_stats(arm, s, ranks) for s in SEEDS]
        for p in per: print(json.dumps({"arm": str(arm), **p}))
        mean = lambda f: round(sum(f(p) for p in per) / len(per), 2)
        summ = {"arm": str(arm), "pin_rows": len(ranks), "pin_rank2_rows": sum(1 for v in ranks.values() if v == 2)}
        for lab in ("all", "r1", "r2"):
            summ[lab] = {k: mean(lambda p: p[lab][k] or 0) for k in ("psms", "peptides", "entrap")}
            n, e = summ[lab]["psms"], summ[lab]["entrap"]
            summ[lab]["fdp"] = round(2 * e / n, 5) if n else None
        for k in ("psms_at_fdp1", "peptides_at_fdp1", "r2_at_fdp1", "r2_same_as_r1_in_scan", "r2_new_peptides"): summ[k] = mean(lambda p: p[k])
        summ["sd_psms"] = round((sum((p["all"]["psms"] - summ["all"]["psms"]) ** 2 for p in per) / (len(per) - 1)) ** 0.5, 2)
        print("MEAN " + json.dumps(summ))


main()
