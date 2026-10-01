"""Entrapment check on the Velos UPS1/yeast files (PXD001819).

The samples hold yeast plus the 48 human UPS1 proteins; the database (hy) holds all reviewed human and
yeast proteins plus contaminants. A PSM whose peptide (I = L) occurs only in human proteins other than
UPS1 is false by construction. FDP estimates per Wen et al. (2025): lower bound N_E / N and combined
N_E * (1 + 1/r) / N, with r = entrapment-only / sample peptides of the searchable tryptic space.
"""
import csv
import json
import re
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
FASTA = ROOT.parent / "prose-andes-reproduction/broad_benchmark/data/hy_td.fasta"
# Sigma UPS1 (48 proteins; ubiquitin P62988 is listed under the current polyubiquitin entries)
UPS1 = """P00915 P00918 P01031 P69905 P68871 P41159 P02768 P0CG48 P0CG47 P62987 P62979 P04040 P00167 P01133 P02144 P15559
P62937 Q06830 P63165 P00709 P06730 P12081 P61626 Q15843 P02753 P16083 P63279 P01008 P61769 P55957 O76070 P08263 P01344
P01127 P10599 P99999 P06396 P09211 P01112 P01579 P02787 P02788 P05413 P10636 P10145 P02741 O00762 P01375 P51965 P08758""".split()
FILES = ["velos_125_R1", "velos_5000_R2", "velos_25000_R3"]


def read_fasta():
    entries, hdr, seq = [], None, []
    for line in FASTA.open():
        if line.startswith(">"):
            if hdr: entries.append((hdr, "".join(seq)))
            hdr, seq = line[1:].strip(), []
        else:
            seq.append(line.strip())
    if hdr: entries.append((hdr, "".join(seq)))
    return [(h, s) for h, s in entries if not h.startswith("DECOY_")]


def classes(entries):
    groups = {"sample": [], "ups1": [], "human": []}
    found = set()
    for h, s in entries:
        acc = h.split("|")[1] if "|" in h else h.split()[0]
        s = s.replace("I", "L")
        if "Cont_" in h or "_YEAST" in h.split()[0]:
            groups["sample"].append(s)
        elif acc in UPS1:
            groups["ups1"].append(s); found.add(acc)
        elif "_HUMAN" in h.split()[0]:
            groups["human"].append(s)
        else:
            groups["sample"].append(s)  # other species (contaminant-like)
    missing = sorted(set(UPS1) - found)
    return {k: "\n".join(v) for k, v in groups.items()}, missing


def digest(s):
    sites = [0] + [i + 1 for i in range(len(s) - 1) if s[i] in "KR" and s[i + 1] != "P"] + [len(s)]
    out = set()
    for i in range(len(sites) - 1):
        for j in range(i + 1, min(i + 4, len(sites))):
            p = s[sites[i]:sites[j]]
            if 7 <= len(p) <= 40: out.add(p)
    return out


def stripped(peptide):
    plain = re.sub(r"\[[^\]]*\]", "", peptide)  # modification masses contain dots
    core = plain.split(".")[1] if plain.count(".") == 2 else plain
    return re.sub(r"[^A-Z]", "", core).replace("I", "L")


def main(arms):
    entries = read_fasta()
    db, missing = classes(entries)
    sample_peps, ent_peps = set(), set()
    for h, s in entries:
        s2 = s.replace("I", "L")
        acc = h.split("|")[1] if "|" in h else ""
        if "_HUMAN" in h.split()[0] and "Cont_" not in h and acc not in UPS1:
            ent_peps |= digest(s2)
        else:
            sample_peps |= digest(s2)
    r = len(ent_peps - sample_peps) / len(sample_peps)
    print(f"UPS1 accessions not in FASTA: {missing}; r = {r:.3f}", file=sys.stderr)
    cache = {}

    def category(pep):
        if pep not in cache:
            cache[pep] = "sample" if pep in db["sample"] else "ups1" if pep in db["ups1"] else "entrapment" if pep in db["human"] else "unknown"
        return cache[pep]

    accepted = {}
    rows = []
    for arm in arms:
        for d in FILES:
            for seed in (1, 42, 137):
                path = ROOT / "results" / d / arm / f"s{seed}" / "target.tsv"
                if not path.exists():
                    continue
                acc = {}
                with path.open() as f:
                    for row in csv.DictReader(f, delimiter="\t"):
                        if float(row["q-value"]) <= 0.01:
                            acc[row["PSMId"].rsplit("_", 1)[0]] = stripped(row["peptide"])
                accepted[(arm, d, seed)] = acc
    for arm in arms:
        tot = {"n": 0, "sample": 0, "ups1": 0, "entrapment": 0, "unknown": 0}
        gained = {"n": 0, "entrapment": 0}
        for d in FILES:
            for seed in (1, 42, 137):
                acc = accepted.get((arm, d, seed))
                if acc is None:
                    continue
                for scan, pep in acc.items():
                    tot["n"] += 1; tot[category(pep)] += 1
                base = accepted.get(("base", d, seed), {})
                for scan, pep in acc.items():
                    if scan not in base:
                        gained["n"] += 1; gained["entrapment"] += category(pep) == "entrapment"
        if not tot["n"]:
            continue
        lb = tot["entrapment"] / tot["n"]
        rows.append(dict(arm=arm, psms_sum_3files_3seeds=tot["n"], sample=tot["sample"], ups1=tot["ups1"], entrapment=tot["entrapment"],
                         unknown=tot["unknown"], fdp_lower=round(100 * lb, 3), fdp_combined=round(100 * lb * (1 + 1 / r), 3),
                         gained_vs_base=gained["n"], gained_entrapment=gained["entrapment"],
                         gained_fdp_combined=round(100 * gained["entrapment"] * (1 + 1 / r) / gained["n"], 2) if gained["n"] else ""))
    w = csv.DictWriter(sys.stdout, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main(sys.argv[1:] or ["base", "pr68", "pr66_priors", "p_zauto", "p_residue", "p_zauto_residue", "pr10335", "pr65_default"])
