"""Entrapment FDP with the frozen #10364 harness logic, for comparison with EntrapmentEvaluator.

  frozen_eval.py velos    <hy_td.fasta> <target.tsv>...   natural entrapment (pr-prototypes-20261001/entrapment_velos.py)
  frozen_eval.py shuffled <target.tsv>...                 r = 1 shuffled entrapment (sage-comparison-20261002/evaluate_entrapment.py)

Each target.tsv is a Percolator --results-psms file. The functions below are copied from those two scripts
(only the FASTA path and the result paths are arguments here). Prints one row per input file: accepted PSMs at
q <= 0.01 (N), entrapment PSMs (N_E) and, for velos, the sample/ups1/unknown categories."""
import csv, re, sys

UPS1 = """P00915 P00918 P01031 P69905 P68871 P41159 P02768 P0CG48 P0CG47 P62987 P62979 P04040 P00167 P01133 P02144 P15559
P62937 Q06830 P63165 P00709 P06730 P12081 P61626 Q15843 P02753 P16083 P63279 P01008 P61769 P55957 O76070 P08263 P01344
P01127 P10599 P99999 P06396 P09211 P01112 P01579 P02787 P02788 P05413 P10636 P10145 P02741 O00762 P01375 P51965 P08758""".split()


# --- entrapment_velos.py (verbatim apart from the FASTA argument) ---
def read_fasta(FASTA):
    entries, hdr, seq = [], None, []
    for line in open(FASTA):
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


def velos(fasta, files):
    entries = read_fasta(fasta)
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
    print(f"# UPS1 accessions not in FASTA: {missing}; r = {r:.6f}; T/E = 1/r = {1 / r:.6f}")
    cache = {}

    def category(pep):
        if pep not in cache:
            cache[pep] = "sample" if pep in db["sample"] else "ups1" if pep in db["ups1"] else "entrapment" if pep in db["human"] else "unknown"
        return cache[pep]

    print("file\tN\tN_E\tsample\tups1\tunknown\tfdp_combined")
    for path in files:
        acc = {}
        with open(path) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                if float(row["q-value"]) <= 0.01:
                    acc[row["PSMId"].rsplit("_", 1)[0]] = stripped(row["peptide"])
        tot = {"n": 0, "sample": 0, "ups1": 0, "entrapment": 0, "unknown": 0}
        for scan, pep in acc.items():
            tot["n"] += 1; tot[category(pep)] += 1
        print(f"{path}\t{tot['n']}\t{tot['entrapment']}\t{tot['sample']}\t{tot['ups1']}\t{tot['unknown']}\t"
              f"{tot['entrapment'] / tot['n'] * (1 + 1 / r) if tot['n'] else 0:.6f}")


# --- evaluate_entrapment.py (verbatim accepted(), PREFIX = ENTRAP_) ---
PREFIX = "ENTRAP_"


def accepted(path):
    out = {}
    for line in open(path).read().splitlines()[1:]:
        f = line.split("\t")
        if float(f[2]) <= 0.01:
            scan = f[0].rsplit("_", 2)[-2]
            out[scan] = (f[4], all(p.startswith(PREFIX) for p in f[5:] if p))
    return out


def shuffled(files):
    print("file\tN\tN_E\tfdp_combined")
    for path in files:
        acc = accepted(path)
        n, ne = len(acc), sum(e for _, e in acc.values())
        print(f"{path}\t{n}\t{ne}\t{2 * ne / n if n else 0:.6f}")


if __name__ == "__main__":
    if sys.argv[1] == "velos": velos(sys.argv[2], sys.argv[3:])
    else: shuffled(sys.argv[2:])
