"""Spread of the entrapment FDP over database realisations, per design and group, from summarize.py's per_group table.
Also the FDP when decoy accessions are ignored in the classification (a PSM is entrapment if all its non-decoy proteins
are ENTRAP_), counted from the Percolator target.tsv files.  realisations.py <per_group.tsv> <runs_dir>"""
import csv, statistics, sys
from pathlib import Path
DESIGNS = {"#10364 make_entrapment.py (C-term K/R/P fixed, N-term free)": ["frozenE", "frozenEs2", "frozenEs3"],
           "DecoyDatabase fix, keepPeptideNTerm true (default)": ["ddfixE1rev", "ddfixE1s2rev", "ddfixE1s3rev"],
           "DecoyDatabase fix, keepPeptideNTerm false": ["ddfixE1noNrev", "ddfixE1noNs2rev", "ddfixE1noNs3rev"],
           "DecoyDatabase b0ea4828 (Wen et al.: C-term + K/R/P fixed, unique)": ["wenE1s1rev", "wenE1s2rev", "wenE1s3rev"],
           "DecoyDatabase final (C-term + K/R/P fixed, identity threshold, no target collisions)": ["v3E1s1rev", "v3E1s2rev", "v3E1s3rev"],
           "andes build_shuffled_entrap.py (whole protein shuffled)": ["andesE1s1", "andesE1s2", "andesE1s3"]}
GROUPS = {"hfx": ["hfx_A2", "hfx_B1", "hfx_B3"], "astral": ["astral_A2", "astral_B1", "astral_B3"],
          "tmtpro": ["eclipse_tmtpro_10855", "eclipse_tmtpro_10858", "eclipse_tmtpro_10863"]}
rows = {(r["group"], r["arm"]): r for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")}
runs = Path(sys.argv[2])


def no_decoy_rule(group, arm, te):
    n = ne = 0
    for d in GROUPS[group]:
        for s in (1, 42, 137):
            for line in open(runs / d / arm / f"s{s}/target.tsv").read().splitlines()[1:]:
                f = line.split("\t")
                if float(f[2]) > 0.01: continue
                n += 1
                prots = [p for p in f[5:] if p and not p.startswith("DECOY_")]
                ne += bool(prots) and all(p.startswith("ENTRAP_") for p in prots)
    return 100 * ne / n * (1 + te)


print("group\tdesign\tFDP % per realisation\tmean\tsd\tFDP % ignoring decoy accessions (mean)\tPSMs per file (mean)")
for g in GROUPS:
    for name, arms in DESIGNS.items():
        if not all((g, a) in rows for a in arms): continue
        f = [float(rows[g, a]["fdp_pct"]) for a in arms]
        nd = [no_decoy_rule(g, a, float(rows[g, a]["te_ratio"])) for a in arms]
        psms = statistics.mean(float(rows[g, a]["psms_per_file"]) for a in arms)
        print(f"{g}\t{name}\t{' / '.join(f'{x:.2f}' for x in f)}\t{statistics.mean(f):.2f}\t{statistics.stdev(f):.2f}\t{statistics.mean(nd):.2f}\t{psms:.0f}")
