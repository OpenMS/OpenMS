"""Entrapment FDP per instrument group and database arm, with EntrapmentEvaluator.

  summarize.py <bin_dir> <runs_dir> <databases_dir> <out_prefix> [arm ...]

For each group/arm: EntrapmentEvaluator on the 3 files x 3 Percolator seeds (s*/target.idXML) with -database (T/E from
the searched database, Trypsin/P as searched, I = L) and on the ProSE in-process Percolator output (inproc.idXML).
Writes <out_prefix>_per_file.tsv (per file: mean PSMs over seeds, entrapment PSMs, FDP of the pooled seeds) and
<out_prefix>_per_group.tsv (mean over files of the per-file seed means; FDP from N_E and N pooled over files and seeds,
as the #10364 entrapment scripts do)."""
import os, statistics, subprocess, sys, tempfile
from pathlib import Path

bin_dir, runs, dbs, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), sys.argv[4]
ONLY = set(sys.argv[5:])  # optional: only these arms
GROUPS = {"velos": ["velos_125_R1", "velos_5000_R2", "velos_25000_R3"], "hfx": ["hfx_A2", "hfx_B1", "hfx_B3"],
          "astral": ["astral_A2", "astral_B1", "astral_B3"], "tmtpro": ["eclipse_tmtpro_10855", "eclipse_tmtpro_10858", "eclipse_tmtpro_10863"]}
ARMS = {  # arm: (database for hye group, for human group, for hy group, label)
    "natural": (None, None, "entrapment/hy_natural_td.fasta", "natural entrapment (#10364), r = 3.376"),
    "frozenE": ("entrapment/hye_entrap_td.fasta", "entrapment/human_entrap_td.fasta", None, "#10364 make_entrapment.py, 1:1"),
    "frozenEs2": ("entrapment/hye_entrap_s2.fasta", "entrapment/human_entrap_s2.fasta", None, "make_entrapment_salted.py salt s2"),
    "frozenEs3": ("entrapment/hye_entrap_s3.fasta", "entrapment/human_entrap_s3.fasta", None, "make_entrapment_salted.py salt s3"),
    "ddE1": ("dd/hye_dd_r1.fasta", "dd/human_dd_r1.fasta", None, "DecoyDatabase -entrapment_ratio 1, PR head"),
    "ddfixE1": ("dd/hye_ddfix_r1.fasta", "dd/human_ddfix_r1.fasta", None, "DecoyDatabase -entrapment_ratio 1, with fix"),
    "ddfixE1rev": ("dd/hye_ddfix_r1_rev.fasta", "dd/human_ddfix_r1_rev.fasta", None,
                   "DecoyDatabase -entrapment_ratio 1, with fix, decoys replaced by full reversals as in #10364"),
    "ddfixE1s2rev": ("dd/hye_ddfix_r1_seed2_rev.fasta", "dd/human_ddfix_r1_seed2_rev.fasta", None,
                     "as ddfixE1rev, DecoyDatabase -seed 2"),
    "ddfixE1s3rev": ("dd/hye_ddfix_r1_seed3_rev.fasta", "dd/human_ddfix_r1_seed3_rev.fasta", None,
                     "as ddfixE1rev, DecoyDatabase -seed 3"),
    "ddfixE1noNs2rev": ("dd/hye_ddfix_r1_noN_seed2_rev.fasta", "dd/human_ddfix_r1_noN_seed2_rev.fasta", None,
                        "as ddfixE1noNrev, -seed 2"),
    "ddfixE1noNs3rev": ("dd/hye_ddfix_r1_noN_seed3_rev.fasta", "dd/human_ddfix_r1_noN_seed3_rev.fasta", None,
                        "as ddfixE1noNrev, -seed 3"),
    "ddfixE1noNrev": ("dd/hye_ddfix_r1_noN_rev.fasta", "dd/human_ddfix_r1_noN_rev.fasta", None,
                      "as ddfixE1rev, -Decoy:keepPeptideNTerm false (N-terminus shuffled too, as in #10364)"),
    "wenE1s1rev": ("dd/hye_wen_s1_rev.fasta", "dd/human_wen_s1_rev.fasta", None, "DecoyDatabase b0ea4828 (Wen et al. shuffle), -seed 1, full-reversal decoys"),
    "andesE1s1": ("andes/hye_andes_s1.fasta", "andes/human_andes_s1.fasta", None, "andes build_shuffled_entrap.py seed 1 (whole-protein shuffle), full-reversal decoys"),
    "wenE1s2rev": ("dd/hye_wen_s2_rev.fasta", "dd/human_wen_s2_rev.fasta", None, "DecoyDatabase b0ea4828 (Wen et al. shuffle), -seed 2, full-reversal decoys"),
    "andesE1s2": ("andes/hye_andes_s2.fasta", "andes/human_andes_s2.fasta", None, "andes build_shuffled_entrap.py seed 2 (whole-protein shuffle), full-reversal decoys"),
    "wenE1s3rev": ("dd/hye_wen_s3_rev.fasta", "dd/human_wen_s3_rev.fasta", None, "DecoyDatabase b0ea4828 (Wen et al. shuffle), -seed 3, full-reversal decoys"),
    "andesE1s3": ("andes/hye_andes_s3.fasta", "andes/human_andes_s3.fasta", None, "andes build_shuffled_entrap.py seed 3 (whole-protein shuffle), full-reversal decoys"),
    "foreignE": ("foreign/hye_arath_td.fasta", "foreign/human_arath_td.fasta", None, "foreign entrapment: reviewed Arabidopsis (UP000006548), full-reversal decoys"),
    "v3E1s1rev": ("dd/hye_v3_s1_rev.fasta", "dd/human_v3_s1_rev.fasta", None, "DecoyDatabase final (C-term + K/R/P fixed, identity threshold, no target collisions), -seed 1, full-reversal decoys"),
    "v3E1s2rev": ("dd/hye_v3_s2_rev.fasta", "dd/human_v3_s2_rev.fasta", None, "DecoyDatabase final (C-term + K/R/P fixed, identity threshold, no target collisions), -seed 2, full-reversal decoys"),
    "v3E1s3rev": ("dd/hye_v3_s3_rev.fasta", "dd/human_v3_s3_rev.fasta", None, "DecoyDatabase final (C-term + K/R/P fixed, identity threshold, no target collisions), -seed 3, full-reversal decoys"),
    "ddfixE2": ("dd/hye_ddfix_r2.fasta", None, None, "DecoyDatabase -entrapment_ratio 2, with fix"),
}


def ee(files, db):
    with tempfile.NamedTemporaryFile(suffix=".tsv") as t:
        subprocess.run([bin_dir / "EntrapmentEvaluator", "-in"] + files + ["-database", db, "-enzyme", "Trypsin/P", "-out", t.name,
                        "-no_progress"], check=True, capture_output=True,
                   env=dict(os.environ, LD_LIBRARY_PATH=str(Path(bin_dir).parent / "lib"), OPENMS_DATA_PATH=str(Path(bin_dir).parent / "share/OpenMS")))
        rows = {}
        lines = open(t.name).read().splitlines()
        head = lines[0].split("\t")
        for line in lines[1:]:
            f = dict(zip(head, line.split("\t")))
            if "level" in f:  # EntrapmentEvaluator since e1f3d743: one row per level
                if f["class"] != "all": continue
                key = f["file"]; r = rows.setdefault(key, {})
                if f["level"] == "psm":
                    r.update(n=int(f["accepted"]), ne=int(f["entrapment"]), fdp=float(f["combined_fdp"]), te=float(f["te_ratio"]))
                else:
                    r.update(pep=int(f["accepted"]), pep_e=int(f["entrapment"]), pep_paired=f["paired_fdp"])
            elif f["class"] == "all":
                rows[f["file"]] = dict(n=int(f["accepted_psms"]), ne=int(f["entrapment_psms"]), fdp=float(f["psm_fdp"]),
                                       pep=int(f["accepted_peptides"]), pep_e=int(f["entrapment_peptides"]), te=float(f["te_ratio"]))
        return rows


per_file, per_group = [], []
for group, datasets in GROUPS.items():
    for arm, (hye, human, hy, label) in ARMS.items():
        db = {"velos": hy, "hfx": hye, "astral": hye, "tmtpro": human}[group]
        if ONLY and arm not in ONLY: continue
        if db is None or not all((runs / d / arm / "s137/target.idXML").exists() for d in datasets): continue
        db = dbs / db
        seeds = [runs / d / arm / f"s{s}/target.idXML" for d in datasets for s in (1, 42, 137)]
        r = ee([str(p) for p in seeds], db)
        means, ne_means = [], []
        for d in datasets:
            rs = [r[str(runs / d / arm / f"s{s}/target.idXML")] for s in (1, 42, 137)]
            n, ne = sum(x["n"] for x in rs), sum(x["ne"] for x in rs)
            te = rs[0]["te"]
            means.append(n / 3); ne_means.append(ne / 3)
            per_file.append(dict(group=group, dataset=d, arm=arm, psms_mean=round(n / 3, 1), entrapment_psms_mean=round(ne / 3, 1),
                                 fdp_pct=round(100 * ne / n * (1 + te), 2), te_ratio=te))
        tot = r["all files"]
        inproc = ee([str(runs / d / arm / "inproc.idXML") for d in datasets], db)["all files"]
        per_group.append(dict(group=group, arm=arm, database=label, psms_per_file=round(statistics.mean(means), 1),
                              entrapment_per_file=round(statistics.mean(ne_means), 1), te_ratio=tot["te"],
                              fdp_pct=round(100 * tot["fdp"], 2), fdp_te1_pct=round(200 * tot["ne"] / tot["n"], 2),
                              peptides_pooled=tot["pep"], peptide_fdp_pct=round(100 * tot["pep_e"] / tot["pep"] * (1 + tot["te"]), 2),
                              inproc_psms_per_file=round(inproc["n"] / 3, 1), inproc_fdp_pct=round(100 * inproc["fdp"], 2)))
        print(per_group[-1], flush=True)
for name, rows in (("per_file", per_file), ("per_group", per_group)):
    with open(f"{out}_{name}.tsv", "w") as f:
        f.write("\t".join(rows[0]) + "\n")
        for row in rows: f.write("\t".join(str(v) for v in row.values()) + "\n")
