"""Class-FDR benchmark runner (PR #10636). One ProSE search per (setting, dataset), then per Percolator seed:
PercolatorAdapter (in-process, SVM score) -> annotate.py -> FalseDiscoveryRate global and -PSM_class arms -> extract.py.

usage: run.py <setting> [datasets...]
"""
import json, os, subprocess, sys, time, shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
BIN = Path(os.environ.get("OPENMS_BIN", "/home/user/OpenMS-build/bin"))
PROTO = Path(os.environ.get("PROTO", "/home/user/bench/proto/z/prose-andes-reproduction"))
DB = Path(os.environ.get("DBDIR", "/home/user/bench/db"))
WORK = Path(os.environ.get("WORK", "/home/user/bench/work"))
SEEDS = (1, 42, 137)
HARD_MODS = ["Oxidation (M)", "Deamidated (N)", "Deamidated (Q)", "Phospho (S)", "Phospho (T)", "Phospho (Y)"]
SETTINGS = {
    # frozen #10364 parameters, r = 1 shuffled entrapment
    "std": dict(db=None, over={}, arms=["mod", "bench_mc", "bench_z"]),
    # more variable modifications: modified PSMs (deamidation, phospho on unenriched samples) become a hard minority
    "hardmods": dict(db=None, over={"modifications:variable": HARD_MODS, "modifications:variable_max_per_peptide": 2,
                           "database:chunk_size": 10000},  # memory: ~300 M peptidoforms
                     arms=["mod"]),
    # andes-like: HYE + Swiss-Prot human isoforms (UniProt 2026_03), class canonical vs isoform-specific
    # semi-specific search: semi-tryptic PSMs are a hard minority
    "semi": dict(db=None, over={"peptide:enzyme_specificity": "semi", "database:chunk_size": 10000}, arms=["bench_term", "mod"]),
    "iso": dict(db="hyeiso", over={}, arms=["bench_iso", "mod"]),
}
# immunopeptidomics (MHCquant2 PXD058436, DN12 liver HLA class I, techRep1, 8,000 MS2 spectra): unspecific cleavage
IMMUNO = {"liver_classI_r1": dict(mzml="/home/user/bench/dl/immuno/liver_classI_r1.mzML", db="human", params={
    "enzyme": "unspecific cleavage", "peptide:enzyme_specificity": "none", "low_memory": "true", "decoy_prefix": "DECOY_",
    "decoys": "auto", "peptide:min_size": 8, "peptide:max_size": 12, "peptide:missed_cleavages": 0,
    "precursor:min_charge": 1, "precursor:max_charge": 4, "precursor:mass_tolerance_lower": 20.0,
    "precursor:mass_tolerance_upper": 20.0, "precursor:mass_tolerance_unit": "ppm", "fragment:mass_tolerance": 30.0,
    "fragment:mass_tolerance_unit": "ppm", "modifications:fixed": [], "modifications:variable": ["Oxidation (M)"],
    "modifications:variable_max_per_peptide": 1, "report:top_hits": 10, "scoring:max_candidates_per_spectrum": 50})}
SETTINGS["immuno"] = dict(db=None, over={}, arms=["bench_len", "bench_z", "mod"])
FIXED_RE = "Carbamidomethyl|TMT6plex|TMTpro"

def sh(cmd, log):
    with open(log, "a") as f:
        f.write("$ " + " ".join(map(str, cmd)) + "\n"); f.flush()
        subprocess.run(list(map(str, cmd)), check=True, stdout=f, stderr=subprocess.STDOUT)

def prose_cmd(ds, setting, out):
    if ds in IMMUNO:
        im = IMMUNO[ds]
        cmd = [BIN / "ProSE", "-in", im["mzml"], "-database", DB / f"{im['db']}_entrap_td.fasta", "-out_idxml", out,
               "-threads", 4, "-no_progress", "-Search:rescore", "false"]
        for k, v in im["params"].items():
            cmd.append(f"-Search:{k}"); cmd.extend(v if isinstance(v, list) else [v])
        return cmd
    job = next(j for j in json.load(open(PROTO / "jobs.json"))
               if j["stage"] == "broad_benchmark" and j["arm"] == "prose_auto" and j["dataset"] == ds)
    man = next(m for m in json.load(open(PROTO / "broad_benchmark/manifest.json")) if m["id"] == ds)
    params = dict(job["parameters"]); params.update(SETTINGS[setting]["over"])
    params["rescore"] = "false"  # rescoring done below with PercolatorAdapter and explicit seeds
    skip = {"precursor:isotopes", "snes_enabled"}  # not ProSE TOPP parameters (frozen harness only)
    db = DB / f"{SETTINGS[setting]['db'] or man['database']}_entrap_td.fasta"
    cmd = [BIN / "ProSE", "-in", PROTO / "broad_benchmark/data" / f"{ds}.mzML", "-database", db, "-out_idxml", out,
           "-threads", 4, "-no_progress"]
    for k, v in params.items():
        if k in skip: continue
        cmd.append(f"-Search:{k}")
        cmd.extend(v if isinstance(v, list) else [v])
    return cmd

def run(setting, ds):
    d = WORK / setting / ds; d.mkdir(parents=True, exist_ok=True); log = d / "log.txt"
    if (d / "done").exists(): return
    search = d / "search.idXML"
    t0 = time.time()
    if not search.exists():
        tmp = d / "search.tmp.idXML"
        sh(prose_cmd(ds, setting, tmp), log); tmp.rename(search)
    t1 = time.time()
    for s in SEEDS:
        sd = d / f"s{s}"; sd.mkdir(exist_ok=True)
        perc, ann = sd / "perc.idXML", sd / "perc_annot.idXML"
        sh([BIN / "PercolatorAdapter", "-in", search, "-out", perc, "-seed", s, "-score_type", "svm",
            "-post_processing_tdc", "-threads", 4], log)
        sh(["python3", "-I", HERE / "annotate.py", perc, ann], log); perc.unlink()
        for arm in ["global"] + SETTINGS[setting]["arms"]:
            out = sd / f"{arm}.idXML"
            cmd = [BIN / "FalseDiscoveryRate", "-in", ann, "-out", out, "-PSM", "true", "-protein", "false",
                   "-FDR:PSM", 1, "-algorithm:add_decoy_peptides", "true", "-threads", 1]
            if arm != "global": cmd += ["-PSM_class", "modification" if arm == "mod" else arm]
            sh(cmd, log)
            sh(["python3", "-I", HERE / "extract.py", out, sd / f"{arm}.tsv.gz", FIXED_RE], log)
            if not (setting == "std" and s == 1 and arm == "global"): out.unlink()
        if not (setting == "std" and s == 1): ann.unlink()
    search.unlink()
    (d / "done").write_text(json.dumps(dict(search_s=round(t1 - t0), total_s=round(time.time() - t0))))

if __name__ == "__main__":
    setting = sys.argv[1]
    dss = sys.argv[2:] or [m["id"] for m in json.load(open(PROTO / "broad_benchmark/manifest.json"))]
    for ds in dss:
        if ds not in IMMUNO and not (PROTO / "broad_benchmark/data" / f"{ds}.mzML").exists():
            print("missing input", ds, flush=True); continue
        print(time.strftime("%H:%M:%S"), setting, ds, flush=True)
        try: run(setting, ds)
        except subprocess.CalledProcessError as e: print("FAILED", setting, ds, e, flush=True)
