"""Benchmark the untangled ProSE PRs against develop on the frozen 20-file suite (OpenMS#10364).

Every arm replays the frozen `native_dedup:<dataset>:default` parameters of the historical
reference and changes only the settings listed in ARMS. Keys a build does not define are
skipped by the driver and logged in search.log. PIN export, Percolator 3.09.0 (-Y -U,
seeds 1/42/137) and native TDC are identical to the 2026-09-30 PR ablation.
"""

import argparse
import concurrent.futures
import csv
import gzip
import hashlib
import json
import os
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BUNDLE = ROOT.parent / "prose-andes-reproduction"
DATA = BUNDLE / "broad_benchmark/data"
sys.path.insert(0, str(BUNDLE / "broad_benchmark"))
from analyze import native  # noqa: E402

ENV = dict(
    os.environ,
    OMP_NUM_THREADS="4",
    OPENMS_DATA_PATH=str(BUNDLE / "nightly/pyopenms/share/OpenMS"),
    PYTHONPATH=str(BUNDLE / "nightly") + os.pathsep + str(BUNDLE / "deps"),
)
MANIFEST = json.loads((BUNDLE / "broad_benchmark/manifest.json").read_text())
INPUTS = json.loads((BUNDLE / "inputs.json").read_text())
DATABASES = json.loads((BUNDLE / "databases.json").read_text())
HISTORICAL = {j["id"]: j for j in json.loads((BUNDLE / "jobs.json").read_text())}

MASS = {"scoring:method": "mass_accuracy"}
CAL = {"calibration:enabled": "true"}
PRIORS = {"annotate:self_trained_ion_priors": "true"}
ENTRAP_FILES = ("astral_A2", "astral_B1", "astral_B3", "tims_plasma_30min", "tims_plasma_50min")
ENTRAP_DB = {"hye": "hye_entrap_td.fasta", "human": "human_entrap_td.fasta"}
FOREIGN_DB = {"hye": "hye_foreign_td.fasta", "human": "human_foreign_td.fasta"}
# arm: (build, parameter delta, files: False/"all", True/"ppm" (high-resolution), "da" (low-resolution) or a tuple of ids)
ARMS = {
    "base": ("base", {}, False),
    "pr68": ("pr10368", {}, False),
    "pr66_priors": ("pr10366", PRIORS, False),
    "pr65_default": ("pr10365", {"calibration:enabled": "auto"}, False),
    "pr65s_mass": ("pr10365_scorer", dict(MASS, **{"calibration:enabled": "false"}), True),
    "pr65s_mass_cal": ("pr10365_scorer", dict(MASS, **{"calibration:enabled": "true"}), True),
    "pr65_mass_cal": ("pr10365", dict(MASS, **{"calibration:enabled": "true"}), True),
    "pr10335": ("pr10335", {}, False),
    # Prototype of the #10379 fixes (377d3c7): the defaults reproduce 16ff034.
    "c_default": ("proto10379", {"calibration:enabled": "auto"}, True),
    "h_prec": ("proto10379", dict(CAL, **{"calibration:apply": "precursor"}), True),
    "h_frag": ("proto10379", dict(CAL, **{"calibration:apply": "fragment"}), True),
    "h_prec_robust": ("proto10379", dict(CAL, **{"calibration:apply": "precursor", "calibration:precursor_window": "robust"}), True),
    "h_prec_robust_stride": ("proto10379", dict(CAL, **{"calibration:apply": "precursor", "calibration:precursor_window": "robust",
                                                       "calibration:subset": "stride"}), True),
    "m_full_noapply": ("proto10379", dict(MASS, **CAL, **{"calibration:apply": "none", "scoring:mass_error_kernel_fit": "full"}), True),
    "m_shift_noapply": ("proto10379", dict(MASS, **CAL, **{"calibration:apply": "none", "scoring:mass_error_kernel_fit": "shift"}), True),
    "m_shift_prec_robust_stride": ("proto10379", dict(MASS, **CAL, **{"calibration:apply": "precursor", "calibration:precursor_window": "robust",
                                                                      "calibration:subset": "stride", "scoring:mass_error_kernel_fit": "shift"}), True),
    "m_shift_frag": ("proto10379", dict(MASS, **CAL, **{"calibration:apply": "fragment", "scoring:mass_error_kernel_fit": "shift"}), True),
    # Final #10379 commit (auto: fragment window only; kernel fit 'shift' by default).
    "f79_default": ("fix10379", {"calibration:enabled": "auto"}, tuple(r["id"] for r in MANIFEST if r["fragment_unit"] == "ppm")
                    + ("velos_125_R1", "lumos_tmt_5058")),
    "f79_mass": ("fix10379", dict(MASS, **{"calibration:enabled": "auto"}), True),
    # #10378 b1c49fd: cross-fitted ion priors (3 spectrum folds), fragment charges from the deisotoping decision.
    "f78x_priors": ("fix10378", PRIORS, False),
    "f78z_check": ("fix10378", PRIORS, ("velos_125_R1", "astral_A2")),  # 523214b vs ef0dcee (f78x_priors)
    # Peak-filter test (ANDES applies its top-20-per-100-Da filter only to TMT/iTRAQ: -14% on Astral otherwise).
    "wt40": ("base", {"peaks:window_top": 40}, "all"),
    "wt100": ("base", {"peaks:window_top": 100}, "all"),
    # Entrapment (entrapment/make_entrapment.py): targets + shuffled paired-mass entrapment (r = 1) + reversed decoys.
    # Dense-spectrum quota (peakfilter branch off develop): dense high-res spectra keep peaks:dense_window_top per window.
    "pf_off": ("peakfilter_bench", {"peaks:dense_window_top": 0}, ("hfx_A2", "astral_A2", "tims_plasma_30min", "velos_125_R1")),
    "pf": ("peakfilter_bench", {}, "all"),
    "pf_l10": ("peakfilter_bench", {"peaks:dense_intensity_loss": 0.1}, "ppm"),
    "pf_l30": ("peakfilter_bench", {"peaks:dense_intensity_loss": 0.3}, "ppm"),
    "pf_E": ("peakfilter_bench", {}, ENTRAP_FILES),
    # Sage comparison prototypes (fragz: develop f5ea2d04 + scoring:max_fragment_charge).
    "fz_off": ("fragz", {}, ("eclipse_tmtpro_10855", "velos_125_R1")),
    "fz2": ("fragz", {"scoring:max_fragment_charge": 2}, "all"),
    "fz3": ("fragz", {"scoring:max_fragment_charge": 3}, "all"),
    "mm4": ("base", {"fragment:min_matched_ions": 4}, "all"),
    # deiso: fragz + fragment:deisotope_start_check (Deisotoper start_intensity_check).
    "tbase_E": ("base", {}, ("eclipse_tmtpro_10855", "eclipse_tmtpro_10858", "eclipse_tmtpro_10863")),
    "tds1_E": ("deiso", {"fragment:deisotope_start_check": 1}, ("eclipse_tmtpro_10855", "eclipse_tmtpro_10858", "eclipse_tmtpro_10863")),
    "hds1_E": ("deiso", {"fragment:deisotope_start_check": 1}, ("hfx_A2", "hfx_B1", "hfx_B3", "astral_A2", "astral_B1", "astral_B3")),
    "hds1_mm4_E": ("deiso", {"fragment:deisotope_start_check": 1, "fragment:min_matched_ions": 4}, ("hfx_A2", "hfx_B1", "hfx_B3", "astral_A2", "astral_B1", "astral_B3")),
    "hbase_E": ("base", {}, ("hfx_A2", "hfx_B1", "hfx_B3")),
    "ds_off": ("deiso", {}, ("eclipse_tmtpro_10855", "hfx_A2")),
    "dfix": ("deisofix_bench", {}, ("eclipse_tmtpro_10855", "astral_A2", "velos_125_R1")),
    "ds1": ("deiso", {"fragment:deisotope_start_check": 1}, "ppm"),
    "ds1_fz2": ("deiso", {"fragment:deisotope_start_check": 1, "scoring:max_fragment_charge": 2}, "all"),
    "ds1_mm4": ("deiso", {"fragment:deisotope_start_check": 1, "fragment:min_matched_ions": 4}, "ppm"),
    "ds1_mm3": ("deiso", {"fragment:deisotope_start_check": 1, "fragment:min_matched_ions": 3}, "ppm"),
    "ds1_fz2_mm4": ("deiso", {"fragment:deisotope_start_check": 1, "scoring:max_fragment_charge": 2, "fragment:min_matched_ions": 4}, "all"),
    "fz2_mm4": ("fragz", {"scoring:max_fragment_charge": 2, "fragment:min_matched_ions": 4}, "all"),
    "base_E": ("base", {}, ENTRAP_FILES),
    # Foreign-like control (entrapment/make_foreign.py): same size as the _E databases, fully shuffled proteins.
    "base_F": ("base", {}, ENTRAP_FILES),
    "pf_F": ("peakfilter_bench", {}, ENTRAP_FILES),
    "f78x_E": ("fix10378", PRIORS, ENTRAP_FILES),
    "wt100_E": ("base", {"peaks:window_top": 100}, ENTRAP_FILES),
    # Final #10378 commit b02b538 (fragment charges from the deisotoping decision, no parameter).
    "f78_priors": ("fix10378", PRIORS, ("velos_125_R1", "velos_5000_R2", "velos_25000_R3", "lumos_tmt_5058", "lumos_tmt_5059",
                   "lumos_tmt_5066", "hfx_A2", "astral_A2")),
    # Prototype of the #10378 fixes (129e4d4): the defaults reproduce 241f9a7.
    "p_default": ("proto10378", PRIORS, ("velos_125_R1", "hfx_A2")),
    "p_zauto": ("proto10378", dict(PRIORS, **{"annotate:ion_prior_fragment_charges": "auto"}), ("velos_125_R1", "velos_5000_R2", "velos_25000_R3",
                "lumos_tmt_5058", "lumos_tmt_5059", "lumos_tmt_5066", "hfx_A2")),
    "p_residue": ("proto10378", dict(PRIORS, **{"annotate:ion_prior_residue_context": "true"}), False),
    "p_zauto_residue": ("proto10378", dict(PRIORS, **{"annotate:ion_prior_fragment_charges": "auto",
                                                      "annotate:ion_prior_residue_context": "true"}), "da"),
}


def selected(row, files):
    """Whether an arm's file selector includes a manifest row."""
    if files in (False, "all"):
        return True
    if files in (True, "ppm"):
        return row["fragment_unit"] == "ppm"
    if files == "da":
        return row["fragment_unit"] == "Da"
    return row["id"] in files


def sha(path):
    with Path(path).open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def command(cmd, out, name):
    started = time.monotonic()
    cmd = list(map(str, cmd))
    log_part = out / f"{name}.log.part"
    with log_part.open("w") as log:
        proc = subprocess.Popen(cmd, env=ENV, stdout=log, stderr=log, cwd=ROOT)
        _, status, usage = os.wait4(proc.pid, 0)
        proc.returncode = os.waitstatus_to_exitcode(status)
    log_part.replace(out / f"{name}.log")
    rec = dict(command=cmd, returncode=proc.returncode, seconds=time.monotonic() - started,
               max_rss_kib=usage.ru_maxrss, user_seconds=usage.ru_utime, system_seconds=usage.ru_stime)
    (out / f"{name}.command.json").write_text(json.dumps(rec, indent=2))
    if proc.returncode:
        raise RuntimeError(f"{out}: {name} failed ({proc.returncode}): " + (out / f"{name}.log").read_text()[-2500:])
    return rec


def rescore(out, seed, check):
    assert sha(out / "input.pin") == check["sha256"]
    folder = out / f"s{seed}"
    folder.mkdir()
    rec = command([BUNDLE / "percolator/usr/bin/percolator", "-Y", "-U", "--seed", seed,
                   "--results-psms", folder / "target.tsv", "--decoy-results-psms", folder / "decoy.tsv",
                   out / "input.pin"], folder, "percolator")
    assert f"Found {check['candidate_rows']} PSMs" in (folder / "percolator.log").read_text()
    with (folder / "target.tsv").open() as f:
        rows = [r for r in csv.DictReader(f, delimiter="\t") if float(r["q-value"]) <= 0.01]
    assert len({r["PSMId"].rsplit("_", 2)[-2] for r in rows}) == len(rows)
    result = dict(seed=seed, accepted=len(rows), target_sha256=sha(folder / "target.tsv"),
                  decoy_sha256=sha(folder / "decoy.tsv"), command=rec)
    (folder / "result.json").write_text(json.dumps(result, indent=2))
    return result


def search_and_score(binary, mzml, fasta, params, out, dataset, arm, export_extra=()):
    (out / "params.json").write_text(json.dumps(params, indent=2))
    search = command([binary, mzml, fasta, out / "params.json", out / "native.tsv", out / "search.idXML"], out, "search")
    command(["python3", ROOT / "export.py", dataset, arm, out, *export_extra], out, "export")
    check = json.loads((out / "input_checks.json").read_text())
    seeds = [rescore(out, seed, check) for seed in [1, 42, 137]]
    ignored = [line.split(": ", 1)[1] for line in (out / "search.log").read_text().splitlines()
               if line.startswith("[harness] ignored parameter")]
    return search, check, seeds, ignored


def finish(out, result):
    for filename in ["search.idXML", "native.pin", "native.tsv"]:
        source = out / filename
        with source.open("rb") as src, gzip.open(str(source) + ".gz.part", "wb", compresslevel=1) as dst:
            shutil.copyfileobj(src, dst)
        Path(str(source) + ".gz.part").replace(Path(str(source) + ".gz"))
        source.unlink()
    (out / "summary.json.part").write_text(json.dumps(result, indent=2))
    (out / "summary.json.part").replace(out / "summary.json")


def build_record(build):
    record = json.loads((ROOT / "build_records.json").read_text())[build]
    binary = ROOT / "builds" / build / "prose"
    assert sha(binary) == record["binary_sha256"] and record["tests_passed"], build
    return record, binary


def run(row, arm, output_root, resume):
    dataset = row["id"]
    out = output_root / dataset / arm
    if resume and (out / "summary.json").exists():
        saved = json.loads((out / "summary.json").read_text())
        # input.pin of finished jobs may have been removed to save space; native.pin.gz regenerates it
        if (out / "input.pin").exists():
            assert sha(out / "input.pin") == saved["input_checks"]["sha256"]
        return saved
    if out.exists() and any(out.iterdir()):
        raise RuntimeError("Refusing to overwrite incomplete output: " + str(out))
    out.mkdir(parents=True, exist_ok=True)
    build, delta, _ = ARMS[arm]
    record, binary = build_record(build)
    mzml = DATA / f"{dataset}.mzML"
    entrapment = arm.endswith(("_E", "_F"))
    fasta = (ROOT / "entrapment" / (ENTRAP_DB if arm.endswith("_E") else FOREIGN_DB)[row["database"]]) if entrapment \
        else DATA / f"{row['database']}_td.fasta"
    assert sha(mzml) == INPUTS[dataset]["selected_sha256"]
    assert entrapment or sha(fasta) == DATABASES[row["database"]]["sha256"]
    previous = HISTORICAL[f"native_dedup:{dataset}:default"]
    params = dict(previous["parameters"])
    params.update(delta)
    search, check, seeds, ignored = search_and_score(binary, mzml, fasta, params, out, dataset, arm)
    assert check["mzml_sha256"] == INPUTS[dataset]["selected_sha256"]
    if entrapment:
        check["database_sha256"] = sha(fasta)  # export.py records the benchmark database; this arm searched the entrapment one
        check["database"] = str(fasta.relative_to(ROOT))
    else:
        assert check["database_sha256"] == DATABASES[row["database"]]["sha256"]
    result = dict(
        dataset=dataset, arm=arm, build=record, parameters=params, ignored_parameters=ignored,
        resolved=json.loads((out / "resolved.json").read_text()), input_checks=check,
        seed_results=seeds, seed_psms=[s["accepted"] for s in seeds],
        mean_psms=statistics.mean(s["accepted"] for s in seeds), native=native(out), search=search,
        native_tsv_sha256=sha(out / "native.tsv"), native_pin_sha256=sha(out / "native.pin"),
    )
    finish(out, result)
    print(json.dumps(dict(dataset=dataset, arm=arm, seed_psms=result["seed_psms"], native=result["native"],
                          seconds=round(search["seconds"], 1), rss_mib=round(search["max_rss_kib"] / 1024))), flush=True)
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("datasets", nargs="*")
    ap.add_argument("--arms", nargs="+", choices=list(ARMS), default=list(ARMS))
    ap.add_argument("--output-root", type=Path, default=ROOT / "results")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--workers", type=int, default=2)
    args = ap.parse_args()
    assert sha(BUNDLE / "percolator/usr/bin/percolator") == "1f067b5d438a3a88be8a88f636844baea824e239fd2c5c053462ae56fd0e7c15"
    jobs = [(r, arm) for r in MANIFEST for arm in args.arms
            if (not args.datasets or r["id"] in args.datasets) and selected(r, ARMS[arm][2])]
    print("Planned searches:", len(jobs), flush=True)
    failures = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run, r, arm, args.output_root.resolve(), args.resume): (r["id"], arm) for r, arm in jobs}
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as error:  # report every failure, never hide one
                failures.append((futures[future], repr(error)[-600:]))
                print("FAILED", futures[future], repr(error)[-600:], flush=True)
    print("Done:", len(jobs) - len(failures), "ok,", len(failures), "failed", flush=True)
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
