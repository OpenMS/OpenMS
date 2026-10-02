"""Sage v0.14.7 on the immutable OpenMS/ANDES 20-file benchmark."""

import argparse, collections, csv, hashlib, json, math, os, re, statistics, subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BUNDLE = ROOT.parent / "prose-andes-reproduction"
DATA = BUNDLE / "broad_benchmark"
ENV = dict(
    os.environ,
    RAYON_NUM_THREADS="4",
    OMP_NUM_THREADS="4",
    SAGE_LOG="info",
    OPENMS_DATA_PATH=str(BUNDLE / "nightly/pyopenms/share/OpenMS"),
    PYTHONPATH=str(BUNDLE / "nightly") + os.pathsep + str(BUNDLE / "deps"),
)
sys.path[:0] = [str(BUNDLE / "nightly"), str(BUNDLE / "deps")]
import pyopenms as p

INPUTS = json.loads((BUNDLE / "inputs.json").read_text())
DBS = json.loads((BUNDLE / "databases.json").read_text())
MANIFEST = json.loads((DATA / "manifest.json").read_text())


def sha(path):
    with Path(path).open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def config(r, out):
    fixed = {
        "C": p.ModificationsDB()
        .getModification("Carbamidomethyl (C)")
        .getDiffMonoMass()
    }
    if r["label"] != "none":
        fixed.update(
            {
                "K": p.ModificationsDB()
                .getModification(r["label"] + " (K)")
                .getDiffMonoMass(),
                "^": p.ModificationsDB()
                .getModification(r["label"] + " (N-term)")
                .getDiffMonoMass(),
            }
        )
    prec = 50 if r["label"] == "TMTpro" else 20
    return {
        "database": {
            "bucket_size": 8192,
            "enzyme": {
                "missed_cleavages": 2,
                "min_len": 7,
                "max_len": 40,
                "cleave_at": "KR",
                "restrict": None,
                "c_terminal": True,
                "semi_enzymatic": False,
            },
            "fragment_min_mz": 150.0,
            "fragment_max_mz": 2000.0,
            "peptide_min_mass": 100.0,
            "peptide_max_mass": 9000.0,
            "ion_kinds": ["b", "y"],
            "min_ion_index": 2,
            "static_mods": fixed,
            "variable_mods": {
                "M": [
                    p.ModificationsDB()
                    .getModification("Oxidation (M)")
                    .getDiffMonoMass()
                ]
            },
            "max_variable_mods": 1,
            "decoy_tag": "DECOY_",
            "generate_decoys": False,
            "fasta": str(DATA / "data" / f"{r['database']}_td.fasta"),
        },
        "precursor_tol": {"ppm": [-prec, prec]},
        "fragment_tol": {
            r["fragment_unit"].lower(): [
                -r["fragment_tolerance"],
                r["fragment_tolerance"],
            ]
        },
        "precursor_charge": [2, 5],
        "isotope_errors": [-1, 2],
        "deisotope": True,
        "chimera": False,
        "wide_window": False,
        "predict_rt": False,
        "min_peaks": 15,
        "max_peaks": 150,
        "min_matched_peaks": 4,
        "max_fragment_charge": None,
        "report_psms": 10,
        "output_directory": str(out),
        "mzml_paths": [str(DATA / "data" / f"{r['id']}.mzML")],
    }


def command(cmd, out, name):
    started = time.monotonic()
    with (out / (name + ".log")).open("w") as f:
        proc = subprocess.Popen(list(map(str, cmd)), stdout=f, stderr=f, env=ENV)
        _, status, usage = os.wait4(proc.pid, 0)
        proc.returncode = os.waitstatus_to_exitcode(status)
    rec = {
        "command": list(map(str, cmd)),
        "returncode": proc.returncode,
        "seconds": time.monotonic() - started,
        "max_rss_kib": usage.ru_maxrss,
        "user_seconds": usage.ru_utime,
        "system_seconds": usage.ru_stime,
    }
    (out / (name + ".command.json")).write_text(json.dumps(rec, indent=2))
    if proc.returncode:
        raise RuntimeError(
            f"{name} failed: " + (out / (name + ".log")).read_text()[-3000:]
        )
    return rec


def parse_sequence(text, label):
    nterm = None
    if text.startswith("["):
        m = re.match(r"^\[([+-]?[0-9.]+)\]-?", text)
        assert m, text
        nterm = float(m[1])
        text = text[m.end() :]
    tokens = list(re.finditer(r"([A-Z])(?:\[([+-]?[0-9.]+)\])?", text))
    assert "".join(m[0] for m in tokens) == text, text
    seq = p.AASequence.fromString("".join(m[1] for m in tokens))
    names = {"C": "Carbamidomethyl (C)", "M": "Oxidation (M)"}
    if label != "none":
        names["K"] = label + " (K)"
    for i, m in enumerate(tokens):
        if m[2] is not None:
            name = names[m[1]]
            mass = p.ModificationsDB().getModification(name).getDiffMonoMass()
            assert abs(float(m[2]) - mass) < 0.001, (m[0], name, mass)
            seq.setModification(i, name)
    if nterm is not None:
        assert label != "none"
        name = label + " (N-term)"
        assert (
            abs(nterm - p.ModificationsDB().getModification(name).getDiffMonoMass())
            < 0.001
        )
        seq.setNTerminalModification(name)
    assert (nterm is not None) == (label != "none")
    return seq


def normalize(r, out):
    spectra = {
        s["scan"]: s
        for s in json.loads((DATA / "provenance" / f"{r['id']}.json").read_text())[
            "scans"
        ]
    }
    with (out / "results.sage.pin").open() as f:
        rd = csv.DictReader(f, delimiter="\t")
        header = rd.fieldnames
        rows = list(rd)
    # Outcomes of Sage's built-in LDA/FDR are not inputs to the second classifier.
    removed = [
        "FileName",
        "retentiontime",
        "ion_mobility",
        "aligned_rt",
        "predicted_rt",
        "sqrt(delta_rt_model)",
        "predicted_mobility",
        "sqrt(delta_mobility)",
        "posterior_error",
    ]
    assert set(removed) <= set(header)
    header = [x for x in header if x not in removed]
    rows.sort(
        key=lambda d: (
            int(d["ScanNr"]),
            int(d["rank"]),
            d["Peptide"],
            int(d["Label"]),
            float(d["isotope_error"]),
        )
    )
    counts = collections.Counter()
    labels = collections.Counter()
    unique = set()
    max_mass_error = 0
    max_obs_error = 0
    for d in rows:
        scan = int(d["ScanNr"])
        s = spectra[scan]
        z = s["charge"]
        seq = parse_sequence(d["Peptide"], r["label"])
        calc = seq.getMonoWeight()
        obs = (s["precursor_mz"] - 1.007276466771) * z
        assert d[f"z={z}"] == "1"
        max_mass_error = max(max_mass_error, abs(calc - float(d["CalcMass"])))
        max_obs_error = max(max_obs_error, abs(obs - float(d["ExpMass"])))
        assert (
            abs(calc - float(d["CalcMass"])) < 0.02
            and abs(obs - float(d["ExpMass"])) < 0.02
        )
        iso = round(float(d["isotope_error"]) / 1.0033548378)
        assert (
            iso in range(-1, 3)
            and abs(float(d["isotope_error"]) - iso * 1.0033548378) < 1e-5
        )
        assert int(d["missed_cleavages"]) == sum(
            a in "KR" for a in seq.toUnmodifiedString()[:-1]
        )
        key = (scan, seq.toString(), int(d["Label"]), z, iso)
        assert key not in unique, key
        unique.add(key)
        prots = d["Proteins"].split(";")
        assert all(x.startswith("DECOY_") for x in prots) == (d["Label"] == "-1"), (
            d["Label"],
            prots,
        )
        d.update(
            SpecId=f"{r['id']}_{scan}_{counts[scan]}",
            ExpMass=obs,
            CalcMass=calc,
            Peptide="-." + seq.toString() + ".-",
        )
        counts[scan] += 1
        labels[d["Label"]] += 1
        assert all(
            math.isfinite(float(d[k]))
            for k in header
            if k not in ["SpecId", "Peptide", "Proteins"]
        )
    with (out / "input.pin.part").open("w") as f:
        wr = csv.DictWriter(f, fieldnames=header, delimiter="\t", extrasaction="ignore")
        wr.writeheader()
        wr.writerows(rows)
    (out / "input.pin.part").replace(out / "input.pin")
    assert len(labels) == 2 and max(counts.values()) <= 10
    checks = {
        "dataset": r["id"],
        "candidate_rows": len(rows),
        "spectra": len(counts),
        "labels": dict(labels),
        "removed_features": removed,
        "features": header,
        "sha256": sha(out / "input.pin"),
        "native_pin_sha256": sha(out / "results.sage.pin"),
        "database_sha256": DBS[r["database"]]["sha256"],
        "mzml_sha256": INPUTS[r["id"]]["selected_sha256"],
        "max_sage_vs_openms_calc_mass_error_Da": max_mass_error,
        "max_sage_vs_shared_observed_mass_error_Da": max_obs_error,
        "redundant_hypotheses": 0,
    }
    (out / "input_checks.json").write_text(json.dumps(checks, indent=2))
    return checks


def native_counts(out):
    with (out / "results.sage.tsv").open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    winners = {}
    for d in rows:
        scan = d["scannr"]
        score = float(d["hyperscore"])
        label = int(d["label"])
        if scan not in winners or (score, -label) > (
            winners[scan][0],
            -winners[scan][1],
        ):
            winners[scan] = (score, label)
    bins = collections.defaultdict(lambda: [0, 0])
    for v, lab in winners.values():
        bins[v][lab == -1] += 1
    t = d = best = 0
    for v, (nt, nd) in sorted(bins.items(), reverse=True):
        t += nt
        d += nd
        if (d + 1) / max(t, 1) <= 0.01:
            best = t
    own = [x for x in rows if x["label"] == "1" and float(x["spectrum_q"]) <= 0.01]
    return {
        "psms": best,
        "feature": "hyperscore",
        "scored_spectra": len(winners),
        "sage_builtin_spectrum_q_target_rows": len(own),
        "sage_builtin_spectrum_q_unique_scans": len(set(x["scannr"] for x in own)),
    }


def rescore(out, seed):
    check = json.loads((out / "input_checks.json").read_text())
    assert sha(out / "input.pin") == check["sha256"]
    target = out / f"s{seed}"
    target.mkdir(exist_ok=True)
    if not all(
        (target / name).exists()
        for name in ["percolator.command.json", "target.tsv", "decoy.tsv"]
    ):
        command(
            [
                BUNDLE / "percolator/usr/bin/percolator",
                "-Y",
                "-U",
                "--seed",
                seed,
                "--results-psms",
                target / "target.tsv",
                "--decoy-results-psms",
                target / "decoy.tsv",
                out / "input.pin",
            ],
            target,
            "percolator",
        )
    assert (
        json.loads((target / "percolator.command.json").read_text())["returncode"] == 0
    )
    assert sha(out / "input.pin") == check["sha256"]
    assert (
        f"Found {check['candidate_rows']} PSMs"
        in (target / "percolator.log").read_text()
    )
    with (target / "target.tsv").open() as f:
        accepted = [
            d for d in csv.DictReader(f, delimiter="\t") if float(d["q-value"]) <= 0.01
        ]
    scans = [int(d["PSMId"].rsplit("_", 2)[-2]) for d in accepted]
    assert len(scans) == len(set(scans))
    info = {
        "seed": seed,
        "accepted": len(accepted),
        "distinct_sequences": len(
            set(re.sub(r"\([^)]*\)|\[[^]]*\]", "", x["peptide"]) for x in accepted)
        ),
        "target_sha256": sha(target / "target.tsv"),
        "decoy_sha256": sha(target / "decoy.tsv"),
    }
    (target / "result.json").write_text(json.dumps(info, indent=2))
    return info


def run(r, output_root):
    dataset = r["id"]
    out = output_root / dataset
    if out.exists() and any(out.iterdir()):
        raise RuntimeError("Refusing to overwrite " + str(out))
    out.mkdir(parents=True, exist_ok=True)
    assert sha(DATA / "data" / f"{dataset}.mzML") == INPUTS[dataset]["selected_sha256"]
    assert (
        sha(DATA / "data" / f"{r['database']}_td.fasta") == DBS[r["database"]]["sha256"]
    )
    cfg = config(r, out)
    (out / "config.json").write_text(json.dumps(cfg, indent=2))
    rec = command(
        [
            ROOT / "bin/sage",
            "--disable-telemetry-i-dont-want-to-improve-sage",
            "--write-pin",
            "--batch-size",
            "1",
            out / "config.json",
        ],
        out,
        "search",
    )
    resolved = json.loads((out / "results.json").read_text())
    assert (
        resolved["database"]["generate_decoys"] is False
        and resolved["predict_rt"] is False
    )
    checks = normalize(r, out)
    native = native_counts(out)
    seeds = [rescore(out, s) for s in [1, 42, 137]]
    result = {
        "dataset": dataset,
        "version": json.loads((ROOT / "version.json").read_text()),
        "config": cfg,
        "resolved": resolved,
        "native": native,
        "seed_results": seeds,
        "seed_psms": [s["accepted"] for s in seeds],
        "mean_psms": statistics.mean(s["accepted"] for s in seeds),
        "input_checks": checks,
        "search": rec,
    }
    (out / "summary.json").write_text(json.dumps(result, indent=2))
    print(
        json.dumps(
            {
                "dataset": dataset,
                "seed_psms": result["seed_psms"],
                "native": native,
                "seconds": rec["seconds"],
            }
        ),
        flush=True,
    )
    return result


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("datasets", nargs="*")
    ap.add_argument("--output-root", type=Path, default=ROOT / "results")
    ap.add_argument("--wait-for-inputs", action="store_true")
    ap.add_argument("--resume", action="store_true")
    a = ap.parse_args()
    version = json.loads((ROOT / "version.json").read_text())
    assert sha(ROOT / "bin/sage") == version["binary_sha256"]
    assert (
        sha(BUNDLE / "percolator/usr/bin/percolator")
        == "1f067b5d438a3a88be8a88f636844baea824e239fd2c5c053462ae56fd0e7c15"
    )
    for r in MANIFEST:
        if a.datasets and r["id"] not in a.datasets:
            continue
        out = a.output_root.resolve() / r["id"]
        if a.resume and (out / "summary.json").exists():
            old = json.loads((out / "summary.json").read_text())
            assert old["input_checks"]["sha256"] == sha(out / "input.pin")
            print(
                json.dumps({"dataset": r["id"], "status": "already completed"}),
                flush=True,
            )
            continue
        while a.wait_for_inputs and not (DATA / "data" / f"{r['id']}.mzML").exists():
            time.sleep(10)
        run(r, a.output_root.resolve())
