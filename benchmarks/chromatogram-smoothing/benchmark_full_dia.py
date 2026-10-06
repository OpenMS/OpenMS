#!/usr/bin/env python3
"""Extract and benchmark real DIA fragment chromatograms from an indexed mzML."""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import statistics
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pyopenms as oms

from benchmark_smoothers import chrom_from, memory, pick_summary, smoother_for

SCRIPT = Path(__file__).resolve()
BENCHMARK_DIR = SCRIPT.parent
RESULTS = BENCHMARK_DIR / "results"
PARAMETERS = RESULTS / "summary.json"
SOURCE_URL = "https://zenodo.org/records/18863866"
SOURCE_DOI = "10.5281/zenodo.18863866"
SOURCE_ACCESSION = "PXD005573"


def log(message):
    print(message, file=sys.stderr, flush=True)


def extract(path, n_targets):
    """Stream one frequent DIA MS2 isolation window into real fragment XICs."""
    experiment = oms.OnDiscMSExperiment()
    if not experiment.openFile(str(path)):
        raise RuntimeError(
            "OpenMS could not open this mzML as indexed mzML. "
            "Use the indexed file from the source archive."
        )

    metadata = experiment.getMetaData()
    by_window = defaultdict(list)
    for index in range(metadata.getNrSpectra()):
        spectrum = metadata.getSpectrum(index)
        if spectrum.getMSLevel() != 2:
            continue
        precursors = spectrum.getPrecursors()
        if not precursors:
            continue
        precursor = precursors[0]
        window = (
            round(precursor.getMZ(), 4),
            round(precursor.getIsolationWindowLowerOffset(), 4),
            round(precursor.getIsolationWindowUpperOffset(), 4),
        )
        by_window[window].append(index)

    if not by_window:
        raise RuntimeError("The mzML contains no MS2 isolation-window metadata.")
    window, indexes = max(by_window.items(), key=lambda item: len(item[1]))
    log(f"Selected MS2 window {window}; {len(indexes)} scans")

    # Rank recurring, high-intensity fragment masses in that window.
    bin_width = 0.02
    max_bin = int(2000 / bin_width) + 1
    scores = np.zeros(max_bin, dtype=np.float64)
    weighted_mz = np.zeros(max_bin, dtype=np.float64)
    for step, index in enumerate(indexes, 1):
        mz, intensity = experiment.getSpectrum(index).get_peaks()
        if not len(mz):
            continue
        valid = (mz >= 150.0) & (mz <= 2000.0) & (intensity > 0)
        valid_indexes = np.flatnonzero(valid)
        if valid_indexes.size == 0:
            continue
        count = min(128, valid_indexes.size)
        if valid_indexes.size > count:
            strongest = np.argpartition(intensity[valid_indexes], -count)[-count:]
            valid_indexes = valid_indexes[strongest]
        bins = np.floor(mz[valid_indexes] / bin_width).astype(np.int32)
        np.add.at(scores, bins, intensity[valid_indexes])
        np.add.at(weighted_mz, bins, mz[valid_indexes] * intensity[valid_indexes])
        if step % 2000 == 0:
            log(f"Scored {step}/{len(indexes)} scans")

    candidates = np.flatnonzero(scores > 0)
    candidates = candidates[np.argsort(scores[candidates])[::-1]]
    bins_selected = []
    for mass_bin in candidates:
        if all(abs(int(mass_bin) - previous) * bin_width >= 0.25 for previous in bins_selected):
            bins_selected.append(int(mass_bin))
            if len(bins_selected) == n_targets:
                break
    targets = np.asarray(
        [weighted_mz[b] / scores[b] for b in bins_selected], dtype=np.float64
    )
    if len(targets) < 10:
        raise RuntimeError(f"Only {len(targets)} recurring fragment masses were found.")

    # Extract each selected fragment XIC across every scan in the same DIA window.
    rt = np.empty(len(indexes), dtype=np.float64)
    traces = np.zeros((len(indexes), len(targets)), dtype=np.float64)
    for row, index in enumerate(indexes):
        spectrum = experiment.getSpectrum(index)
        rt[row] = spectrum.getRT()
        mz, intensity = spectrum.get_peaks()
        if len(mz) > 1 and np.any(mz[1:] < mz[:-1]):
            order = np.argsort(mz)
            mz, intensity = mz[order], intensity[order]
        for column, target in enumerate(targets):
            tolerance = max(0.005, target * 10e-6)
            left = np.searchsorted(mz, target - tolerance)
            right = np.searchsorted(mz, target + tolerance)
            if right > left:
                traces[row, column] = float(np.sum(intensity[left:right], dtype=np.float64))
        if (row + 1) % 2000 == 0:
            log(f"Extracted {row + 1}/{len(indexes)} scans")

    order = np.argsort(rt)
    rt, traces = rt[order], traces[order]
    intervals = np.diff(rt)
    if not np.isfinite(rt).all() or np.any(intervals <= 0):
        raise RuntimeError('Extracted retention times are not finite and strictly increasing.')
    if not np.isfinite(traces).all():
        raise RuntimeError('Extracted chromatograms contain non-finite intensities.')
    median_interval = float(np.median(intervals)) if len(intervals) else 0.0
    interval_cv = (
        float(np.std(intervals) / median_interval) if median_interval > 0 else None
    )
    info = {
        "ms2_window_center_mz": window[0],
        "isolation_lower_offset_mz": window[1],
        "isolation_upper_offset_mz": window[2],
        "spectra_in_window": len(indexes),
        "rt_start_seconds": float(rt[0]),
        "rt_end_seconds": float(rt[-1]),
        "median_rt_interval_seconds": median_interval,
        "rt_interval_cv": interval_cv,
        "fragment_mz_targets": targets.tolist(),
        "point_count_per_trace": int(len(rt)),
    }
    return rt, targets, traces, info


def percent_change(value, baseline):
    if value == "" or baseline == "" or baseline is None or baseline == 0:
        return ""
    return 100.0 * (float(value) - float(baseline)) / abs(float(baseline))


def memory_worker(method, archive, parameters, repeats):
    data = np.load(archive)
    rt, targets, traces = data["rt"], data["targets"], data["traces"]
    before = memory()
    rows, timings = [], []
    for index, target in enumerate(targets):
        intensity = traces[:, index]
        raw = chrom_from(rt, intensity)
        raw_count, raw_area, raw_width = pick_summary(raw)
        filt = smoother_for(method, parameters)
        per_trace_times = []
        smoothed = None
        for _ in range(repeats):
            candidate = copy.copy(raw)
            start = time.perf_counter()
            filt.filter(candidate)
            per_trace_times.append(1000.0 * (time.perf_counter() - start))
            smoothed = candidate
        timings.extend(per_trace_times)
        smoothed_count, smoothed_area, smoothed_width = pick_summary(smoothed)
        values = np.asarray(smoothed.get_peaks()[1], dtype=np.float64)
        if values.shape != intensity.shape:
            raise RuntimeError(f'{method} changed the chromatogram point count.')
        if not np.isfinite(values).all():
            raise RuntimeError(f'{method} produced non-finite intensities.')
        raw_auc = float(np.trapezoid(intensity, rt))
        smoothed_auc = float(np.trapezoid(values, rt))
        rows.append({
            "fragment_mz": float(target),
            "point_count": int(len(rt)),
            "method": method,
            "raw_picked_peak_count": int(raw_count),
            "picked_peak_count": int(smoothed_count),
            "picked_peak_count_delta": int(smoothed_count - raw_count),
            "raw_picked_integrated_intensity": raw_area,
            "picked_integrated_intensity": smoothed_area,
            "picked_integrated_intensity_change_percent": percent_change(smoothed_area, raw_area),
            "raw_median_picked_fwhm_seconds": raw_width,
            "median_picked_fwhm_seconds": smoothed_width,
            "median_picked_fwhm_change_percent": percent_change(smoothed_width, raw_width),
            "raw_trace_auc": raw_auc,
            "smoothed_trace_auc": smoothed_auc,
            "trace_auc_change_percent": percent_change(smoothed_auc, raw_auc),
            "median_filter_ms": float(statistics.median(per_trace_times)),
        })
    after = memory()
    return {
        "method": method,
        "rows": rows,
        "repeats_per_trace": repeats,
        "median_filter_ms_per_trace": float(np.median([r["median_filter_ms"] for r in rows])),
        "p95_filter_ms_per_trace": float(np.percentile([r["median_filter_ms"] for r in rows], 95)),
        "total_filter_ms": float(sum(timings)),
        "working_set_before_mb": before["working_set_mb"],
        "working_set_after_mb": after["working_set_mb"],
        "peak_working_set_mb": after["peak_working_set_mb"],
        "private_bytes_before_mb": before["private_mb"],
        "private_bytes_after_mb": after["private_mb"],
        "peak_private_bytes_mb": after["peak_private_mb"],
    }


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def file_hashes(path):
    md5, sha256 = hashlib.md5(), hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            md5.update(chunk)
            sha256.update(chunk)
    return {"md5": md5.hexdigest(), "sha256": sha256.hexdigest()}


def load_chromatograms(archive):
    """Load and validate the compact, pre-extracted DIA chromatogram input."""
    with np.load(archive, allow_pickle=False) as data:
        rt = np.asarray(data["rt"], dtype=np.float64)
        targets = np.asarray(data["targets"], dtype=np.float64)
        traces = np.asarray(data["traces"], dtype=np.float64)
    if rt.ndim != 1 or targets.ndim != 1 or traces.shape != (len(rt), len(targets)):
        raise ValueError("DIA archive must contain 1D rt/targets and a matching 2D traces array.")
    if len(rt) < 2 or len(targets) < 1:
        raise ValueError("DIA archive must contain at least two time points and one target.")
    if not np.isfinite(rt).all() or np.any(np.diff(rt) <= 0):
        raise ValueError("DIA archive retention times must be finite and strictly increasing.")
    if not np.isfinite(targets).all() or not np.isfinite(traces).all():
        raise ValueError("DIA archive targets and intensities must be finite.")
    return rt, targets, traces


def run_worker(archive, method, params, repeats):
    command = [
        sys.executable, str(SCRIPT), "--worker", str(archive),
        "--method", method, "--parameters", json.dumps(params),
        "--repeats", str(repeats),
    ]
    result = subprocess.run(command, cwd=BENCHMARK_DIR, text=True, capture_output=True, check=True)
    return json.loads(result.stdout.strip().splitlines()[-1])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path)
    parser.add_argument(
        "--chromatograms", type=Path,
        help="Use a pre-extracted DIA NPZ archive and its adjacent JSON provenance file.",
    )
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--targets", type=int, default=64)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--worker", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--method", choices=("modified_sinc", "savitzky_golay"), help=argparse.SUPPRESS)
    parser.add_argument("--parameters", help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args.worker:
        print(json.dumps(memory_worker(
            args.method, args.worker, json.loads(args.parameters), args.repeats
        )))
        return

    if (args.input is None) == (args.chromatograms is None):
        parser.error('provide exactly one of --input or --chromatograms')
    RESULTS.mkdir(parents=True, exist_ok=True)
    if args.input is not None:
        if args.output_dir is None:
            parser.error('--output-dir is required when extracting from --input')
        args.output_dir.mkdir(parents=True, exist_ok=True)
        rt, targets, traces, extraction = extract(args.input, args.targets)
        archive = args.output_dir / "extracted_fragment_chromatograms.npz"
        np.savez_compressed(archive, rt=rt, targets=targets, traces=traces)
        source_file = {
            "name": args.input.name,
            "size_bytes": args.input.stat().st_size,
            "hashes": file_hashes(args.input),
        }
        archive_metadata = {
            "dataset": SOURCE_ACCESSION,
            "dataset_description": "HeLa DIA acquisition, 4-hour run",
            "source_url": SOURCE_URL,
            "source_doi": SOURCE_DOI,
            "source_file": source_file,
            "published_md5": "19b814e1bcc9b67afbdac6624428eb31",
            "extraction": {
                **extraction,
                "method": "OnDiscMSExperiment; most frequent MS2 isolation window; "
                          "highest cumulative-intensity nonadjacent fragment m/z targets; "
                          "+/-10 ppm, minimum 0.005 Da",
                "archive_file": archive.name,
            },
            "archive": {
                "file_name": archive.name,
                "size_bytes": archive.stat().st_size,
                "sha256": file_hashes(archive)["sha256"],
            },
        }
        archive.with_suffix(".json").write_text(
            json.dumps(archive_metadata, indent=2) + "\n", encoding="utf-8"
        )
    else:
        archive = args.chromatograms
        rt, targets, traces = load_chromatograms(archive)
        metadata_path = archive.with_suffix(".json")
        if not metadata_path.is_file():
            parser.error(f'archive provenance file not found: {metadata_path}')
        archive_metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        archive_info = archive_metadata.get("archive", {})
        if archive_info.get("size_bytes") != archive.stat().st_size:
            parser.error(f'archive size does not match provenance: {archive}')
        if archive_info.get("sha256") != file_hashes(archive)["sha256"]:
            parser.error(f'archive checksum does not match provenance: {archive}')
        source_file = archive_metadata["source_file"]
        extraction = archive_metadata["extraction"]
    archive = archive.resolve()

    with PARAMETERS.open(encoding="utf-8") as handle:
        selected = json.load(handle)["selected_parameters"]
    method_results = [
        run_worker(archive, method, selected[method], args.repeats)
        for method in ("modified_sinc", "savitzky_golay")
    ]
    rows = [row for result in method_results for row in result.pop("rows")]
    csv_path = RESULTS / "full_dia_chromatograms.csv"
    write_csv(csv_path, rows)
    summary = {
        "dataset": archive_metadata["dataset"],
        "dataset_description": archive_metadata["dataset_description"],
        "source_url": archive_metadata["source_url"],
        "source_doi": archive_metadata["source_doi"],
        "input_mode": "raw_mzML" if args.input is not None else "pre_extracted_npz",
        "input_file": source_file["name"],
        "input_size_bytes": source_file["size_bytes"],
        "input_hashes": source_file["hashes"],
        "published_md5": archive_metadata["published_md5"],
        "chromatogram_archive": archive_metadata["archive"],
        "openms_pyopenms": oms.__version__,
        "pyopenms_source_commit": "5d5cbff4053b281763a1a79bf69e81c27967cfdf",
        "selected_parameters_from_synthetic_training": selected,
        "peak_picker_parameters": {
            "method": "corrected", "signal_to_noise": 1.0, "use_gauss": False,
            "sgolay_frame_length": 3, "sgolay_polynomial_order": 2,
        },
        "extraction": {
            **extraction,
            "target_count": int(len(targets)),
        },
        "method_performance": method_results,
        "csv_results": csv_path.relative_to(SCRIPT.parents[2]).as_posix(),
        "interpretation_limits": [
            "One run and one DIA isolation window; not a multi-instrument benchmark.",
            "Targets are selected by cumulative intensity, so weak real fragments are not sampled in equal proportion.",
            "Real-data peak counts and areas are compared with the unsmoothed trace; "
            "there is no ground truth for the true peak shape.",
            "The selected parameters were tuned on synthetic traces before this external-data evaluation.",
            "Process memory includes Python, pyOpenMS, NumPy, and extracted XICs; it is not allocation attribution.",
        ],
    }
    summary_path = RESULTS / "full_dia_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "input_file": source_file["name"],
        "chromatogram_archive": archive.name,
        "extracted_traces": int(len(targets)),
        "points_per_trace": int(len(rt)),
        "summary": str(summary_path),
        "csv": str(csv_path),
    }, indent=2))


if __name__ == "__main__":
    main()
