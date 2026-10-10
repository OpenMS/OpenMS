#!/usr/bin/env python3
"""Reproducible OpenMS chromatogram smoother benchmark for issue #10425."""

from __future__ import annotations
import argparse
import copy
import csv
import hashlib
import itertools
import json
import math
import os
import platform
import statistics
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pyopenms as oms

SCRIPT = Path(__file__).resolve()
REPO = SCRIPT.parents[2]
FIXTURE_DIR = REPO / "src/tests/topp"
OUT = SCRIPT.parent / "results"
FIXTURES = (
    "NoiseFilterSGolay_2_input.chrom.mzML",
    "OpenSwathWorkflow_17_output.chrom.mzML",
    "OpenSwathWorkflow_22_output.chrom.mzML",
    "OpenSwathWorkflow_23_output.chrom.mzML",
)
TRAIN = (101, 202, 303, 404)
VALIDATE = tuple(range(1001, 1013))
SHAPES = ("narrow", "broad", "overlapping", "low_intensity")
SYNTHETIC_RT_STEP_SECONDS = 5.26686


def synthetic(shape, seed):
    """Create seeded Gaussian chromatograms sampled at the measured DIA RT cadence."""
    rng = np.random.default_rng(seed)
    rt_step = SYNTHETIC_RT_STEP_SECONDS
    point_index = np.arange(401, dtype=np.float64)
    rt = point_index * rt_step
    baseline = 100.0
    # Peak centers and widths are defined in points so filter behavior is unchanged.
    configs = {
        "narrow": ([(200.0, 1000.0, 2.6)], 24.0),
        "broad": ([(200.0, 1000.0, 16.0)], 24.0),
        "overlapping": ([(188.0, 950.0, 6.5), (212.0, 720.0, 8.5)], 22.0),
        "low_intensity": ([(200.0, 75.0, 6.0)], 18.0),
    }
    components, noise_sd = configs[shape]
    clean = np.full(rt.shape, baseline, dtype=np.float64)
    for center, height, sigma in components:
        clean += height * np.exp(-0.5 * ((point_index - center) / sigma) ** 2)
    components_seconds = [
        (center * rt_step, height, sigma * rt_step)
        for center, height, sigma in components
    ]
    return rt, clean + rng.normal(0.0, noise_sd, rt.size), clean, components_seconds


def chrom_from(rt, intensity):
    chrom = oms.MSChromatogram()
    chrom.set_peaks(np.asarray(rt, dtype=np.float64), np.asarray(intensity, dtype=np.float64))
    return chrom


def smoother_for(method, values):
    smoother = oms.ModifiedSincSmoother() if method == "modified_sinc" else oms.SavitzkyGolayFilter()
    params = smoother.getParameters()
    for key, value in values.items():
        params.setValue(key, value)
    smoother.setParameters(params)
    return smoother


def filter_values(method, params, rt, intensity):
    chrom = chrom_from(rt, intensity)
    smoother_for(method, params).filter(chrom)
    output = np.asarray(chrom.get_peaks()[1], dtype=np.float64)
    if output.shape != np.asarray(intensity).shape:
        raise RuntimeError(f"{method} changed the chromatogram point count")
    if not np.isfinite(output).all():
        raise RuntimeError(f"{method} produced non-finite intensity values")
    return output


def candidates():
    sinc = [
        {"is_ms1": is_ms1, "degree": degree, "m": m}
        for is_ms1, degree, m in itertools.product((False, True), (4, 6, 8, 10), (2, 3, 4, 5, 7, 9, 11, 13, 15, 19, 23, 27))
        if m >= degree // 2 + (1 if is_ms1 else 2)
    ]
    sg = [
        {"frame_length": frame, "polynomial_order": order}
        for frame, order in itertools.product((5, 7, 9, 11, 15, 21, 31), (2, 3, 4))
        if order < frame
    ]
    return {"modified_sinc": sinc, "savitzky_golay": sg}


def tune():
    training = [synthetic(shape, seed) for shape in SHAPES for seed in TRAIN]
    selected, rows = {}, []
    for method, parameter_sets in candidates().items():
        scores = []
        for params in parameter_sets:
            errors = []
            for rt, noisy, clean, peaks in training:
                scale = max(amp for _, amp, _ in peaks)
                predicted = filter_values(method, params, rt, noisy)
                errors.append(float(np.sqrt(np.mean((predicted - clean) ** 2)) / scale))
            score = float(np.mean(errors))
            scores.append((score, params))
            rows.append({"method": method, "parameters": json.dumps(params, sort_keys=True),
                         "training_mean_normalized_rmse": score, "selected": False})
        score, params = min(scores, key=lambda item: item[0])
        selected[method] = params
        for row in rows:
            if row["method"] == method and row["parameters"] == json.dumps(params, sort_keys=True):
                row["selected"] = True
                row["training_mean_normalized_rmse"] = score
                break
    return selected, rows


def crossing(x0, y0, x1, y1, target):
    return (x0 + x1) / 2.0 if y1 == y0 else x0 + (target - y0) * (x1 - x0) / (y1 - y0)



def peak_regions(rt, peaks):
    """Partition a synthetic trace so overlapping peaks have separate local ROIs."""
    regions = []
    centers = [center for center, _, _ in peaks]
    for index, (center, _, sigma) in enumerate(peaks):
        left = centers[index - 1] + (center - centers[index - 1]) / 2.0 if index else center - 5.0 * sigma
        right = center + (centers[index + 1] - center) / 2.0 if index + 1 < len(peaks) else center + 5.0 * sigma
        lo = max(0, int(np.searchsorted(rt, left, side="left")))
        hi = min(rt.size, int(np.searchsorted(rt, right, side="right")))
        regions.append((lo, hi))
    return regions


def measured_region_fwhm(rt, values, lo, hi):
    """Measure half-prominence width inside one peak's assigned region."""
    if hi - lo < 3:
        return None
    apex = lo + int(np.argmax(values[lo:hi]))
    baseline = max(float(values[lo]), float(values[hi - 1]))
    if values[apex] <= baseline:
        return None
    half = baseline + (float(values[apex]) - baseline) / 2.0
    left = apex
    while left > lo and values[left] > half:
        left -= 1
    right = apex
    while right + 1 < hi and values[right] > half:
        right += 1
    if left == lo or right == hi - 1 or values[left] > half or values[right] > half:
        return None
    xleft = crossing(rt[left], values[left], rt[left + 1], values[left + 1], half)
    xright = crossing(rt[right - 1], values[right - 1], rt[right], values[right], half)
    return float(xright - xleft)


def shape_metrics(rt, clean, original, values, peaks):
    scale = max(amp for _, amp, _ in peaks)
    rmse = float(np.sqrt(np.mean((values - clean) ** 2)) / scale)
    background = np.ones(rt.size, dtype=bool)
    for center, _, sigma in peaks:
        background &= np.abs(rt - center) > 4.0 * sigma
    if np.count_nonzero(background) < 2:
        raise RuntimeError("Synthetic trace has too few background samples for noise measurement")
    input_noise = float(np.sqrt(np.mean((original[background] - clean[background]) ** 2)))
    output_noise = float(np.sqrt(np.mean((values[background] - clean[background]) ** 2)))
    noise_reduction = 100.0 * (1.0 - output_noise / max(input_noise, 1e-12))
    rt_errors, height_errors, width_errors = [], [], []
    for lo, hi in peak_regions(rt, peaks):
        clean_apex = lo + int(np.argmax(clean[lo:hi]))
        apex = lo + int(np.argmax(values[lo:hi]))
        rt_errors.append(abs(float(rt[apex]) - float(rt[clean_apex])))
        clean_height = max(float(clean[clean_apex]) - 100.0, 1e-12)
        height_errors.append(abs(float(values[apex]) - float(clean[clean_apex])) / clean_height * 100.0)
        expected_width = measured_region_fwhm(rt, clean, lo, hi)
        observed_width = measured_region_fwhm(rt, values, lo, hi)
        if expected_width is not None and observed_width is not None:
            width_errors.append(abs(observed_width - expected_width) / expected_width * 100.0)
    roi = (rt >= min(c - 4*s for c, _, s in peaks)) & (rt <= max(c + 4*s for c, _, s in peaks))
    area = float(np.trapezoid(values[roi] - 100.0, rt[roi]))
    expected_area = sum(h*s*math.sqrt(2.0*math.pi) for _, h, s in peaks)
    return {
        "normalized_rmse": rmse,
        "background_noise_rms": output_noise,
        "noise_reduction_percent": noise_reduction,
        "apex_rt_error_seconds_mean": float(np.mean(rt_errors)),
        "apex_height_error_percent_mean": float(np.mean(height_errors)),
        "area_error_percent": abs(area - expected_area) / expected_area * 100.0,
        "fwhm_error_percent": float(np.mean(width_errors)) if width_errors else "",
    }


def picker():
    obj = oms.PeakPickerChromatogram()
    params = obj.getParameters()
    params.setValue("method", "corrected")
    params.setValue("signal_to_noise", 1.0)
    params.setValue("use_gauss", False)
    params.setValue("sgolay_frame_length", 3)
    params.setValue("sgolay_polynomial_order", 2)
    obj.setParameters(params)
    return obj


def pick_summary(input_chrom):
    picked = oms.MSChromatogram()
    # The three-argument overload overwrites its third chromatogram argument.
    # Use the two-argument overload so the pre-smoothed input is actually picked.
    picker().pickChromatogram(input_chrom, picked)
    arrays = {a.getName(): np.asarray(list(a), dtype=np.float64) for a in picked.getFloatDataArrays()}
    areas, widths = arrays.get("IntegratedIntensity"), arrays.get("FWHM")
    area_sum = float(np.nansum(areas)) if areas is not None and areas.size else ""
    width_median = float(np.nanmedian(widths)) if widths is not None and widths.size else ""
    return picked.size(), area_sum, width_median


def load_chroms(path):
    exp = oms.MSExperiment()
    oms.MzMLFile().load(str(path), exp)
    return list(exp.getChromatograms())


def csv_write(path, rows):
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def synthetic_benchmark(params):
    rows, grouped = [], {}
    for shape in SHAPES:
        grouped[shape] = {}
        for seed in VALIDATE:
            rt, noisy, clean, peaks = synthetic(shape, seed)
            raw = chrom_from(rt, noisy)
            count, picked_area, picked_width = pick_summary(raw)
            metrics = shape_metrics(rt, clean, noisy, noisy, peaks)
            rows.append({"shape": shape, "seed": seed, "method": "raw", **metrics,
                         "picked_peak_count": count, "picked_integrated_intensity": picked_area,
                         "picked_median_fwhm_seconds": picked_width})
            grouped[shape].setdefault("raw", []).append(metrics["normalized_rmse"])
            for method in params:
                values = filter_values(method, params[method], rt, noisy)
                smooth = chrom_from(rt, values)
                count, picked_area, picked_width = pick_summary(smooth)
                metrics = shape_metrics(rt, clean, noisy, values, peaks)
                rows.append({"shape": shape, "seed": seed, "method": method, **metrics,
                             "picked_peak_count": count, "picked_integrated_intensity": picked_area,
                             "picked_median_fwhm_seconds": picked_width})
                grouped[shape].setdefault(method, []).append(metrics["normalized_rmse"])
    summary = {
        shape: {method: {"median_normalized_rmse": float(np.median(v)),
                         "p90_normalized_rmse": float(np.percentile(v, 90))}
                for method, v in methods.items()}
        for shape, methods in grouped.items()
    }
    return rows, summary


def real_benchmark(params):
    rows = []
    for filename in FIXTURES:
        for index, raw in enumerate(load_chroms(FIXTURE_DIR / filename)):
            rt, intensity = raw.get_peaks()
            rt, intensity = np.asarray(rt, dtype=float), np.asarray(intensity, dtype=float)
            if rt.size < 3:
                continue
            intervals = np.diff(rt)
            median_dt = float(np.median(intervals))
            cv = float(np.std(intervals) / median_dt) if median_dt else float("nan")
            for method in ("raw", *params.keys()):
                smoothed = copy.copy(raw)
                if method != "raw":
                    smoother_for(method, params[method]).filter(smoothed)
                picked, area, width = pick_summary(smoothed)
                values = np.asarray(smoothed.get_peaks()[1], dtype=float)
                rows.append({
                    "fixture": filename, "chromatogram_index": index, "point_count": rt.size,
                    "median_rt_interval_seconds": median_dt, "rt_interval_cv": cv,
                    "method": method, "picked_peak_count": picked,
                    "picked_integrated_intensity_sum": area, "median_picked_fwhm_seconds": width,
                    "trace_auc": float(np.trapezoid(values, rt)),
                })
    return rows


def memory():
    if os.name == "nt":
        import ctypes
        from ctypes import wintypes
        class Counters(ctypes.Structure):
            _fields_ = [
                ("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
                ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
                ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t), ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t),
            ]
        c = Counters()
        c.cb = ctypes.sizeof(c)
        kernel = ctypes.WinDLL("kernel32", use_last_error=True)
        psapi = ctypes.WinDLL("psapi", use_last_error=True)
        kernel.GetCurrentProcess.restype = wintypes.HANDLE
        psapi.GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.POINTER(Counters), wintypes.DWORD]
        psapi.GetProcessMemoryInfo.restype = wintypes.BOOL
        if not psapi.GetProcessMemoryInfo(kernel.GetCurrentProcess(), ctypes.byref(c), c.cb):
            raise ctypes.WinError(ctypes.get_last_error())
        mib = 1024.0**2
        return {"working_set_mb": c.WorkingSetSize/mib, "peak_working_set_mb": c.PeakWorkingSetSize/mib,
                "private_mb": c.PagefileUsage/mib, "peak_private_mb": c.PeakPagefileUsage/mib}
    import resource
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    mib = 1024.0**2 if platform.system() == "Darwin" else 1024.0
    return {"working_set_mb": None, "peak_working_set_mb": peak/mib,
            "private_mb": None, "peak_private_mb": None}


def worker(data):
    method, source, params = data["method"], data["source"], data["parameters"]
    if source == "large_synthetic_100k":
        rt = np.arange(100_000, dtype=float)
        intensity = 100 + 1000*np.exp(-0.5*((rt-50_000)/40)**2)
        intensity += np.random.default_rng(917).normal(0, 20, rt.size)
        template = chrom_from(rt, intensity)
    else:
        template = max((c for c in load_chroms(FIXTURE_DIR / source) if c.size() >= 3), key=lambda c:c.size())
    filt = smoother_for(method, params)
    for _ in range(2):
        filt.filter(copy.copy(template))
    before = memory()
    n = 6 if template.size() >= 50_000 else 40
    timings = []
    for _ in range(n):
        candidate = copy.copy(template)
        start = time.perf_counter()
        filt.filter(candidate)
        timings.append(1000*(time.perf_counter()-start))
    after = memory()
    return {"source": source, "method": method, "parameters": params, "point_count": template.size(),
            "repeats": n, "median_ms": statistics.median(timings),
            "p95_ms": float(np.percentile(timings, 95)), "before": before, "after": after}


def performance(params):
    rows = []
    for source, method in itertools.product((*FIXTURES, "large_synthetic_100k"), params):
        request = {"source": source, "method": method, "parameters": params[method]}
        process = subprocess.run([sys.executable, str(SCRIPT), "--worker", json.dumps(request)],
                                 cwd=REPO, text=True, capture_output=True, check=True)
        result = json.loads(process.stdout.strip().splitlines()[-1])
        rows.append({
            "source": source, "method": method, "parameters": json.dumps(params[method], sort_keys=True),
            "point_count": result["point_count"], "repeat_count": result["repeats"],
            "median_filter_ms": result["median_ms"], "p95_filter_ms": result["p95_ms"],
            "working_set_before_mb": result["before"]["working_set_mb"],
            "working_set_after_mb": result["after"]["working_set_mb"],
            "peak_working_set_mb": result["after"]["peak_working_set_mb"],
            "private_bytes_before_mb": result["before"]["private_mb"],
            "private_bytes_after_mb": result["after"]["private_mb"],
            "peak_private_bytes_mb": result["after"]["peak_private_mb"],
        })
    return rows


def implementation_sources_match_wheel():
    paths = [
        "src/openms/source/PROCESSING/SMOOTHING/ModifiedSincSmoother.cpp",
        "src/openms/include/OpenMS/PROCESSING/SMOOTHING/ModifiedSincSmoother.h",
        "src/openms/source/PROCESSING/SMOOTHING/SavitzkyGolayFilter.cpp",
        "src/openms/include/OpenMS/PROCESSING/SMOOTHING/SavitzkyGolayFilter.h",
        "src/openms/source/ANALYSIS/OPENSWATH/PeakPickerChromatogram.cpp",
        "src/openms/include/OpenMS/ANALYSIS/OPENSWATH/PeakPickerChromatogram.h",
    ]
    result = subprocess.run(
        ["git", "diff", "--quiet", "5d5cbff4053b281763a1a79bf69e81c27967cfdf", "HEAD", "--", *paths],
        cwd=REPO,
    )
    if result.returncode not in (0, 1):
        raise RuntimeError("Could not compare wheel and checkout implementation sources")
    return result.returncode == 0


def main():
    cli = argparse.ArgumentParser()
    cli.add_argument("--worker", help=argparse.SUPPRESS)
    args = cli.parse_args()
    if args.worker:
        print(json.dumps(worker(json.loads(args.worker))))
        return

    OUT.mkdir(parents=True, exist_ok=True)
    params, search = tune()
    synth_rows, synth_summary = synthetic_benchmark(params)
    real_rows = real_benchmark(params)
    perf_rows = performance(params)
    csv_write(OUT / "parameter_search.csv", search)
    csv_write(OUT / "synthetic_validation.csv", synth_rows)
    csv_write(OUT / "real_chromatograms.csv", real_rows)
    csv_write(OUT / "performance_memory.csv", perf_rows)
    metadata = {
        "issue": "OpenMS/OpenMS#10425",
        "source_commit": subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO, text=True,
                                        capture_output=True, check=True).stdout.strip(),
        "python": sys.version, "platform": platform.platform(),
        "pyopenms": oms.__version__,
        "pyopenms_source_commit": "5d5cbff4053b281763a1a79bf69e81c27967cfdf",
        "smoother_and_picker_sources_match_wheel_build": implementation_sources_match_wheel(),
        "numpy": np.__version__,
        "selected_parameters": params, "training_seeds_per_shape": len(TRAIN),
        "validation_seeds_per_shape": len(VALIDATE), "synthetic_rt_step_seconds": SYNTHETIC_RT_STEP_SECONDS,
        "synthetic_validation_summary": synth_summary,
        "fixture_sha256": {n: hashlib.sha256((FIXTURE_DIR/n).read_bytes()).hexdigest() for n in FIXTURES},
        "fixtures": list(FIXTURES),
        "peak_picker_parameters": {"method": "corrected", "signal_to_noise": 1.0, "use_gauss": False,
                                   "sgolay_frame_length": 3, "sgolay_polynomial_order": 2},
        "limitations": [
            "Tracked fixtures are small regression chromatograms, not full production DIA runs.",
            "Peak picking uses identical fixed settings on unsmoothed and smoothed inputs.",
            "Memory is process-level Python/OpenMS working-set/private-byte measurement, not allocation attribution.",
            "The 100k-point synthetic input is a stress check, not a replacement for full-file profiling.",
        ],
    }
    (OUT/"summary.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output_dir": str(OUT), "selected_parameters": params,
                      "synthetic_summary": synth_summary, "real_rows": len(real_rows),
                      "performance_rows": len(perf_rows)}, indent=2))


if __name__ == "__main__":
    main()
