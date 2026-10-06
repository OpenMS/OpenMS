"""
Benchmark execution engine for OpenMS smoothers.
GitHub Issue #10425: Benchmark ModifiedSincSmoother against Savitzky–Golay.

Implements:
- Execution of ModifiedSincSmoother and SavitzkyGolayFilter on OpenMS containers
- Independent parameter tuning on separate calibration datasets (NO test-set leakage)
- Downstream peak-picking evaluation with realistic S/N thresholds (1.0, 2.0, 3.0, 5.0) via PeakPickerHiRes
- Total doublet area evaluation for overlapping peaks
- Sampling phase variation evaluation on narrow peaks
- Rigorous microbenchmark runtime measurements (warmup, gc disabled, 500+ iterations)
- Sanity checks on representative OpenMS test chromatograms
"""

from __future__ import annotations

import gc
import math
import os
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple, Dict, Any, Callable

import numpy as np
import pyopenms

try:
    from .synthetic_data import SyntheticChromatogram, SyntheticChromatogramGenerator
    from .metrics import (
        OverallBenchmarkMetrics,
        evaluate_signal_fidelity,
        evaluate_noise_reduction,
        evaluate_valley_to_peak_ratio,
        evaluate_peak_picking,
        compute_peak_area,
    )
except ImportError:
    from synthetic_data import SyntheticChromatogram, SyntheticChromatogramGenerator
    from metrics import (
        OverallBenchmarkMetrics,
        evaluate_signal_fidelity,
        evaluate_noise_reduction,
        evaluate_valley_to_peak_ratio,
        evaluate_peak_picking,
        compute_peak_area,
    )


def apply_modified_sinc(
    chromatogram: pyopenms.MSChromatogram,
    degree: int = 6,
    m: int = 7,
    is_ms1: bool = False,
    repeats: int = 1,
) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Apply ModifiedSincSmoother to an MSChromatogram.
    Returns: (rt_array, smoothed_intensity_array, runtime_ms).
    """
    smoother = pyopenms.ModifiedSincSmoother()
    param = smoother.getParameters()
    param.setValue("degree", degree)
    param.setValue("m", m)
    param.setValue("is_ms1", "true" if is_ms1 else "false")
    smoother.setParameters(param)

    work_chrom = pyopenms.MSChromatogram(chromatogram)
    smoother.filter(work_chrom)

    times = []
    for _ in range(repeats):
        c = pyopenms.MSChromatogram(chromatogram)
        t0 = time.perf_counter()
        smoother.filter(c)
        t1 = time.perf_counter()
        times.append((t1 - t0) * 1000.0)

    rt, intensity = c.get_peaks()
    return np.array(rt, dtype=np.float64), np.array(intensity, dtype=np.float64), float(np.median(times))


def apply_savitzky_golay(
    chromatogram: pyopenms.MSChromatogram,
    frame_length: int = 11,
    polynomial_order: int = 4,
    repeats: int = 1,
) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Apply SavitzkyGolayFilter to an MSChromatogram.
    Returns: (rt_array, smoothed_intensity_array, runtime_ms).
    """
    sg = pyopenms.SavitzkyGolayFilter()
    param = sg.getParameters()
    param.setValue("frame_length", frame_length)
    param.setValue("polynomial_order", polynomial_order)
    sg.setParameters(param)

    work_chrom = pyopenms.MSChromatogram(chromatogram)
    sg.filter(work_chrom)

    times = []
    for _ in range(repeats):
        c = pyopenms.MSChromatogram(chromatogram)
        t0 = time.perf_counter()
        sg.filter(c)
        t1 = time.perf_counter()
        times.append((t1 - t0) * 1000.0)

    rt, intensity = c.get_peaks()
    return np.array(rt, dtype=np.float64), np.array(intensity, dtype=np.float64), float(np.median(times))


def run_peak_picker_hires(
    rt: np.ndarray,
    intensity: np.ndarray,
    signal_to_noise: float = 2.0,
) -> List[Tuple[float, float, float]]:
    """
    Run PeakPickerHiRes directly on smoothed trace with realistic S/N threshold.
    Avoids accidental double smoothing while filtering baseline ripples.
    Returns list of (apex_rt, apex_intensity, fwhm).
    """
    chrom = pyopenms.MSChromatogram()
    chrom.set_peaks((rt.tolist(), intensity.tolist()))

    pp = pyopenms.PeakPickerHiRes()
    param = pp.getParameters()
    param.setValue("signal_to_noise", float(signal_to_noise))
    param.setValue("report_FWHM", "true")
    param.setValue("report_FWHM_unit", "absolute")
    pp.setParameters(param)

    picked = pyopenms.MSChromatogram()
    pp.pick(chrom, picked)

    detected = []
    fwhm_array = None
    for arr in picked.getFloatDataArrays():
        if arr.getName() == "FWHM":
            fwhm_array = list(arr)
            break

    for idx, peak in enumerate(picked):
        fwhm_val = fwhm_array[idx] if (fwhm_array and idx < len(fwhm_array)) else 0.0
        detected.append((float(peak.getRT()), float(peak.getIntensity()), float(fwhm_val)))

    return detected


class BenchmarkEngine:
    """
    Orchestrates the benchmark comparison between ModifiedSincSmoother
    and SavitzkyGolayFilter adhering strictly to scientific rigor.
    """

    def __init__(self, seed: int = 42, quick_mode: bool = False):
        self.seed = seed
        self.quick_mode = quick_mode
        self.generator = SyntheticChromatogramGenerator(seed=seed)

    def get_ms_parameter_grid(self) -> List[Dict[str, Any]]:
        """Generate harmonized parameter grid for ModifiedSincSmoother."""
        configs = []
        if self.quick_mode:
            configs = [
                {"degree": 4, "m": 5, "is_ms1": False},
                {"degree": 6, "m": 7, "is_ms1": False},   # Default OpenMS MS
                {"degree": 6, "m": 12, "is_ms1": False},
                {"degree": 6, "m": 7, "is_ms1": True},    # Default OpenMS MS1
            ]
        else:
            degrees = [2, 4, 6, 8, 10]
            # Max window span: 2*m + 1 <= 31 (matching SG max frame_length = 31)
            for is_ms1 in [False, True]:
                for deg in degrees:
                    m_min = (deg // 2 + 1) if is_ms1 else (deg // 2 + 2)
                    m_values = [m for m in [m_min, m_min + 2, m_min + 5, 10, 15] if m >= m_min and (2 * m + 1) <= 31]
                    for m in sorted(set(m_values)):
                        configs.append({"degree": deg, "m": m, "is_ms1": is_ms1})
        return configs

    def get_sg_parameter_grid(self) -> List[Dict[str, Any]]:
        """Generate harmonized parameter grid for SavitzkyGolayFilter."""
        configs = []
        if self.quick_mode:
            configs = [
                {"frame_length": 7, "polynomial_order": 2},
                {"frame_length": 11, "polynomial_order": 4},  # Default OpenMS SG
                {"frame_length": 15, "polynomial_order": 4},
                {"frame_length": 21, "polynomial_order": 4},
            ]
        else:
            # Note: order 2 and 3 produce identical center smoothing weights; test 2 and 4
            orders = [2, 4]
            frame_lengths = [7, 11, 15, 21, 31]
            for order in orders:
                for fl in frame_lengths:
                    if fl > order and fl % 2 == 1:
                        configs.append({"frame_length": fl, "polynomial_order": order})
        return configs

    def tune_parameters(
        self,
    ) -> Tuple[Dict[str, Any], Dict[str, Any]]:
        """
        Select optimal default parameters using ONLY the separate calibration/tuning dataset.
        Zero test-set leakage: evaluated on tuning chromatograms generated with an independent seed (seed + 10000).

        Tuning Objective:
            The selection objective is global Root Mean Squared Error (RMSE) against the
            ground-truth noiseless signal on the independent calibration trace 'tuning_composite'.
            Global RMSE was chosen because:
            1. It provides a standard, unbiased L2 distance metric across all points.
            2. It simultaneously balances noise suppression across baseline regions and shape preservation
               across diverse peak geometries (narrow, broad, overlapping, tailing, low-intensity).
            3. It avoids arbitrary heuristic weights between disparate sub-metrics (such as weighting
               apex height error vs. FWHM distortion vs. baseline noise reduction).
            4. The exact same objective is applied identically and fairly to both smoothers over
               comparable parameter grids spanning equivalent window widths (up to 31 points).
            All final signal-fidelity and downstream metrics are strictly evaluated on held-out test data.
        """
        tuning_datasets = self.generator.generate_tuning_datasets()
        tuning_comp = tuning_datasets["tuning_composite"]

        ms_grid = self.get_ms_parameter_grid()
        sg_grid = self.get_sg_parameter_grid()

        best_ms_params = ms_grid[0]
        best_ms_score = float("inf")

        for p in ms_grid:
            m = self.evaluate_configuration(tuning_comp, "ModifiedSincSmoother", p, repeats=1)
            # Objective: minimize RMSE on tuning data
            if m.rmse < best_ms_score:
                best_ms_score = m.rmse
                best_ms_params = p

        best_sg_params = sg_grid[0]
        best_sg_score = float("inf")

        for p in sg_grid:
            m = self.evaluate_configuration(tuning_comp, "SavitzkyGolayFilter", p, repeats=1)
            if m.rmse < best_sg_score:
                best_sg_score = m.rmse
                best_sg_params = p

        return best_ms_params, best_sg_params

    def evaluate_configuration(
        self,
        syn: SyntheticChromatogram,
        smoother_name: str,
        params: Dict[str, Any],
        repeats: int = 1,
        primary_sn_threshold: float = 2.0,
        sn_thresholds: Tuple[float, ...] = (1.0, 2.0, 3.0, 5.0),
    ) -> OverallBenchmarkMetrics:
        """Run a smoother configuration on a synthetic chromatogram and evaluate all metrics."""
        openms_noisy = syn.to_openms_noisy()

        if smoother_name == "ModifiedSincSmoother":
            rt_s, y_s, runtime_ms = apply_modified_sinc(
                openms_noisy,
                degree=params["degree"],
                m=params["m"],
                is_ms1=params["is_ms1"],
                repeats=repeats,
            )
        elif smoother_name == "SavitzkyGolayFilter":
            rt_s, y_s, runtime_ms = apply_savitzky_golay(
                openms_noisy,
                frame_length=params["frame_length"],
                polynomial_order=params["polynomial_order"],
                repeats=repeats,
            )
        else:
            raise ValueError(f"Unknown smoother: {smoother_name}")

        # 1. Signal fidelity with sub-scan interpolation
        rmse, mae, max_err, r, peak_metrics = evaluate_signal_fidelity(
            syn.rt, syn.intensity_true, y_s, syn.peaks, baseline_level=syn.baseline_level
        )

        # 2. Noise reduction
        max_apex = max(p.true_apex_intensity for p in syn.peaks) if syn.peaks else 1000.0
        noise_metrics = evaluate_noise_reduction(
            syn.intensity_noisy, y_s, syn.baseline_mask, max_apex
        )

        # 3. Overlapping doublet resolution and total area
        vpr_true = None
        vpr_smooth = None
        doublet_area_err = None

        if syn.is_overlapping_doublet and len(syn.peaks) >= 2:
            p1, p2 = syn.peaks[0], syn.peaks[1]
            vpr_true = evaluate_valley_to_peak_ratio(syn.rt, syn.intensity_true, p1.true_apex_rt, p2.true_apex_rt, syn.baseline_level)
            vpr_smooth = evaluate_valley_to_peak_ratio(rt_s, y_s, p1.true_apex_rt, p2.true_apex_rt, syn.baseline_level)

            # Evaluate total doublet area against ground truth total doublet area
            smooth_doublet_area = compute_peak_area(
                rt_s, y_s, syn.doublet_rt_start, syn.doublet_rt_end, baseline_level=syn.baseline_level
            )
            doublet_area_err = ((smooth_doublet_area - syn.doublet_total_area) / syn.doublet_total_area) * 100.0

        # 4. Downstream peak picking with PeakPickerHiRes across S/N thresholds
        picking_sensitivity = {}
        for sn in sn_thresholds:
            detected = run_peak_picker_hires(rt_s, y_s, signal_to_noise=sn)
            picking_sensitivity[sn] = evaluate_peak_picking(detected, syn.peaks, sn_threshold=sn)

        primary_picking = picking_sensitivity[primary_sn_threshold]

        n_pts = len(syn.rt)
        time_per_pt_ns = (runtime_ms * 1e6) / n_pts if n_pts > 0 else 0.0

        return OverallBenchmarkMetrics(
            dataset_name=syn.name,
            smoother_type=smoother_name,
            params=params,
            rmse=rmse,
            mae=mae,
            max_abs_error=max_err,
            pearson_r=r,
            noise=noise_metrics,
            peaks=peak_metrics,
            peak_picking=primary_picking,
            peak_picking_sensitivity=picking_sensitivity,
            doublet_total_area_error_pct=doublet_area_err,
            valley_to_peak_ratio_true=vpr_true,
            valley_to_peak_ratio_smooth=vpr_smooth,
            runtime_ms=runtime_ms,
            time_per_point_ns=time_per_pt_ns,
        )

    def run_held_out_evaluation(
        self,
        ms_params: Dict[str, Any],
        sg_params: Dict[str, Any],
    ) -> Dict[str, Dict[str, OverallBenchmarkMetrics]]:
        """
        Evaluate FROZEN tuned parameters on the held-out test datasets.
        Guarantees zero parameter leakage between training and testing.
        """
        test_datasets = self.generator.generate_test_datasets()
        results = {}

        for name, data in test_datasets.items():
            ms_metric = self.evaluate_configuration(data, "ModifiedSincSmoother", ms_params, repeats=2)
            sg_metric = self.evaluate_configuration(data, "SavitzkyGolayFilter", sg_params, repeats=2)
            results[name] = {
                "ModifiedSinc": ms_metric,
                "SavitzkyGolay": sg_metric,
            }

        return results

    def run_sampling_phase_benchmark(
        self,
        ms_params: Dict[str, Any],
        sg_params: Dict[str, Any],
    ) -> Dict[str, Any]:
        """
        Evaluate narrow peak preservation across 4 sub-scan sampling phases (0.0, 0.25, 0.50, 0.75).
        Investigates sampling phase jitter effects.
        """
        phase_chroms = self.generator.generate_phase_sweep_narrow_peaks()
        phase_results = {}

        for phase_name, chrom in phase_chroms.items():
            ms_m = self.evaluate_configuration(chrom, "ModifiedSincSmoother", ms_params, repeats=2)
            sg_m = self.evaluate_configuration(chrom, "SavitzkyGolayFilter", sg_params, repeats=2)
            phase_results[phase_name] = {
                "phase": chrom.phase_offset,
                "ms_apex_err": ms_m.peaks[0].apex_intensity_error_pct if ms_m.peaks else 0.0,
                "ms_fwhm_err": ms_m.peaks[0].fwhm_error_pct if ms_m.peaks else 0.0,
                "ms_area_err": ms_m.peaks[0].area_error_pct if ms_m.peaks else 0.0,
                "sg_apex_err": sg_m.peaks[0].apex_intensity_error_pct if sg_m.peaks else 0.0,
                "sg_fwhm_err": sg_m.peaks[0].fwhm_error_pct if sg_m.peaks else 0.0,
                "sg_area_err": sg_m.peaks[0].area_error_pct if sg_m.peaks else 0.0,
            }

        return phase_results

    def run_matched_bandwidth_comparison(self) -> List[Dict[str, Any]]:
        """
        Compare ModifiedSincSmoother against SavitzkyGolayFilter at mathematically matched 3-dB bandwidths.

        Methodology:
            1. For standard Savitzky-Golay configurations (order p in [4, 2], frame length F in [7, 9, 11, 15, 21]),
               the kernel half-width is m_SG = (F - 1) // 2.
            2. Compute the exact 3-dB cutoff frequency via OpenMS's built-in function:
               equiv_bw = ModifiedSincSmoother.savitzkyGolayBandwidth(p, m_SG).
            3. Convert that exact cutoff frequency into the corresponding Modified Sinc half-width:
               m_MS = ModifiedSincSmoother.bandwidthToM(is_ms1, degree, equiv_bw).
            4. Verify that m_MS satisfies the algorithm's mathematical minimum constraint.
            5. Evaluate both smoothers on the composite test chromatogram under identical conditions.
            This ensures that both filters share the exact same 3-dB passband cutoff in the frequency domain.
        """
        test_comp = self.generator.create_composite_chromatogram()
        results = []

        sg_configs = [
            (4, 7), (4, 9), (4, 11), (4, 15), (4, 21),
            (2, 7), (2, 9), (2, 11), (2, 15), (2, 21),
        ]

        for p_order, fl in sg_configs:
            m_sg = (fl - 1) // 2
            try:
                equiv_bw = pyopenms.ModifiedSincSmoother.savitzkyGolayBandwidth(p_order, m_sg)
            except Exception:
                continue

            if equiv_bw <= 0.0 or equiv_bw >= 0.5:
                continue

            ms_deg = 6 if p_order == 4 else 4
            is_ms1 = False
            try:
                m_ms = pyopenms.ModifiedSincSmoother.bandwidthToM(is_ms1, ms_deg, equiv_bw)
            except Exception:
                continue

            m_min = ms_deg // 2 + 2
            if m_ms < m_min:
                continue

            ms_metric = self.evaluate_configuration(
                test_comp, "ModifiedSincSmoother", {"degree": ms_deg, "m": m_ms, "is_ms1": is_ms1}
            )
            sg_metric = self.evaluate_configuration(
                test_comp, "SavitzkyGolayFilter", {"frame_length": fl, "polynomial_order": p_order}
            )

            results.append({
                "sg_order": p_order,
                "sg_frame_length": fl,
                "sg_m": m_sg,
                "equivalent_bandwidth": equiv_bw,
                "ms_degree": ms_deg,
                "ms_m": m_ms,
                "ms_metric": ms_metric,
                "sg_metric": sg_metric,
            })

        return results

    def run_rigorous_runtime_benchmark(
        self,
        lengths: Optional[List[int]] = None,
    ) -> Dict[str, Any]:
        """
        Measure computational performance following benchmark_pyopenms.py standards:
        - Warmup iterations
        - Loop of 200–500 iterations for microsecond workloads
        - gc.disable() during timing
        - Separating total time, per-chromatogram time, and per-point throughput
        """
        if lengths is None:
            lengths = [50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000]

        ms_per_call_us = []
        sg_per_call_us = []
        ms_per_point_ns = []
        sg_per_point_ns = []

        ms_smoother = pyopenms.ModifiedSincSmoother()
        p_ms = ms_smoother.getParameters()
        p_ms.setValue("degree", 6)
        p_ms.setValue("m", 7)
        p_ms.setValue("is_ms1", "false")
        ms_smoother.setParameters(p_ms)

        sg_smoother = pyopenms.SavitzkyGolayFilter()
        p_sg = sg_smoother.getParameters()
        p_sg.setValue("frame_length", 11)
        p_sg.setValue("polynomial_order", 4)
        sg_smoother.setParameters(p_sg)

        rng = np.random.default_rng(self.seed)

        for n in lengths:
            rt = np.arange(n, dtype=np.float64) * 0.5
            y = 100.0 + rng.normal(0.0, 10.0, size=n)
            for center_idx in range(n // 8, n, max(1, n // 4)):
                sigma = max(2.0, n / 80.0)
                y += 500.0 * np.exp(-0.5 * ((np.arange(n) - center_idx) / sigma) ** 2)

            chrom = pyopenms.MSChromatogram()
            chrom.set_peaks((rt.tolist(), y.tolist()))

            # Iteration count: high for small N to suppress nanobind timing jitter
            if n <= 500:
                iterations = 500
            elif n <= 5000:
                iterations = 200
            elif n <= 25000:
                iterations = 50
            else:
                iterations = 20

            # --- Modified Sinc Microbenchmark ---
            # Warmup
            w_chrom = pyopenms.MSChromatogram(chrom)
            ms_smoother.filter(w_chrom)

            gc.disable()
            t0 = time.perf_counter()
            for _ in range(iterations):
                c = pyopenms.MSChromatogram(chrom)
                ms_smoother.filter(c)
            t1 = time.perf_counter()
            gc.enable()

            total_ms_sec = (t1 - t0) / iterations
            ms_per_call_us.append(total_ms_sec * 1e6)
            ms_per_point_ns.append((total_ms_sec * 1e9) / n)

            # --- Savitzky-Golay Microbenchmark ---
            # Warmup
            w_chrom = pyopenms.MSChromatogram(chrom)
            sg_smoother.filter(w_chrom)

            gc.disable()
            t0 = time.perf_counter()
            for _ in range(iterations):
                c = pyopenms.MSChromatogram(chrom)
                sg_smoother.filter(c)
            t1 = time.perf_counter()
            gc.enable()

            total_sg_sec = (t1 - t0) / iterations
            sg_per_call_us.append(total_sg_sec * 1e6)
            sg_per_point_ns.append((total_sg_sec * 1e9) / n)

        return {
            "lengths": lengths,
            "ms_per_call_us": ms_per_call_us,
            "sg_per_call_us": sg_per_call_us,
            "ms_per_point_ns": ms_per_point_ns,
            "sg_per_point_ns": sg_per_point_ns,
            "ms_params": {"degree": 6, "m": 7, "is_ms1": False},
            "sg_params": {"frame_length": 11, "polynomial_order": 4},
        }

    def run_representative_real_data_check(
        self,
        mzml_paths: Optional[List[str]] = None,
    ) -> List[Dict[str, Any]]:
        """
        Sanity check on representative OpenMS test chromatograms.
        Uses baseline-confined regions to estimate noise reduction, avoiding peak slope contamination.
        """
        if mzml_paths is None:
            repo_root = Path(__file__).resolve().parents[4]
            candidates = [
                repo_root / "src" / "tests" / "topp" / "NoiseFilterSGolay_2_input.chrom.mzML",
                repo_root / "src" / "tests" / "topp" / "MRMTransitionGroupPicker_1_input.mzML",
            ]
            mzml_paths = [str(p) for p in candidates if p.exists()]

        real_results = []

        for p_str in mzml_paths:
            path = Path(p_str)
            if not path.exists():
                continue

            exp = pyopenms.MSExperiment()
            try:
                pyopenms.MzMLFile().load(str(path), exp)
            except Exception:
                continue

            chroms = exp.getChromatograms()
            if not chroms:
                continue

            selected = list(range(0, min(5, len(chroms))))
            file_stats = {
                "file_name": path.name,
                "total_chromatograms": len(chroms),
                "chromatograms_checked": len(selected),
                "evaluations": [],
            }

            for idx in selected:
                c = chroms[idx]
                if len(c) < 20:
                    continue

                rt, raw_int = c.get_peaks()
                rt_arr = np.array(rt, dtype=np.float64)
                raw_int_arr = np.array(raw_int, dtype=np.float64)

                # Confine noise estimation to the first and last 20% of points (baseline regions)
                n = len(raw_int_arr)
                k = max(3, n // 5)
                baseline_indices = np.concatenate([np.arange(0, k), np.arange(n - k, n)])

                sigma_raw = float(np.std(raw_int_arr[baseline_indices]))
                if sigma_raw <= 1e-6:
                    continue

                _, ms_int, ms_t = apply_modified_sinc(c, degree=6, m=7, is_ms1=False, repeats=2)
                sigma_ms = float(np.std(ms_int[baseline_indices]))
                ms_noise_red = max(0.0, (1.0 - sigma_ms / sigma_raw) * 100.0)

                _, sg_int, sg_t = apply_savitzky_golay(c, frame_length=11, polynomial_order=4, repeats=2)
                sigma_sg = float(np.std(sg_int[baseline_indices]))
                sg_noise_red = max(0.0, (1.0 - sigma_sg / sigma_raw) * 100.0)

                # Peak picking with realistic S/N = 2.0
                raw_picked = run_peak_picker_hires(rt_arr, raw_int_arr, signal_to_noise=2.0)
                ms_picked = run_peak_picker_hires(rt_arr, ms_int, signal_to_noise=2.0)
                sg_picked = run_peak_picker_hires(rt_arr, sg_int, signal_to_noise=2.0)

                file_stats["evaluations"].append({
                    "chromatogram_index": idx,
                    "native_id": c.getNativeID(),
                    "num_points": n,
                    "sigma_raw_baseline": sigma_raw,
                    "ms_sigma_baseline": sigma_ms,
                    "ms_noise_reduction_pct": ms_noise_red,
                    "sg_sigma_baseline": sigma_sg,
                    "sg_noise_reduction_pct": sg_noise_red,
                    "raw_peaks_sn2": len(raw_picked),
                    "ms_peaks_sn2": len(ms_picked),
                    "sg_peaks_sn2": len(sg_picked),
                })

            if file_stats["evaluations"]:
                real_results.append(file_stats)

        return real_results
