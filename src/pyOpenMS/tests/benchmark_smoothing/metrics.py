"""
Evaluation metrics module for OpenMS smoothing benchmark.
GitHub Issue #10425: Benchmark ModifiedSincSmoother against Savitzky–Golay.

Calculates:
- Sub-scan interpolated signal fidelity (parabolic apex RT & intensity, flank-interpolated FWHM)
- Total doublet area preservation for overlapping peaks (preventing artificial +34% error)
- Non-rectified baseline noise reduction and true SNR gain
- Downstream peak-picking evaluation with configurable S/N thresholds
- Overlapping peak valley-to-peak ratio (VPR)
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, asdict
from typing import List, Optional, Tuple, Dict, Any

import numpy as np

# Numerical trapezoid integration compatible with numpy 1.x and 2.x
trapz_fn = getattr(np, "trapezoid", getattr(np, "trapz", None))


@dataclass
class PeakFidelityMetric:
    """Fidelity metrics evaluated for an individual chromatographic peak using sub-scan interpolation."""
    peak_type: str
    true_rt: float
    true_intensity: float
    true_fwhm: float
    true_area: float
    smooth_apex_rt: float
    smooth_apex_intensity: float
    smooth_fwhm: float
    smooth_area: float
    apex_rt_shift: float             # delta RT = smooth - true
    apex_intensity_error_pct: float  # (smooth - true) / true * 100%
    area_error_pct: float            # (smooth - true) / true * 100%
    fwhm_error_pct: float            # (smooth - true) / true * 100%

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class NoiseMetrics:
    """Noise reduction and SNR metrics."""
    noise_sigma_raw: float
    noise_sigma_smooth: float
    noise_reduction_pct: float       # (1 - sigma_smooth / sigma_raw) * 100%
    snr_raw: float
    snr_smooth: float
    snr_gain_factor: float           # snr_smooth / snr_raw

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class PeakPickingMetrics:
    """Downstream peak picking performance."""
    sn_threshold: float
    total_true_peaks: int
    total_detected_peaks: int
    true_positives: int
    false_positives: int
    false_negatives: int
    precision: float
    recall: float
    f1_score: float
    mean_rt_error: float
    mean_fwhm_error: float

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class OverallBenchmarkMetrics:
    """Complete summary of metrics for a single smoother configuration on a chromatogram."""
    dataset_name: str
    smoother_type: str
    params: Dict[str, Any]
    rmse: float
    mae: float
    max_abs_error: float
    pearson_r: float
    noise: Optional[NoiseMetrics] = None
    peaks: List[PeakFidelityMetric] = field(default_factory=list)
    peak_picking: Optional[PeakPickingMetrics] = None
    peak_picking_sensitivity: Dict[float, PeakPickingMetrics] = field(default_factory=dict)
    doublet_total_area_error_pct: Optional[float] = None
    valley_to_peak_ratio_true: Optional[float] = None
    valley_to_peak_ratio_smooth: Optional[float] = None
    runtime_ms: float = 0.0
    time_per_point_ns: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        if "peak_picking_sensitivity" in d:
            d["peak_picking_sensitivity"] = {
                str(k): v for k, v in d["peak_picking_sensitivity"].items()
            }
        return d


def parabolic_apex_interpolation(sub_rt: np.ndarray, sub_int: np.ndarray) -> Tuple[float, float]:
    """
    Sub-scan 3-point parabolic interpolation around the maximum.
    Avoids discrete grid quantization bias on narrow peaks.
    """
    max_idx = int(np.argmax(sub_int))
    if max_idx == 0 or max_idx == len(sub_int) - 1:
        return float(sub_rt[max_idx]), float(sub_int[max_idx])

    alpha = float(sub_int[max_idx - 1])
    beta = float(sub_int[max_idx])
    gamma = float(sub_int[max_idx + 1])

    denom = alpha - 2.0 * beta + gamma
    if abs(denom) < 1e-12:
        return float(sub_rt[max_idx]), float(sub_int[max_idx])

    p = 0.5 * (alpha - gamma) / denom
    p = float(np.clip(p, -0.5, 0.5))

    dt = float(sub_rt[max_idx] - sub_rt[max_idx - 1])
    apex_rt = float(sub_rt[max_idx] + p * dt)
    apex_int = float(beta - 0.25 * (alpha - gamma) * p)
    return apex_rt, apex_int


def interpolated_fwhm(
    sub_rt: np.ndarray,
    sub_int: np.ndarray,
    apex_int: float,
    baseline: float = 0.0,
) -> float:
    """
    Continuous linear flank interpolation for FWHM at half-maximum.
    Avoids discrete integer index width quantization.
    """
    if len(sub_rt) < 3 or apex_int <= baseline:
        return 0.0

    half_max = baseline + (apex_int - baseline) / 2.0
    max_idx = int(np.argmax(sub_int))

    # Left flank crossing
    left_rt = float(sub_rt[0])
    for i in range(max_idx, 0, -1):
        if sub_int[i - 1] <= half_max <= sub_int[i]:
            denom = float(sub_int[i] - sub_int[i - 1])
            frac = float(half_max - sub_int[i - 1]) / denom if abs(denom) > 1e-12 else 0.0
            left_rt = float(sub_rt[i - 1] + frac * (sub_rt[i] - sub_rt[i - 1]))
            break

    # Right flank crossing
    right_rt = float(sub_rt[-1])
    for i in range(max_idx, len(sub_int) - 1):
        if sub_int[i] >= half_max >= sub_int[i + 1]:
            denom = float(sub_int[i] - sub_int[i + 1])
            frac = float(sub_int[i] - half_max) / denom if abs(denom) > 1e-12 else 0.0
            right_rt = float(sub_rt[i] + frac * (sub_rt[i + 1] - sub_rt[i]))
            break

    return max(0.0, float(right_rt - left_rt))


def compute_fwhm_and_apex(
    rt: np.ndarray,
    intensity: np.ndarray,
    rt_start: float,
    rt_end: float,
    baseline_level: float = 0.0,
) -> Tuple[float, float, float]:
    """
    Compute peak apex RT, apex intensity, and FWHM within the [rt_start, rt_end] window
    using sub-scan parabolic and flank interpolation.
    """
    window = (rt >= rt_start) & (rt <= rt_end)
    if not np.any(window):
        return float(rt[0]), 0.0, 0.0

    sub_rt = rt[window]
    sub_int = intensity[window]

    apex_rt, apex_int = parabolic_apex_interpolation(sub_rt, sub_int)
    fwhm = interpolated_fwhm(sub_rt, sub_int, apex_int, baseline=baseline_level)

    # Return baseline-subtracted apex height
    net_apex_height = max(0.0, apex_int - baseline_level)
    return apex_rt, net_apex_height, fwhm


def compute_peak_area(
    rt: np.ndarray,
    intensity: np.ndarray,
    rt_start: float,
    rt_end: float,
    baseline_level: float = 0.0,
) -> float:
    """Compute baseline-corrected peak area within [rt_start, rt_end] using trapezoidal integration."""
    window = (rt >= rt_start) & (rt <= rt_end)
    if not np.any(window):
        return 0.0
    sub_rt = rt[window]
    sub_int = np.maximum(0.0, intensity[window] - baseline_level)
    return float(trapz_fn(sub_int, sub_rt))


def evaluate_signal_fidelity(
    rt: np.ndarray,
    y_true: np.ndarray,
    y_smooth: np.ndarray,
    peak_metadata_list: List[Any],
    baseline_level: float = 0.0,
) -> Tuple[float, float, float, float, List[PeakFidelityMetric]]:
    """Compute global fidelity (RMSE, MAE, MaxErr, Pearson R) and per-peak metrics."""
    diff = y_smooth - y_true
    rmse = float(np.sqrt(np.mean(diff ** 2)))
    mae = float(np.mean(np.abs(diff)))
    max_abs_err = float(np.max(np.abs(diff)))

    std_true = float(np.std(y_true))
    std_smooth = float(np.std(y_smooth))
    if std_true > 1e-12 and std_smooth > 1e-12:
        r = float(np.corrcoef(y_true, y_smooth)[0, 1])
    else:
        r = 1.0 if rmse < 1e-6 else 0.0

    peak_metrics: List[PeakFidelityMetric] = []
    for meta in peak_metadata_list:
        apex_rt_s, apex_int_s, fwhm_s = compute_fwhm_and_apex(
            rt, y_smooth, meta.rt_start, meta.rt_end, baseline_level=baseline_level
        )
        area_s = compute_peak_area(rt, y_smooth, meta.rt_start, meta.rt_end, baseline_level=baseline_level)

        rt_shift = apex_rt_s - meta.true_apex_rt
        apex_int_err = ((apex_int_s - meta.true_apex_intensity) / meta.true_apex_intensity) * 100.0 if meta.true_apex_intensity > 0 else 0.0
        area_err = ((area_s - meta.true_area) / meta.true_area) * 100.0 if meta.true_area > 0 else 0.0
        fwhm_err = ((fwhm_s - meta.true_fwhm) / meta.true_fwhm) * 100.0 if meta.true_fwhm > 0 else 0.0

        peak_metrics.append(PeakFidelityMetric(
            peak_type=meta.peak_type,
            true_rt=meta.true_apex_rt,
            true_intensity=meta.true_apex_intensity,
            true_fwhm=meta.true_fwhm,
            true_area=meta.true_area,
            smooth_apex_rt=apex_rt_s,
            smooth_apex_intensity=apex_int_s,
            smooth_fwhm=fwhm_s,
            smooth_area=area_s,
            apex_rt_shift=rt_shift,
            apex_intensity_error_pct=apex_int_err,
            area_error_pct=area_err,
            fwhm_error_pct=fwhm_err,
        ))

    return rmse, mae, max_abs_err, r, peak_metrics


def evaluate_noise_reduction(
    y_raw: np.ndarray,
    y_smooth: np.ndarray,
    baseline_mask: Optional[np.ndarray],
    signal_apex_height: float,
) -> Optional[NoiseMetrics]:
    """
    Measure noise standard deviation before and after smoothing in baseline regions.
    Uses non-rectified baseline samples.
    """
    if baseline_mask is not None and np.sum(baseline_mask) >= 10:
        base_raw = y_raw[baseline_mask]
        base_smooth = y_smooth[baseline_mask]
        sigma_raw = float(np.std(base_raw))
        sigma_smooth = float(np.std(base_smooth))
    else:
        diff_raw = np.diff(y_raw) / math.sqrt(2.0)
        diff_smooth = np.diff(y_smooth) / math.sqrt(2.0)
        sigma_raw = float(np.std(diff_raw))
        sigma_smooth = float(np.std(diff_smooth))

    if sigma_raw <= 1e-9:
        return None

    noise_red_pct = max(0.0, (1.0 - sigma_smooth / sigma_raw) * 100.0)
    snr_raw = signal_apex_height / sigma_raw if sigma_raw > 0 else float("inf")
    snr_smooth = signal_apex_height / sigma_smooth if sigma_smooth > 0 else float("inf")
    snr_gain = snr_smooth / snr_raw if snr_raw > 0 and not math.isinf(snr_raw) else 1.0

    return NoiseMetrics(
        noise_sigma_raw=sigma_raw,
        noise_sigma_smooth=sigma_smooth,
        noise_reduction_pct=noise_red_pct,
        snr_raw=snr_raw,
        snr_smooth=snr_smooth,
        snr_gain_factor=snr_gain,
    )


def evaluate_valley_to_peak_ratio(
    rt: np.ndarray,
    intensity: np.ndarray,
    peak1_rt: float,
    peak2_rt: float,
    baseline_level: float = 0.0,
) -> Optional[float]:
    """
    Compute valley-to-peak ratio for overlapping doublet peaks:
    VPR = (I_valley - baseline) / min(I_peak1 - baseline, I_peak2 - baseline).
    """
    rt_lo, rt_hi = min(peak1_rt, peak2_rt), max(peak1_rt, peak2_rt)
    window = (rt >= rt_lo) & (rt <= rt_hi)
    if np.sum(window) < 3:
        return None

    sub_int = intensity[window] - baseline_level
    valley_val = float(np.min(sub_int))

    w1 = (rt >= rt_lo - 2.0) & (rt <= rt_lo + 2.0)
    w2 = (rt >= rt_hi - 2.0) & (rt <= rt_hi + 2.0)
    h1 = float(np.max(intensity[w1]) - baseline_level) if np.any(w1) else 1.0
    h2 = float(np.max(intensity[w2]) - baseline_level) if np.any(w2) else 1.0

    min_apex = max(1e-6, min(h1, h2))
    return float(valley_val / min_apex)


def evaluate_peak_picking(
    detected_peaks: List[Tuple[float, float, float]],
    true_peaks: List[Any],
    rt_tolerance_factor: float = 0.5,
    sn_threshold: float = 2.0,
) -> PeakPickingMetrics:
    """
    Compare detected chromatographic peaks against ground-truth peaks.
    Match condition: |RT_detected - RT_true| <= rt_tolerance_factor * FWHM_true.
    """
    matched_true = set()
    matched_det = set()
    rt_errors = []
    fwhm_errors = []

    for det_idx, (det_rt, det_int, det_fwhm) in enumerate(detected_peaks):
        best_true_idx = None
        best_dist = float("inf")

        for true_idx, true_p in enumerate(true_peaks):
            if true_idx in matched_true:
                continue
            tol = max(0.5, rt_tolerance_factor * true_p.true_fwhm)
            dist = abs(det_rt - true_p.true_apex_rt)
            if dist <= tol and dist < best_dist:
                best_dist = dist
                best_true_idx = true_idx

        if best_true_idx is not None:
            matched_true.add(best_true_idx)
            matched_det.add(det_idx)
            true_p = true_peaks[best_true_idx]
            rt_errors.append(abs(det_rt - true_p.true_apex_rt))
            if det_fwhm > 0 and true_p.true_fwhm > 0:
                fwhm_errors.append(abs(det_fwhm - true_p.true_fwhm) / true_p.true_fwhm * 100.0)

    tp = len(matched_det)
    fp = len(detected_peaks) - tp
    fn = len(true_peaks) - len(matched_true)

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = (2.0 * precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0

    mean_rt_err = float(np.mean(rt_errors)) if rt_errors else 0.0
    mean_fwhm_err = float(np.mean(fwhm_errors)) if fwhm_errors else 0.0

    return PeakPickingMetrics(
        sn_threshold=sn_threshold,
        total_true_peaks=len(true_peaks),
        total_detected_peaks=len(detected_peaks),
        true_positives=tp,
        false_positives=fp,
        false_negatives=fn,
        precision=precision,
        recall=recall,
        f1_score=f1,
        mean_rt_error=mean_rt_err,
        mean_fwhm_error=mean_fwhm_err,
    )
