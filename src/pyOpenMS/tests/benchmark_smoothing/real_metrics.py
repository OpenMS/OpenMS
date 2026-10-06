"""
Real-data empirical evaluation metrics module for OpenMS smoothing benchmark.
GitHub Issue #10425: Benchmark ModifiedSincSmoother against Savitzky–Golay on chromatograms.

CRITICAL SCIENTIFIC PRINCIPLES:
- No synthetic ground-truth RMSE is calculated on real chromatograms (no known 'true' signal exists).
- Signal preservation is assessed via relative distortion:
  * Apex RT shift (sub-scan parabolic interpolation)
  * Apex intensity change & ratio
  * Peak area ratio (numerical integration)
  * FWHM ratio & width change (continuous flank interpolation)
- Noise is estimated strictly on peak-free flank / baseline regions via Median Absolute Deviation (MAD)
  of first differences divided by sqrt(2), avoiding confounding peak slopes with noise.
- Overlapping peak resolution is assessed via local maxima counts, valley-to-peak ratio (VPR),
  and peak separation.
- Downstream detection is evaluated with PeakPickerHiRes at S/N = 2.0.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, asdict
from typing import List, Dict, Any, Optional, Tuple

import numpy as np
import pyopenms

try:
    from .metrics import parabolic_apex_interpolation, interpolated_fwhm, trapz_fn
except ImportError:
    from metrics import parabolic_apex_interpolation, interpolated_fwhm, trapz_fn


def estimate_flank_noise(
    y: np.ndarray,
    flank_fraction: float = 0.20,
) -> float:
    """
    Estimate baseline noise strictly from peak-free flank regions.
    Extracts the outer flank regions (first and last flank_fraction of points)
    and computes the robust MAD of first differences / sqrt(2).
    """
    n = len(y)
    if n < 6:
        return 0.0

    k = max(2, int(n * flank_fraction))
    flank_indices = np.concatenate([np.arange(0, k), np.arange(n - k, n)])
    flank_y = y[flank_indices]

    dy = np.diff(flank_y)
    if len(dy) < 2:
        return float(np.std(flank_y))

    med_dy = np.median(dy)
    mad = np.median(np.abs(dy - med_dy))
    scale = 0.6745 * math.sqrt(2.0)
    return float(mad / scale) if scale > 0 else 0.0


def compute_valley_to_peak_ratio_real(y: np.ndarray) -> Optional[float]:
    """
    Compute valley-to-peak ratio (VPR) if trace contains co-eluting doublet/peaks.
    VPR = valley_intensity / min(apex1_intensity, apex2_intensity).
    Returns None if trace does not exhibit multiple distinct peaks.
    """
    if len(y) < 7:
        return None

    # Detect local maxima
    peaks = []
    for i in range(1, len(y) - 1):
        if y[i] > y[i - 1] and y[i] >= y[i + 1]:
            peaks.append(i)

    if len(peaks) < 2:
        return None

    # Sort peaks by intensity and take top 2
    peaks.sort(key=lambda idx: y[idx], reverse=True)
    p1, p2 = sorted(peaks[:2])

    # Find minimum valley between p1 and p2
    if p2 - p1 < 2:
        return None

    valley_idx = p1 + int(np.argmin(y[p1:p2 + 1]))
    valley_val = float(y[valley_idx])
    min_apex = float(min(y[p1], y[p2]))

    baseline = float(np.percentile(y, 10))
    net_valley = max(0.0, valley_val - baseline)
    net_apex = max(1e-9, min_apex - baseline)

    return float(net_valley / net_apex)


@dataclass
class RealSmootherMetrics:
    """Metrics for a single smoother applied to a real chromatogram."""
    smoother_type: str
    params: Dict[str, Any]
    apex_rt: float
    apex_intensity: float
    apex_rt_shift: float             # RT_smooth - RT_raw (s)
    apex_intensity_ratio: float      # I_smooth / I_raw
    apex_intensity_change_pct: float # (I_smooth - I_raw) / I_raw * 100%
    peak_area: float
    peak_area_ratio: float           # Area_smooth / Area_raw
    fwhm: float
    fwhm_ratio: float                # FWHM_smooth / FWHM_raw
    fwhm_change_pct: float           # (FWHM_smooth - FWHM_raw) / FWHM_raw * 100%
    flank_noise_sigma: float
    noise_reduction_pct: float       # (1 - sigma_smooth / sigma_raw) * 100%
    snr: float
    snr_gain_factor: float           # SNR_smooth / SNR_raw
    valley_to_peak_ratio: Optional[float] = None
    num_peaks_detected_sn2: int = 0
    runtime_us: float = 0.0
    ns_per_point: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class RealTraceEvaluation:
    """Comparison of Modified Sinc and Savitzky-Golay on a single real chromatogram."""
    trace_id: str
    native_id: str
    source_file: str
    num_points: int
    raw_apex_rt: float
    raw_apex_intensity: float
    raw_peak_area: float
    raw_fwhm: float
    raw_flank_noise_sigma: float
    raw_snr: float
    raw_valley_to_peak_ratio: Optional[float]
    raw_peaks_detected_sn2: int
    intensity_stratum: str
    width_stratum: str
    modified_sinc: RealSmootherMetrics
    savitzky_golay: RealSmootherMetrics

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class RealDatasetBenchmarkResult:
    """Aggregated empirical benchmark results for a real LC-MS dataset."""
    dataset_accession: str
    dataset_name: str
    total_traces_evaluated: int
    ms_params: Dict[str, Any]
    sg_params: Dict[str, Any]
    traces: List[RealTraceEvaluation]

    # Summary Statistics (Mean ± Std, Median, and IQR = Q75 - Q25)
    # 1. Apex RT Shift (seconds)
    ms_mean_apex_shift: float = 0.0
    ms_median_apex_shift: float = 0.0
    ms_iqr_apex_shift: float = 0.0
    sg_mean_apex_shift: float = 0.0
    sg_median_apex_shift: float = 0.0
    sg_iqr_apex_shift: float = 0.0

    # 2. Apex Intensity Ratio (I_smooth / I_raw)
    ms_mean_apex_ratio: float = 0.0
    ms_median_apex_ratio: float = 0.0
    ms_iqr_apex_ratio: float = 0.0
    sg_mean_apex_ratio: float = 0.0
    sg_median_apex_ratio: float = 0.0
    sg_iqr_apex_ratio: float = 0.0

    # 3. Peak Area Ratio (Area_smooth / Area_raw)
    ms_mean_area_ratio: float = 0.0
    ms_median_area_ratio: float = 0.0
    ms_iqr_area_ratio: float = 0.0
    sg_mean_area_ratio: float = 0.0
    sg_median_area_ratio: float = 0.0
    sg_iqr_area_ratio: float = 0.0

    # 4. FWHM Ratio (FWHM_smooth / FWHM_raw)
    ms_mean_fwhm_ratio: float = 0.0
    ms_median_fwhm_ratio: float = 0.0
    ms_iqr_fwhm_ratio: float = 0.0
    sg_mean_fwhm_ratio: float = 0.0
    sg_median_fwhm_ratio: float = 0.0
    sg_iqr_fwhm_ratio: float = 0.0

    # 5. Baseline Noise Reduction (%)
    ms_mean_noise_red_pct: float = 0.0
    ms_median_noise_red_pct: float = 0.0
    ms_iqr_noise_red_pct: float = 0.0
    sg_mean_noise_red_pct: float = 0.0
    sg_median_noise_red_pct: float = 0.0
    sg_iqr_noise_red_pct: float = 0.0

    # 6. Peak Picking Agreement at S/N = 2.0
    mean_raw_peaks_sn2: float = 0.0
    mean_ms_peaks_sn2: float = 0.0
    mean_sg_peaks_sn2: float = 0.0
    peak_count_exact_match_pct: float = 0.0

    # 7. Smoothing Runtime (microseconds per call)
    ms_mean_runtime_us: float = 0.0
    sg_mean_runtime_us: float = 0.0
    speed_ratio_ms_sg: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        d["traces"] = [t.to_dict() for t in self.traces]
        return d


def evaluate_single_real_trace(
    chrom: pyopenms.MSChromatogram,
    trace_id: str,
    source_file: str,
    intensity_stratum: str,
    width_stratum: str,
    ms_params: Dict[str, Any],
    sg_params: Dict[str, Any],
    ms_smoother_instance: pyopenms.ModifiedSincSmoother,
    sg_smoother_instance: pyopenms.SavitzkyGolayFilter,
) -> RealTraceEvaluation:
    """
    Evaluate both smoothers on a single real chromatogram under identical raw inputs.
    Uses local peak window tracking around the raw apex to prevent peak jumping across distant baseline noise.
    """
    rt_raw, int_raw = chrom.get_peaks()
    rt_arr = np.array(rt_raw, dtype=np.float64)
    y_raw = np.array(int_raw, dtype=np.float64)
    n_pts = len(y_raw)

    # 1. Raw baseline metrics
    raw_apex_rt, raw_apex_int = parabolic_apex_interpolation(rt_arr, y_raw)
    raw_baseline = float(np.percentile(y_raw, 10))
    raw_fwhm = interpolated_fwhm(rt_arr, y_raw, raw_apex_int, raw_baseline)
    raw_area = float(trapz_fn(y_raw - raw_baseline, rt_arr))
    raw_noise = estimate_flank_noise(y_raw)
    raw_snr = max(0.0, (raw_apex_int - raw_baseline) / raw_noise) if raw_noise > 1e-9 else 0.0
    raw_vpr = compute_valley_to_peak_ratio_real(y_raw)
    raw_idx = int(np.argmin(np.abs(rt_arr - raw_apex_rt)))

    # Window around raw apex for consistent local peak tracking (prevents jump to distant noise)
    dt = float(np.median(np.diff(rt_arr))) if len(rt_arr) > 1 else 1.0
    pts_in_fwhm = max(2, int(raw_fwhm / dt)) if dt > 0 else 5
    k_win = max(8, int(pts_in_fwhm * 1.5))
    w_start = max(0, raw_idx - k_win)
    w_end = min(len(rt_arr), raw_idx + k_win + 1)

    # Downstream peak picking on raw
    raw_peaks = pyopenms.MSSpectrum()
    pp_raw = pyopenms.PeakPickerHiRes()
    pp_raw_params = pp_raw.getParameters()
    pp_raw_params.setValue("signal_to_noise", 2.0)
    pp_raw.setParameters(pp_raw_params)
    spec_raw = pyopenms.MSSpectrum()
    spec_raw.set_peaks((rt_arr.tolist(), y_raw.tolist()))
    pp_raw.pick(spec_raw, raw_peaks)
    raw_n_peaks = len(raw_peaks)

    # 2. Modified Sinc smoothing
    c_ms = pyopenms.MSChromatogram(chrom)
    ms_smoother_instance.filter(c_ms)
    _, int_ms = c_ms.get_peaks()
    y_ms = np.array(int_ms, dtype=np.float64)

    ms_apex_rt, ms_apex_int = parabolic_apex_interpolation(rt_arr[w_start:w_end], y_ms[w_start:w_end])
    ms_baseline = float(np.percentile(y_ms, 10))
    ms_fwhm = interpolated_fwhm(rt_arr, y_ms, ms_apex_int, ms_baseline)
    ms_area = float(trapz_fn(y_ms - ms_baseline, rt_arr))
    ms_noise = estimate_flank_noise(y_ms)
    ms_snr = max(0.0, (ms_apex_int - ms_baseline) / ms_noise) if ms_noise > 1e-9 else 0.0
    ms_noise_red = max(0.0, (1.0 - ms_noise / raw_noise) * 100.0) if raw_noise > 1e-9 else 0.0
    ms_vpr = compute_valley_to_peak_ratio_real(y_ms)

    spec_ms = pyopenms.MSSpectrum()
    spec_ms.set_peaks((rt_arr.tolist(), y_ms.tolist()))
    peaks_ms = pyopenms.MSSpectrum()
    pp_raw.pick(spec_ms, peaks_ms)
    ms_n_peaks = len(peaks_ms)

    ms_metrics = RealSmootherMetrics(
        smoother_type="ModifiedSincSmoother",
        params=ms_params,
        apex_rt=ms_apex_rt,
        apex_intensity=ms_apex_int,
        apex_rt_shift=ms_apex_rt - raw_apex_rt,
        apex_intensity_ratio=ms_apex_int / raw_apex_int if raw_apex_int > 0 else 1.0,
        apex_intensity_change_pct=((ms_apex_int - raw_apex_int) / raw_apex_int * 100.0) if raw_apex_int > 0 else 0.0,
        peak_area=ms_area,
        peak_area_ratio=ms_area / raw_area if raw_area > 0 else 1.0,
        fwhm=ms_fwhm,
        fwhm_ratio=ms_fwhm / raw_fwhm if raw_fwhm > 0 else 1.0,
        fwhm_change_pct=((ms_fwhm - raw_fwhm) / raw_fwhm * 100.0) if raw_fwhm > 0 else 0.0,
        flank_noise_sigma=ms_noise,
        noise_reduction_pct=ms_noise_red,
        snr=ms_snr,
        snr_gain_factor=ms_snr / raw_snr if raw_snr > 0 else 1.0,
        valley_to_peak_ratio=ms_vpr,
        num_peaks_detected_sn2=ms_n_peaks,
    )

    # 3. Savitzky-Golay smoothing
    c_sg = pyopenms.MSChromatogram(chrom)
    sg_smoother_instance.filter(c_sg)
    _, int_sg = c_sg.get_peaks()
    y_sg = np.array(int_sg, dtype=np.float64)

    sg_apex_rt, sg_apex_int = parabolic_apex_interpolation(rt_arr[w_start:w_end], y_sg[w_start:w_end])
    sg_baseline = float(np.percentile(y_sg, 10))
    sg_fwhm = interpolated_fwhm(rt_arr, y_sg, sg_apex_int, sg_baseline)
    sg_area = float(trapz_fn(y_sg - sg_baseline, rt_arr))
    sg_noise = estimate_flank_noise(y_sg)
    sg_snr = max(0.0, (sg_apex_int - sg_baseline) / sg_noise) if sg_noise > 1e-9 else 0.0
    sg_noise_red = max(0.0, (1.0 - sg_noise / raw_noise) * 100.0) if raw_noise > 1e-9 else 0.0
    sg_vpr = compute_valley_to_peak_ratio_real(y_sg)

    spec_sg = pyopenms.MSSpectrum()
    spec_sg.set_peaks((rt_arr.tolist(), y_sg.tolist()))
    peaks_sg = pyopenms.MSSpectrum()
    pp_raw.pick(spec_sg, peaks_sg)
    sg_n_peaks = len(peaks_sg)

    sg_metrics = RealSmootherMetrics(
        smoother_type="SavitzkyGolayFilter",
        params=sg_params,
        apex_rt=sg_apex_rt,
        apex_intensity=sg_apex_int,
        apex_rt_shift=sg_apex_rt - raw_apex_rt,
        apex_intensity_ratio=sg_apex_int / raw_apex_int if raw_apex_int > 0 else 1.0,
        apex_intensity_change_pct=((sg_apex_int - raw_apex_int) / raw_apex_int * 100.0) if raw_apex_int > 0 else 0.0,
        peak_area=sg_area,
        peak_area_ratio=sg_area / raw_area if raw_area > 0 else 1.0,
        fwhm=sg_fwhm,
        fwhm_ratio=sg_fwhm / raw_fwhm if raw_fwhm > 0 else 1.0,
        fwhm_change_pct=((sg_fwhm - raw_fwhm) / raw_fwhm * 100.0) if raw_fwhm > 0 else 0.0,
        flank_noise_sigma=sg_noise,
        noise_reduction_pct=sg_noise_red,
        snr=sg_snr,
        snr_gain_factor=sg_snr / raw_snr if raw_snr > 0 else 1.0,
        valley_to_peak_ratio=sg_vpr,
        num_peaks_detected_sn2=sg_n_peaks,
    )

    native_id = chrom.getNativeID() if hasattr(chrom, "getNativeID") and chrom.getNativeID() else trace_id

    return RealTraceEvaluation(
        trace_id=trace_id,
        native_id=native_id,
        source_file=source_file,
        num_points=n_pts,
        raw_apex_rt=raw_apex_rt,
        raw_apex_intensity=raw_apex_int,
        raw_peak_area=raw_area,
        raw_fwhm=raw_fwhm,
        raw_flank_noise_sigma=raw_noise,
        raw_snr=raw_snr,
        raw_valley_to_peak_ratio=raw_vpr,
        raw_peaks_detected_sn2=raw_n_peaks,
        intensity_stratum=intensity_stratum,
        width_stratum=width_stratum,
        modified_sinc=ms_metrics,
        savitzky_golay=sg_metrics,
    )


def summarize_real_dataset_results(
    evaluations: List[RealTraceEvaluation],
    dataset_accession: str,
    dataset_name: str,
    ms_params: Dict[str, Any],
    sg_params: Dict[str, Any],
    ms_time_us: float = 0.0,
    sg_time_us: float = 0.0,
) -> RealDatasetBenchmarkResult:
    """
    Compute aggregate empirical statistics across all evaluated real chromatograms.
    Includes Mean, Median, and Interquartile Range (IQR = Q75 - Q25) for robust reporting.
    """
    if not evaluations:
        return RealDatasetBenchmarkResult(
            dataset_accession=dataset_accession,
            dataset_name=dataset_name,
            total_traces_evaluated=0,
            ms_params=ms_params,
            sg_params=sg_params,
            traces=[],
        )

    def calc_iqr(vals: List[float]) -> float:
        if not vals:
            return 0.0
        return float(np.percentile(vals, 75) - np.percentile(vals, 25))

    ms_shifts = [e.modified_sinc.apex_rt_shift for e in evaluations]
    sg_shifts = [e.savitzky_golay.apex_rt_shift for e in evaluations]

    ms_apex_ratios = [e.modified_sinc.apex_intensity_ratio for e in evaluations]
    sg_apex_ratios = [e.savitzky_golay.apex_intensity_ratio for e in evaluations]

    ms_area_ratios = [e.modified_sinc.peak_area_ratio for e in evaluations]
    sg_area_ratios = [e.savitzky_golay.peak_area_ratio for e in evaluations]

    ms_fwhm_ratios = [e.modified_sinc.fwhm_ratio for e in evaluations]
    sg_fwhm_ratios = [e.savitzky_golay.fwhm_ratio for e in evaluations]

    ms_noise_reds = [e.modified_sinc.noise_reduction_pct for e in evaluations]
    sg_noise_reds = [e.savitzky_golay.noise_reduction_pct for e in evaluations]

    raw_p_counts = [e.raw_peaks_detected_sn2 for e in evaluations]
    ms_p_counts = [e.modified_sinc.num_peaks_detected_sn2 for e in evaluations]
    sg_p_counts = [e.savitzky_golay.num_peaks_detected_sn2 for e in evaluations]

    exact_matches = sum(1 for ms_c, sg_c in zip(ms_p_counts, sg_p_counts) if ms_c == sg_c)
    exact_match_pct = (exact_matches / len(evaluations)) * 100.0

    return RealDatasetBenchmarkResult(
        dataset_accession=dataset_accession,
        dataset_name=dataset_name,
        total_traces_evaluated=len(evaluations),
        ms_params=ms_params,
        sg_params=sg_params,
        traces=evaluations,
        ms_mean_apex_shift=float(np.mean(ms_shifts)),
        ms_median_apex_shift=float(np.median(ms_shifts)),
        ms_iqr_apex_shift=calc_iqr(ms_shifts),
        sg_mean_apex_shift=float(np.mean(sg_shifts)),
        sg_median_apex_shift=float(np.median(sg_shifts)),
        sg_iqr_apex_shift=calc_iqr(sg_shifts),
        ms_mean_apex_ratio=float(np.mean(ms_apex_ratios)),
        ms_median_apex_ratio=float(np.median(ms_apex_ratios)),
        ms_iqr_apex_ratio=calc_iqr(ms_apex_ratios),
        sg_mean_apex_ratio=float(np.mean(sg_apex_ratios)),
        sg_median_apex_ratio=float(np.median(sg_apex_ratios)),
        sg_iqr_apex_ratio=calc_iqr(sg_apex_ratios),
        ms_mean_area_ratio=float(np.mean(ms_area_ratios)),
        ms_median_area_ratio=float(np.median(ms_area_ratios)),
        ms_iqr_area_ratio=calc_iqr(ms_area_ratios),
        sg_mean_area_ratio=float(np.mean(sg_area_ratios)),
        sg_median_area_ratio=float(np.median(sg_area_ratios)),
        sg_iqr_area_ratio=calc_iqr(sg_area_ratios),
        ms_mean_fwhm_ratio=float(np.mean(ms_fwhm_ratios)),
        ms_median_fwhm_ratio=float(np.median(ms_fwhm_ratios)),
        ms_iqr_fwhm_ratio=calc_iqr(ms_fwhm_ratios),
        sg_mean_fwhm_ratio=float(np.mean(sg_fwhm_ratios)),
        sg_median_fwhm_ratio=float(np.median(sg_fwhm_ratios)),
        sg_iqr_fwhm_ratio=calc_iqr(sg_fwhm_ratios),
        ms_mean_noise_red_pct=float(np.mean(ms_noise_reds)),
        ms_median_noise_red_pct=float(np.median(ms_noise_reds)),
        ms_iqr_noise_red_pct=calc_iqr(ms_noise_reds),
        sg_mean_noise_red_pct=float(np.mean(sg_noise_reds)),
        sg_median_noise_red_pct=float(np.median(sg_noise_reds)),
        sg_iqr_noise_red_pct=calc_iqr(sg_noise_reds),
        mean_raw_peaks_sn2=float(np.mean(raw_p_counts)),
        mean_ms_peaks_sn2=float(np.mean(ms_p_counts)),
        mean_sg_peaks_sn2=float(np.mean(sg_p_counts)),
        peak_count_exact_match_pct=exact_match_pct,
        ms_mean_runtime_us=ms_time_us,
        sg_mean_runtime_us=sg_time_us,
        speed_ratio_ms_sg=ms_time_us / sg_time_us if sg_time_us > 0 else 1.0,
    )



@dataclass
class QCReplicateConsistencyResult:
    """
    Empirical consistency and reproducibility across repeated injections (e.g. MTBLS404 pooled QCs).
    Computes Coefficient of Variation (CV% = std / mean * 100%) for apex RT, apex intensity, and area.
    """
    feature_id: str
    num_replicates: int
    raw_apex_rt_cv_pct: float
    ms_apex_rt_cv_pct: float
    sg_apex_rt_cv_pct: float
    raw_apex_int_cv_pct: float
    ms_apex_int_cv_pct: float
    sg_apex_int_cv_pct: float
    raw_area_cv_pct: float
    ms_area_cv_pct: float
    sg_area_cv_pct: float

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def evaluate_qc_replicates(
    replicate_chroms: List[pyopenms.MSChromatogram],
    feature_id: str = "qc_feature",
    ms_params: Optional[Dict[str, Any]] = None,
    sg_params: Optional[Dict[str, Any]] = None,
) -> QCReplicateConsistencyResult:
    """
    Evaluate empirical replicate consistency across repeated injections.
    Measures CV% of apex RT, intensity, and area for Raw vs Modified Sinc vs Savitzky-Golay.
    """
    if ms_params is None:
        ms_params = {"degree": 6, "m": 12, "is_ms1": False}
    if sg_params is None:
        sg_params = {"frame_length": 15, "polynomial_order": 4}

    ms_smoother = pyopenms.ModifiedSincSmoother()
    p_ms = ms_smoother.getParameters()
    p_ms.setValue("degree", ms_params.get("degree", 6))
    p_ms.setValue("m", ms_params.get("m", 12))
    p_ms.setValue("is_ms1", "true" if ms_params.get("is_ms1", False) else "false")
    ms_smoother.setParameters(p_ms)

    sg_smoother = pyopenms.SavitzkyGolayFilter()
    p_sg = sg_smoother.getParameters()
    p_sg.setValue("frame_length", sg_params.get("frame_length", 15))
    p_sg.setValue("polynomial_order", sg_params.get("polynomial_order", 4))
    sg_smoother.setParameters(p_sg)

    raw_rts, raw_ints, raw_areas = [], [], []
    ms_rts, ms_ints, ms_areas = [], [], []
    sg_rts, sg_ints, sg_areas = [], [], []

    for c in replicate_chroms:
        rt_raw, int_raw = c.get_peaks()
        rt_arr = np.array(rt_raw, dtype=np.float64)
        y_raw = np.array(int_raw, dtype=np.float64)

        if len(y_raw) < 5:
            continue

        raw_rt, raw_int = parabolic_apex_interpolation(rt_arr, y_raw)
        base_raw = float(np.percentile(y_raw, 10))
        raw_area = float(trapz_fn(y_raw - base_raw, rt_arr))
        raw_rts.append(raw_rt)
        raw_ints.append(raw_int)
        raw_areas.append(raw_area)

        c_ms = pyopenms.MSChromatogram(c)
        ms_smoother.filter(c_ms)
        _, int_ms = c_ms.get_peaks()
        y_ms = np.array(int_ms, dtype=np.float64)
        ms_rt, ms_int = parabolic_apex_interpolation(rt_arr, y_ms)
        base_ms = float(np.percentile(y_ms, 10))
        ms_area = float(trapz_fn(y_ms - base_ms, rt_arr))
        ms_rts.append(ms_rt)
        ms_ints.append(ms_int)
        ms_areas.append(ms_area)

        c_sg = pyopenms.MSChromatogram(c)
        sg_smoother.filter(c_sg)
        _, int_sg = c_sg.get_peaks()
        y_sg = np.array(int_sg, dtype=np.float64)
        sg_rt, sg_int = parabolic_apex_interpolation(rt_arr, y_sg)
        base_sg = float(np.percentile(y_sg, 10))
        sg_area = float(trapz_fn(y_sg - base_sg, rt_arr))
        sg_rts.append(sg_rt)
        sg_ints.append(sg_int)
        sg_areas.append(sg_area)

    def calc_cv(vals: List[float]) -> float:
        if len(vals) < 2:
            return 0.0
        m = float(np.mean(vals))
        s = float(np.std(vals, ddof=1))
        return (s / m * 100.0) if abs(m) > 1e-9 else 0.0

    return QCReplicateConsistencyResult(
        feature_id=feature_id,
        num_replicates=len(raw_rts),
        raw_apex_rt_cv_pct=calc_cv(raw_rts),
        ms_apex_rt_cv_pct=calc_cv(ms_rts),
        sg_apex_rt_cv_pct=calc_cv(sg_rts),
        raw_apex_int_cv_pct=calc_cv(raw_ints),
        ms_apex_int_cv_pct=calc_cv(ms_ints),
        sg_apex_int_cv_pct=calc_cv(sg_ints),
        raw_area_cv_pct=calc_cv(raw_areas),
        ms_area_cv_pct=calc_cv(ms_areas),
        sg_area_cv_pct=calc_cv(sg_areas),
    )

