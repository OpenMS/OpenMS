"""
Plotting and visualization module for OpenMS smoothing benchmark.
GitHub Issue #10425: Benchmark ModifiedSincSmoother against Savitzky–Golay.

Generates publication-quality figures:
1. chromatogram_profiles_comparison.png: Side-by-side profile overlays for different peak types
2. sampling_phase_sensitivity.png: Apex and width distortion across sub-scan sampling phases
3. sn_threshold_peak_picking.png: Downstream PeakPickerHiRes F1-score across S/N thresholds
4. overlapping_peak_resolution.png: Valley preservation for co-eluting doublet peaks
5. runtime_scaling_vs_length.png: Microbenchmark throughput (ns/pt) and time per call
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import List, Dict, Any, Optional

import matplotlib
matplotlib.use("Agg")  # Headless rendering
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.titlesize": 14,
    "lines.linewidth": 1.7,
    "axes.grid": True,
    "grid.alpha": 0.4,
    "grid.linestyle": ":",
})


def plot_chromatogram_profiles(
    profiles_dict: Dict[str, Dict[str, Any]],
    output_path: str,
):
    """
    Generate 4-panel comparison figure showing Ground Truth, Raw Noisy,
    Modified Sinc, and Savitzky-Golay for:
    - Narrow Peak
    - Broad Peak
    - Overlapping Doublet
    - Low-Intensity Peak
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()

    scenarios = [
        ("narrow", "Narrow Peak (~3 pts FWHM, UHPLC-like)"),
        ("broad", "Broad Peak (~30 pts FWHM, Well Sampled)"),
        ("overlapping", "Overlapping Doublet (Co-eluting Peaks)"),
        ("low_intensity", "Low-Intensity Peak (Near LOD, SNR ~ 4)"),
    ]

    for idx, (key, title) in enumerate(scenarios):
        ax = axes[idx]
        if key not in profiles_dict:
            ax.text(0.5, 0.5, f"Data {key} not found", ha="center", va="center")
            continue

        p = profiles_dict[key]
        rt = p["rt"]
        y_true = p["y_true"]
        y_noisy = p["y_noisy"]
        y_ms = p["y_ms"]
        y_sg = p["y_sg"]

        ax.plot(rt, y_noisy, color="#94a3b8", alpha=0.6, label="Raw Noisy Trace", lw=1.2)
        ax.plot(rt, y_true, color="#0f172a", linestyle="--", label="Ground Truth", lw=2.0)
        ax.plot(rt, y_ms, color="#2563eb", label="Modified Sinc (Tuned)", lw=1.8)
        ax.plot(rt, y_sg, color="#dc2626", label="Savitzky-Golay (Tuned)", lw=1.8)

        ax.set_title(title, fontweight="bold")
        ax.set_xlabel("Retention Time (s)")
        ax.set_ylabel("Intensity")
        ax.legend(loc="upper right", framealpha=0.9)

    plt.tight_layout()
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300)
    plt.close()


def plot_sampling_phase_sensitivity(
    phase_results: Dict[str, Any],
    output_path: str,
):
    """
    Plot apex attenuation error % and FWHM error % across sub-scan sampling phases (0.0 to 0.75).
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    phases = [phase_results[k]["phase"] for k in phase_results]
    ms_apex = [phase_results[k]["ms_apex_err"] for k in phase_results]
    sg_apex = [phase_results[k]["sg_apex_err"] for k in phase_results]

    ms_fwhm = [phase_results[k]["ms_fwhm_err"] for k in phase_results]
    sg_fwhm = [phase_results[k]["sg_fwhm_err"] for k in phase_results]

    # Apex Intensity Error
    ax1.plot(phases, ms_apex, "o-", color="#2563eb", label="Modified Sinc", lw=2)
    ax1.plot(phases, sg_apex, "s-", color="#dc2626", label="Savitzky-Golay", lw=2)
    ax1.set_title("Narrow Peak Apex Error vs Sampling Phase", fontweight="bold")
    ax1.set_xlabel("Sub-scan Sampling Phase Offset (Fraction of dt)")
    ax1.set_ylabel("Apex Intensity Error (%) [Closer to 0% is Better]")
    ax1.set_xticks(phases)
    ax1.legend()

    # FWHM Error
    ax2.plot(phases, ms_fwhm, "o-", color="#2563eb", label="Modified Sinc", lw=2)
    ax2.plot(phases, sg_fwhm, "s-", color="#dc2626", label="Savitzky-Golay", lw=2)
    ax2.set_title("Narrow Peak FWHM Distortion vs Sampling Phase", fontweight="bold")
    ax2.set_xlabel("Sub-scan Sampling Phase Offset (Fraction of dt)")
    ax2.set_ylabel("FWHM Distortion (%) [Closer to 0% is Better]")
    ax2.set_xticks(phases)
    ax2.legend()

    plt.tight_layout()
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300)
    plt.close()


def plot_sn_threshold_peak_picking(
    ms_metrics: Any,
    sg_metrics: Any,
    output_path: str,
):
    """
    Plot downstream peak picking F1 score across S/N thresholds (1.0, 2.0, 3.0, 5.0).
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    thresholds = sorted(ms_metrics.peak_picking_sensitivity.keys())
    ms_f1 = [ms_metrics.peak_picking_sensitivity[t].f1_score for t in thresholds]
    sg_f1 = [sg_metrics.peak_picking_sensitivity[t].f1_score for t in thresholds]

    x = np.arange(len(thresholds))
    width = 0.35

    ax.bar(x - width / 2, ms_f1, width, label="Modified Sinc", color="#2563eb", alpha=0.85)
    ax.bar(x + width / 2, sg_f1, width, label="Savitzky-Golay", color="#dc2626", alpha=0.85)

    ax.set_title("PeakPickerHiRes F1-Score across S/N Thresholds", fontweight="bold")
    ax.set_xlabel("PeakPickerHiRes S/N Threshold")
    ax.set_ylabel("F1 Score")
    ax.set_xticks(x)
    ax.set_xticklabels([f"S/N = {t}" for t in thresholds])
    ax.set_ylim(0.0, 1.05)
    ax.legend()

    plt.tight_layout()
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300)
    plt.close()


def plot_overlapping_resolution(
    valley_data: Dict[str, Any],
    output_path: str,
):
    """
    Bar chart showing valley-to-peak ratio for overlapping doublet.
    Lower VPR = sharper valley between peaks.
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    labels = ["Ground Truth", "Raw Noisy", "Modified Sinc", "Savitzky-Golay"]
    keys = ["true", "noisy", "ms", "sg"]
    values = [
        float(valley_data[k]) if valley_data.get(k) is not None else float("nan")
        for k in keys
    ]
    colors = ["#0f172a", "#94a3b8", "#2563eb", "#dc2626"]

    bars = ax.bar(labels, values, color=colors, width=0.55, edgecolor="black", linewidth=0.6)
    ax.set_ylabel("Valley-to-Peak Ratio (VPR)")
    ax.set_title("Overlapping Doublet Resolvability (Lower VPR = Sharper Valley)", fontweight="bold")

    for bar, val in zip(bars, values):
        missing = np.isnan(val)
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            0.01 if missing else val + 0.01,
            "N/A" if missing else f"{val:.3f}",
            ha="center",
            va="bottom",
            fontweight="bold",
        )

    plt.tight_layout()
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300)
    plt.close()


def plot_runtime_scaling(
    scaling_data: Dict[str, Any],
    output_path: str,
):
    """
    Plot computational runtime scaling O(N) vs chromatogram point count N.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    lengths = scaling_data["lengths"]
    ms_us = scaling_data["ms_per_call_us"]
    sg_us = scaling_data["sg_per_call_us"]
    ms_ns = scaling_data["ms_per_point_ns"]
    sg_ns = scaling_data["sg_per_point_ns"]

    # 1. Total runtime per call (log-log)
    ax1.plot(lengths, ms_us, "o-", color="#2563eb", lw=2, label="Modified Sinc Smoother")
    ax1.plot(lengths, sg_us, "s-", color="#dc2626", lw=2, label="Savitzky-Golay Filter")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel("Chromatogram Points (N)")
    ax1.set_ylabel("Execution Time per Call (microseconds)")
    ax1.set_title("Runtime Scaling vs Length N (Log-Log)", fontweight="bold")
    ax1.legend()

    # 2. Per-point throughput (ns/point)
    ax2.plot(lengths, ms_ns, "o-", color="#2563eb", lw=2, label="Modified Sinc")
    ax2.plot(lengths, sg_ns, "s-", color="#dc2626", lw=2, label="Savitzky-Golay")
    ax2.set_xscale("log")
    ax2.set_xlabel("Chromatogram Points (N)")
    ax2.set_ylabel("Time per Point (nanoseconds)")
    ax2.set_title("Throughput per Point (ns/point)", fontweight="bold")
    ax2.legend()

    plt.tight_layout()
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300)
    plt.close()


def plot_real_chromatogram_overlays(
    overlay_traces: List[Dict[str, Any]],
    output_path: str,
):
    """
    Generate multi-panel figure showing raw real traces vs Modified Sinc vs Savitzky-Golay.
    """
    n_plots = min(4, len(overlay_traces))
    if n_plots == 0:
        return

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()

    for idx in range(4):
        ax = axes[idx]
        if idx >= len(overlay_traces):
            ax.set_visible(False)
            continue

        item = overlay_traces[idx]
        rt = item["rt"]
        y_raw = item["raw"]
        y_ms = item["ms"]
        y_sg = item["sg"]
        title = item.get("title", f"Trace {idx+1}")

        ax.plot(rt, y_raw, color="#94a3b8", alpha=0.7, label="Raw LC-MS Trace", lw=1.2)
        ax.plot(rt, y_ms, color="#2563eb", label="Modified Sinc", lw=1.8)
        ax.plot(rt, y_sg, color="#dc2626", label="Savitzky-Golay", lw=1.8)

        ax.set_title(title, fontweight="bold")
        ax.set_xlabel("Retention Time (s)")
        ax.set_ylabel("Intensity")
        ax.legend(loc="upper right", framealpha=0.9)

    plt.tight_layout()
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300)
    plt.close()


def plot_real_data_distributions(
    evaluations: List[Any],
    output_path: str,
):
    """
    Plot comparative boxplots/distributions of empirical metrics on real chromatograms:
    1. Apex RT Shift (s)
    2. Apex Intensity Ratio
    3. Peak Area Ratio
    4. Baseline Noise Reduction (%)
    """
    if not evaluations:
        return

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    ms_shift = [e.modified_sinc.apex_rt_shift for e in evaluations]
    sg_shift = [e.savitzky_golay.apex_rt_shift for e in evaluations]

    ms_apex_r = [e.modified_sinc.apex_intensity_ratio for e in evaluations]
    sg_apex_r = [e.savitzky_golay.apex_intensity_ratio for e in evaluations]

    ms_area_r = [e.modified_sinc.peak_area_ratio for e in evaluations]
    sg_area_r = [e.savitzky_golay.peak_area_ratio for e in evaluations]

    ms_noise_red = [e.modified_sinc.noise_reduction_pct for e in evaluations]
    sg_noise_red = [e.savitzky_golay.noise_reduction_pct for e in evaluations]

    labels = ["Modified Sinc", "Savitzky-Golay"]

    # 1. Apex RT Shift
    axes[0, 0].boxplot([ms_shift, sg_shift], tick_labels=labels, patch_artist=True)
    axes[0, 0].axhline(0.0, color="gray", linestyle="--", alpha=0.6)
    axes[0, 0].set_title("Apex Retention Time Shift (s)", fontweight="bold")
    axes[0, 0].set_ylabel("Delta RT = Smooth - Raw (s)")

    # 2. Apex Intensity Ratio
    axes[0, 1].boxplot([ms_apex_r, sg_apex_r], tick_labels=labels, patch_artist=True)
    axes[0, 1].axhline(1.0, color="gray", linestyle="--", alpha=0.6)
    axes[0, 1].set_title("Apex Intensity Ratio (I_smooth / I_raw)", fontweight="bold")
    axes[0, 1].set_ylabel("Ratio")

    # 3. Peak Area Ratio
    axes[1, 0].boxplot([ms_area_r, sg_area_r], tick_labels=labels, patch_artist=True)
    axes[1, 0].axhline(1.0, color="gray", linestyle="--", alpha=0.6)
    axes[1, 0].set_title("Peak Area Ratio (Area_smooth / Area_raw)", fontweight="bold")
    axes[1, 0].set_ylabel("Ratio")

    # 4. Flank Baseline Noise Reduction
    axes[1, 1].boxplot([ms_noise_red, sg_noise_red], tick_labels=labels, patch_artist=True)
    axes[1, 1].set_title("Baseline Flank Noise Reduction (%)", fontweight="bold")
    axes[1, 1].set_ylabel("Noise Reduction (%) [Higher is Better]")

    plt.tight_layout()
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300)
    plt.close()


def plot_real_noise_vs_distortion(
    evaluations: List[Any],
    output_path: str,
):
    """
    Scatter plot of Flank Baseline Noise Reduction (%) vs Apex Attenuation (%) on real chromatograms.
    """
    if not evaluations:
        return

    fig, ax = plt.subplots(figsize=(8, 6))

    ms_x = [e.modified_sinc.noise_reduction_pct for e in evaluations]
    ms_y = [abs(e.modified_sinc.apex_intensity_change_pct) for e in evaluations]

    sg_x = [e.savitzky_golay.noise_reduction_pct for e in evaluations]
    sg_y = [abs(e.savitzky_golay.apex_intensity_change_pct) for e in evaluations]

    ax.scatter(ms_x, ms_y, color="#2563eb", alpha=0.7, label="Modified Sinc", edgecolors="none", s=40)
    ax.scatter(sg_x, sg_y, color="#dc2626", alpha=0.7, label="Savitzky-Golay", edgecolors="none", s=40)

    ax.set_title("Real Chromatograms: Noise Reduction vs Apex Attenuation", fontweight="bold")
    ax.set_xlabel("Flank Baseline Noise Reduction (%) [Higher is Better]")
    ax.set_ylabel("Apex Attenuation |Delta I| (%) [Lower is Better]")
    ax.legend()

    plt.tight_layout()
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300)
    plt.close()
