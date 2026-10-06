"""
Synthetic chromatogram generation module for OpenMS smoothing benchmark.
GitHub Issue #10425: Benchmark ModifiedSincSmoother against Savitzky–Golay.

Generates deterministic chromatograms with known ground truth peak shapes:
- Narrow peaks (few points across FWHM) with sampling phase offsets (0.0, 0.25, 0.5, 0.75)
- Broad peaks (many points across FWHM)
- Overlapping peaks (doublets with known valley depths and total area)
- Low-intensity peaks (near limit of detection)
- Asymmetric tailing peaks (Exponentially Modified Gaussian, EMG)
- Realistic composite chromatograms combining multiple peak types
- Physical detector noise models (Gaussian and heteroscedastic) without artificial zero-clipping
- Explicit separation of tuning (training) and held-out evaluation datasets
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Dict, Any

import numpy as np
from scipy.special import erfc
import pyopenms

# Numerical trapezoid integration compatible with numpy 1.x and 2.x
trapz_fn = getattr(np, "trapezoid", getattr(np, "trapz", None))


@dataclass
class PeakMetadata:
    """Ground truth metadata for an individual peak."""
    peak_type: str
    true_apex_rt: float
    true_apex_intensity: float
    true_fwhm: float
    true_area: float
    rt_start: float
    rt_end: float


@dataclass
class SyntheticChromatogram:
    """Container holding ground-truth and noisy chromatograms."""
    name: str
    rt: np.ndarray
    intensity_true: np.ndarray
    intensity_noisy: np.ndarray
    noise_sigma: float
    noise_type: str
    baseline_level: float = 0.0
    phase_offset: float = 0.0
    peaks: List[PeakMetadata] = field(default_factory=list)
    baseline_mask: Optional[np.ndarray] = None  # True where signal is pure baseline
    is_overlapping_doublet: bool = False
    doublet_total_area: float = 0.0
    doublet_rt_start: float = 0.0
    doublet_rt_end: float = 0.0

    def to_openms_true(self) -> pyopenms.MSChromatogram:
        """Convert ground truth signal to pyopenms.MSChromatogram."""
        chrom = pyopenms.MSChromatogram()
        chrom.setName(f"{self.name}_true")
        chrom.set_peaks((self.rt.tolist(), self.intensity_true.tolist()))
        return chrom

    def to_openms_noisy(self) -> pyopenms.MSChromatogram:
        """Convert noisy signal to pyopenms.MSChromatogram."""
        chrom = pyopenms.MSChromatogram()
        chrom.setName(f"{self.name}_noisy")
        chrom.set_peaks((self.rt.tolist(), self.intensity_noisy.tolist()))
        return chrom


def generate_gaussian_peak(
    rt: np.ndarray,
    center_rt: float,
    height: float,
    fwhm: float,
) -> Tuple[np.ndarray, PeakMetadata]:
    """Generate a single Gaussian peak on the RT axis."""
    sigma = fwhm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    intensity = height * np.exp(-0.5 * ((rt - center_rt) / sigma) ** 2)
    true_area = height * sigma * math.sqrt(2.0 * math.pi)

    # Boundaries where peak drops to ~0.3% of height (~ 3.5 sigma)
    rt_start = center_rt - 3.5 * sigma
    rt_end = center_rt + 3.5 * sigma

    meta = PeakMetadata(
        peak_type="gaussian",
        true_apex_rt=center_rt,
        true_apex_intensity=height,
        true_fwhm=fwhm,
        true_area=true_area,
        rt_start=rt_start,
        rt_end=rt_end,
    )
    return intensity, meta


def generate_emg_peak(
    rt: np.ndarray,
    center_rt: float,
    height: float,
    fwhm: float,
    tau: float,
) -> Tuple[np.ndarray, PeakMetadata]:
    """
    Generate an Exponentially Modified Gaussian (EMG) peak.
    Models chromatographic peak tailing analytically.
    """
    sigma = fwhm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    z = (sigma / tau - (rt - center_rt) / sigma) / math.sqrt(2.0)
    exp_arg = 0.5 * (sigma / tau) ** 2 - (rt - center_rt) / tau
    exp_arg = np.clip(exp_arg, -500.0, 500.0)
    profile = (sigma / tau) * math.sqrt(math.pi / 2.0) * np.exp(exp_arg) * erfc(z)

    max_val = np.max(profile)
    if max_val > 0:
        profile = profile * (height / max_val)

    apex_idx = int(np.argmax(profile))
    actual_apex_rt = float(rt[apex_idx])
    actual_apex_height = float(profile[apex_idx])

    # Interpolate half-maximum crossings for true FWHM
    half_max = actual_apex_height / 2.0
    above_half = np.where(profile >= half_max)[0]
    if len(above_half) > 1:
        actual_fwhm = float(rt[above_half[-1]] - rt[above_half[0]])
    else:
        actual_fwhm = fwhm

    actual_area = float(trapz_fn(profile, rt))

    rt_start = center_rt - 3.5 * sigma
    rt_end = center_rt + 3.5 * sigma + 4.5 * tau

    meta = PeakMetadata(
        peak_type="emg_tailing",
        true_apex_rt=actual_apex_rt,
        true_apex_intensity=actual_apex_height,
        true_fwhm=actual_fwhm,
        true_area=actual_area,
        rt_start=rt_start,
        rt_end=rt_end,
    )
    return profile, meta


def add_noise(
    signal: np.ndarray,
    noise_sigma: float,
    noise_type: str = "gaussian",
    seed: int = 42,
) -> np.ndarray:
    """
    Add zero-mean deterministic noise without artificial zero-clipping.
    Physical non-negativity is achieved naturally via realistic detector baseline offset,
    avoiding half-normal distribution distortion.
    """
    rng = np.random.default_rng(seed)
    if noise_type == "gaussian":
        noise = rng.normal(0.0, noise_sigma, size=len(signal))
        return signal + noise
    elif noise_type == "heteroscedastic":
        # Variance proportional to signal intensity + baseline electronic noise
        local_sigma = np.sqrt(noise_sigma ** 2 + 0.02 * np.maximum(0.0, signal))
        noise = rng.normal(0.0, 1.0, size=len(signal)) * local_sigma
        return signal + noise
    else:
        raise ValueError(f"Unknown noise_type: {noise_type}")


class SyntheticChromatogramGenerator:
    """
    Generator for reproducible benchmark chromatograms.
    Supports isolated narrow/broad peaks, overlapping doublets/triplets,
    low-intensity peaks, sampling phase shifts, and held-out test datasets.
    """

    def __init__(self, seed: int = 42):
        self.seed = seed

    def create_narrow_peak(
        self,
        dt: float = 1.0,
        n_points: int = 150,
        height: float = 1000.0,
        fwhm: float = 3.0,  # ~3 points across FWHM
        phase_offset: float = 0.0,  # offset as fraction of dt
        noise_level: float = 0.03,
        noise_type: str = "gaussian",
        baseline_level: float = 150.0,
    ) -> SyntheticChromatogram:
        """Create a narrow chromatographic peak with optional sub-scan sampling phase offset."""
        rt = np.arange(n_points, dtype=np.float64) * dt
        center_rt = (n_points // 2 + phase_offset) * dt
        peak_intensity, meta = generate_gaussian_peak(rt, center_rt, height, fwhm)

        intensity_true = baseline_level + peak_intensity
        sigma_noise = noise_level * height
        intensity_noisy = add_noise(intensity_true, sigma_noise, noise_type, seed=self.seed)

        baseline_mask = (rt < meta.rt_start) | (rt > meta.rt_end)

        return SyntheticChromatogram(
            name=f"narrow_peak_phase_{phase_offset:.2f}",
            rt=rt,
            intensity_true=intensity_true,
            intensity_noisy=intensity_noisy,
            noise_sigma=sigma_noise,
            noise_type=noise_type,
            baseline_level=baseline_level,
            phase_offset=phase_offset,
            peaks=[meta],
            baseline_mask=baseline_mask,
        )

    def create_broad_peak(
        self,
        dt: float = 1.0,
        n_points: int = 200,
        height: float = 1000.0,
        fwhm: float = 30.0,  # ~30 points across FWHM
        noise_level: float = 0.03,
        noise_type: str = "gaussian",
        baseline_level: float = 150.0,
    ) -> SyntheticChromatogram:
        """Create a well-sampled broad chromatographic peak."""
        rt = np.arange(n_points, dtype=np.float64) * dt
        center_rt = float(rt[n_points // 2])
        peak_intensity, meta = generate_gaussian_peak(rt, center_rt, height, fwhm)

        intensity_true = baseline_level + peak_intensity
        sigma_noise = noise_level * height
        intensity_noisy = add_noise(intensity_true, sigma_noise, noise_type, seed=self.seed + 1)
        baseline_mask = (rt < meta.rt_start) | (rt > meta.rt_end)

        return SyntheticChromatogram(
            name="broad_peak",
            rt=rt,
            intensity_true=intensity_true,
            intensity_noisy=intensity_noisy,
            noise_sigma=sigma_noise,
            noise_type=noise_type,
            baseline_level=baseline_level,
            peaks=[meta],
            baseline_mask=baseline_mask,
        )

    def create_overlapping_peaks(
        self,
        dt: float = 1.0,
        n_points: int = 250,
        height1: float = 1000.0,
        height2: float = 800.0,
        fwhm: float = 12.0,
        separation_factor: float = 1.3,
        noise_level: float = 0.02,
        noise_type: str = "gaussian",
        baseline_level: float = 150.0,
    ) -> SyntheticChromatogram:
        """
        Create overlapping doublet peaks with ground truth total doublet area
        to ensure scientifically valid area evaluation.
        """
        rt = np.arange(n_points, dtype=np.float64) * dt
        center1 = float(rt[n_points // 2]) - (separation_factor * fwhm) / 2.0
        center2 = center1 + separation_factor * fwhm

        i1, m1 = generate_gaussian_peak(rt, center1, height1, fwhm)
        i2, m2 = generate_gaussian_peak(rt, center2, height2, fwhm)
        peak_intensity = i1 + i2

        intensity_true = baseline_level + peak_intensity
        sigma_noise = noise_level * height1
        intensity_noisy = add_noise(intensity_true, sigma_noise, noise_type, seed=self.seed + 2)

        doublet_start = min(m1.rt_start, m2.rt_start)
        doublet_end = max(m1.rt_end, m2.rt_end)
        midpoint = (center1 + center2) / 2.0
        m1.rt_end = midpoint
        m2.rt_start = midpoint
        baseline_mask = (rt < doublet_start) | (rt > doublet_end)

        return SyntheticChromatogram(
            name="overlapping_peaks",
            rt=rt,
            intensity_true=intensity_true,
            intensity_noisy=intensity_noisy,
            noise_sigma=sigma_noise,
            noise_type=noise_type,
            baseline_level=baseline_level,
            peaks=[m1, m2],
            baseline_mask=baseline_mask,
            is_overlapping_doublet=True,
            doublet_total_area=m1.true_area + m2.true_area,
            doublet_rt_start=doublet_start,
            doublet_rt_end=doublet_end,
        )

    def create_low_intensity_peak(
        self,
        dt: float = 1.0,
        n_points: int = 150,
        height: float = 60.0,  # small peak near LOD
        fwhm: float = 10.0,
        noise_sigma: float = 15.0,  # SNR ~ 4
        noise_type: str = "gaussian",
        baseline_level: float = 80.0,
    ) -> SyntheticChromatogram:
        """Create a low-intensity peak near the limit of detection (LOD)."""
        rt = np.arange(n_points, dtype=np.float64) * dt
        center_rt = float(rt[n_points // 2])
        peak_intensity, meta = generate_gaussian_peak(rt, center_rt, height, fwhm)

        intensity_true = baseline_level + peak_intensity
        intensity_noisy = add_noise(intensity_true, noise_sigma, noise_type, seed=self.seed + 3)
        baseline_mask = (rt < meta.rt_start) | (rt > meta.rt_end)

        return SyntheticChromatogram(
            name="low_intensity_peak",
            rt=rt,
            intensity_true=intensity_true,
            intensity_noisy=intensity_noisy,
            noise_sigma=noise_sigma,
            noise_type=noise_type,
            baseline_level=baseline_level,
            peaks=[meta],
            baseline_mask=baseline_mask,
        )

    def create_tailing_peak(
        self,
        dt: float = 1.0,
        n_points: int = 200,
        height: float = 1000.0,
        fwhm: float = 10.0,
        tau: float = 8.0,
        noise_level: float = 0.02,
        noise_type: str = "gaussian",
        baseline_level: float = 150.0,
    ) -> SyntheticChromatogram:
        """Create an asymmetric tailing peak (Exponentially Modified Gaussian)."""
        rt = np.arange(n_points, dtype=np.float64) * dt
        center_rt = float(rt[n_points // 3])
        peak_intensity, meta = generate_emg_peak(rt, center_rt, height, fwhm, tau)

        intensity_true = baseline_level + peak_intensity
        sigma_noise = noise_level * height
        intensity_noisy = add_noise(intensity_true, sigma_noise, noise_type, seed=self.seed + 4)
        baseline_mask = (rt < meta.rt_start) | (rt > meta.rt_end)

        return SyntheticChromatogram(
            name="tailing_peak",
            rt=rt,
            intensity_true=intensity_true,
            intensity_noisy=intensity_noisy,
            noise_sigma=sigma_noise,
            noise_type=noise_type,
            baseline_level=baseline_level,
            peaks=[meta],
            baseline_mask=baseline_mask,
        )

    def create_composite_chromatogram(
        self,
        dt: float = 0.5,
        total_time: float = 300.0,
        noise_level: float = 0.02,
        noise_type: str = "heteroscedastic",
        baseline_level: float = 100.0,
        peak_configs: Optional[Dict[str, Dict[str, float]]] = None,
    ) -> SyntheticChromatogram:
        """
        Create a comprehensive realistic chromatogram containing:
        - 1 narrow peak
        - 1 overlapping doublet
        - 1 asymmetric tailing peak
        - 1 low-intensity peak
        - 1 broad high-intensity peak
        """
        rt = np.arange(0.0, total_time, dt, dtype=np.float64)
        peak_intensity = np.zeros_like(rt)
        peaks: List[PeakMetadata] = []

        configs = {
            "narrow": {"center_rt": 40.0, "height": 2500.0, "fwhm": 3.0},
            "d1": {"center_rt": 95.0, "height": 3000.0, "fwhm": 10.0},
            "d2": {"center_rt": 107.0, "height": 1800.0, "fwhm": 10.0},
            "emg": {"center_rt": 160.0, "height": 2000.0, "fwhm": 8.0, "tau": 6.0},
            "low": {"center_rt": 215.0, "height": 140.0, "fwhm": 7.0},
            "broad": {"center_rt": 260.0, "height": 4500.0, "fwhm": 25.0},
        }
        if peak_configs:
            for k, v in peak_configs.items():
                if k in configs:
                    configs[k].update(v)

        # 1. Narrow peak
        c_n = configs["narrow"]
        i_narrow, m_narrow = generate_gaussian_peak(
            rt, center_rt=c_n["center_rt"], height=c_n["height"], fwhm=c_n["fwhm"]
        )
        peak_intensity += i_narrow
        peaks.append(m_narrow)

        # 2. Overlapping doublet
        c_d1 = configs["d1"]
        c_d2 = configs["d2"]
        i_d1, m_d1 = generate_gaussian_peak(
            rt, center_rt=c_d1["center_rt"], height=c_d1["height"], fwhm=c_d1["fwhm"]
        )
        i_d2, m_d2 = generate_gaussian_peak(
            rt, center_rt=c_d2["center_rt"], height=c_d2["height"], fwhm=c_d2["fwhm"]
        )
        midpoint_d = (c_d1["center_rt"] + c_d2["center_rt"]) / 2.0
        m_d1.rt_end = midpoint_d
        m_d2.rt_start = midpoint_d
        doublet_rt_start = min(m_d1.rt_start, m_d2.rt_start)
        doublet_rt_end = max(m_d1.rt_end, m_d2.rt_end)
        peak_intensity += i_d1 + i_d2
        peaks.extend([m_d1, m_d2])

        # 3. Asymmetric tailing peak (EMG)
        c_emg = configs["emg"]
        i_emg, m_emg = generate_emg_peak(
            rt,
            center_rt=c_emg["center_rt"],
            height=c_emg["height"],
            fwhm=c_emg["fwhm"],
            tau=c_emg["tau"],
        )
        peak_intensity += i_emg
        peaks.append(m_emg)

        # 4. Low-intensity peak near LOD
        c_low = configs["low"]
        i_low, m_low = generate_gaussian_peak(
            rt, center_rt=c_low["center_rt"], height=c_low["height"], fwhm=c_low["fwhm"]
        )
        peak_intensity += i_low
        peaks.append(m_low)

        # 5. Broad high-intensity peak
        c_broad = configs["broad"]
        i_broad, m_broad = generate_gaussian_peak(
            rt, center_rt=c_broad["center_rt"], height=c_broad["height"], fwhm=c_broad["fwhm"]
        )
        peak_intensity += i_broad
        peaks.append(m_broad)

        intensity_true = baseline_level + peak_intensity

        baseline_mask = np.ones(len(rt), dtype=bool)
        for p in peaks:
            peak_region = (rt >= p.rt_start) & (rt <= p.rt_end)
            baseline_mask[peak_region] = False

        sigma_noise = noise_level * 3000.0
        intensity_noisy = add_noise(intensity_true, sigma_noise, noise_type, seed=self.seed + 5)

        return SyntheticChromatogram(
            name="composite_chromatogram",
            rt=rt,
            intensity_true=intensity_true,
            intensity_noisy=intensity_noisy,
            noise_sigma=sigma_noise,
            noise_type=noise_type,
            baseline_level=baseline_level,
            peaks=peaks,
            baseline_mask=baseline_mask,
            is_overlapping_doublet=True,
            doublet_total_area=m_d1.true_area + m_d2.true_area,
            doublet_rt_start=doublet_rt_start,
            doublet_rt_end=doublet_rt_end,
        )

    def generate_tuning_datasets(self) -> Dict[str, SyntheticChromatogram]:
        """
        Generate calibration/tuning datasets.
        Uses a separate deterministic seed (self.seed + 10000) and distinct peak locations and shapes
        to ensure parameter tuning is completely isolated from test data.
        """
        tuning_gen = SyntheticChromatogramGenerator(seed=self.seed + 10000)
        tuning_peaks = {
            "narrow": {"center_rt": 45.0, "height": 2400.0, "fwhm": 3.4},
            "d1": {"center_rt": 100.0, "height": 2800.0, "fwhm": 11.0},
            "d2": {"center_rt": 114.0, "height": 1900.0, "fwhm": 11.0},
            "emg": {"center_rt": 172.0, "height": 2100.0, "fwhm": 9.0, "tau": 5.5},
            "low": {"center_rt": 222.0, "height": 150.0, "fwhm": 7.5},
            "broad": {"center_rt": 268.0, "height": 4200.0, "fwhm": 24.0},
        }
        return {
            "tuning_narrow": tuning_gen.create_narrow_peak(n_points=160, fwhm=3.5, height=1200.0),
            "tuning_broad": tuning_gen.create_broad_peak(n_points=220, fwhm=28.0, height=1200.0),
            "tuning_composite": tuning_gen.create_composite_chromatogram(
                dt=0.5,
                total_time=320.0,
                peak_configs=tuning_peaks,
            ),
        }

    def generate_test_datasets(self) -> Dict[str, SyntheticChromatogram]:
        """Generate held-out evaluation test datasets."""
        return {
            "narrow": self.create_narrow_peak(phase_offset=0.0),
            "broad": self.create_broad_peak(),
            "overlapping": self.create_overlapping_peaks(),
            "low_intensity": self.create_low_intensity_peak(),
            "tailing": self.create_tailing_peak(),
            "composite": self.create_composite_chromatogram(),
        }

    def generate_phase_sweep_narrow_peaks(self) -> Dict[str, SyntheticChromatogram]:
        """Generate narrow peaks across 4 sub-scan sampling phases: 0.0, 0.25, 0.50, 0.75."""
        phases = [0.0, 0.25, 0.50, 0.75]
        res = {}
        for p in phases:
            res[f"phase_{p:.2f}"] = self.create_narrow_peak(phase_offset=p, dt=1.0, fwhm=3.0)
        return res
