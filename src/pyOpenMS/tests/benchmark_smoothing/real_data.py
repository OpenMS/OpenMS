"""
Real-data chromatogram extraction and selection module for OpenMS smoothing benchmark.
GitHub Issue #10425: Benchmark ModifiedSincSmoother against Savitzky–Golay on chromatograms.

Provides:
- Metadata specifications for public reference datasets:
  * PASS00779 (PeptideAtlas / OpenSWATH tutorial, AB SCIEX TripleTOF 5600, SWATH-MS / DIA)
  * MTBLS404 (MetaboLights / Sacurine, Thermo LTQ-Orbitrap Discovery, LC-HRMS metabolomics)
  * Local OpenMS reference test chromatograms (packaged fixtures in src/tests/topp/)
- Deterministic chromatogram extraction from spectral mzML and chromatogram mzML
- Non-cherry-picking stratified selection based strictly on RAW signal properties
- Manifest generation recording seeds, selection criteria, rejection counts, and trace metadata
"""

from __future__ import annotations

import gzip
import json
import math
import os
import shutil
import tempfile
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple

import numpy as np
import pyopenms


@dataclass
class DatasetInfo:
    """Metadata specification for a real LC-MS benchmark dataset."""
    accession: str
    name: str
    repository: str
    official_url: str
    instrument: str
    acquisition_mode: str
    ionization_mode: str
    organism_or_matrix: str
    description: str
    recommended_files: List[str]
    is_public: bool = True
    unverified_notes: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


PASS00779_INFO = DatasetInfo(
    accession="PASS00779",
    name="M. tuberculosis OpenSWATH Reference Dataset",
    repository="PeptideAtlas / PASSEL",
    official_url="http://www.peptideatlas.org/PASS/PASS00779",
    instrument="AB SCIEX TripleTOF 5600",
    acquisition_mode="SWATH-MS / DIA (32 x 25 Da precursor isolation windows)",
    ionization_mode="Positive ESI",
    organism_or_matrix="Mycobacterium tuberculosis (Wayne dormancy model)",
    description=(
        "Standard OpenSWATH reference dataset containing 3 technical replicate SWATH-MS runs "
        "(R1, R2, R3) with associated assay library (Mtb_TubercuList-R27_iRT_UPS) and SWATH window files."
    ),
    recommended_files=[
        "olgas_K121026_001_SW_Wayne_R1_d00.mzML.gz",
        "olgas_K121026_007_SW_Wayne_R2_d00.mzML.gz",
        "olgas_K121026_013_SW_Wayne_R3_d00.mzML.gz",
        "SWATHwindows_analysis.tsv",
        "iRTassays.TraML",
    ],
    unverified_notes=[
        "Full mzML.gz files are ~3.3 GB each compressed (>10 GB uncompressed); full download is not suitable for routine CI.",
        "Raw vendor files require SCIEX Wiff reader or conversion via ProteoWizard msconvert.",
    ],
)

MTBLS404_INFO = DatasetInfo(
    accession="MTBLS404",
    name="Sacurine Adult Urinary Metabolome Dataset",
    repository="MetaboLights (EBI)",
    official_url="https://www.ebi.ac.uk/metabolights/MTBLS404",
    instrument="Thermo Fisher Scientific LTQ-Orbitrap Discovery",
    acquisition_mode="LC-HRMS Full-Scan MS1 (m/z 50–1000)",
    ionization_mode="Negative ESI",
    organism_or_matrix="Human urine (cohort of 184 adult volunteers + repeated pooled-QC injections)",
    description=(
        "Comprehensive clinical metabolomics benchmark containing 234 LC-HRMS negative mode files, "
        "including repeated pooled quality control (QC) injections throughout the analytical sequence."
    ),
    recommended_files=[
        "QC01.mzML",
        "QC02.mzML",
        "QC03.mzML",
        "QC04.mzML",
    ],
    unverified_notes=[
        "Dataset contains 234 total injections; only a small subset of 3-6 pooled-QC files should be downloaded.",
        "Centroided mzML files are provided alongside raw vendor .RAW profile files.",
    ],
)

LOCAL_FIXTURES_INFO = DatasetInfo(
    accession="LOCAL_FIXTURES",
    name="OpenMS In-Tree Representative Chromatograms",
    repository="OpenMS Repository (src/tests/topp/)",
    official_url="https://github.com/OpenMS/OpenMS",
    instrument="Various (SCIEX TripleTOF / QTRAP / Thermo)",
    acquisition_mode="MRM / SWATH / DDA chromatograms",
    ionization_mode="Mixed",
    organism_or_matrix="Representative OpenMS TOPP test traces",
    description="In-tree real chromatogram mzML files packaged with OpenMS for deterministic offline testing.",
    recommended_files=[
        "NoiseFilterSGolay_2_input.chrom.mzML",
        "MRMTransitionGroupPicker_1_input.mzML",
        "OpenSwathWorkflow_1_output.chrom.mzML",
        "OpenSwathWorkflow_13_output.chrom.mzML",
    ],
)


@dataclass
class RawTraceProfile:
    """Quantitative characteristics of a raw chromatogram before smoothing."""
    trace_id: str
    native_id: str
    source_file: str
    num_points: int
    rt_start: float
    rt_end: float
    apex_rt: float
    apex_intensity: float
    baseline_level: float
    noise_sigma_mad: float
    snr: float
    approx_fwhm: float
    intensity_stratum: str  # 'low', 'medium', 'high'
    width_stratum: str      # 'narrow', 'medium', 'broad'
    is_valid: bool = True
    rejection_reason: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class SelectionManifest:
    """Documenting deterministic chromatogram selection without cherry-picking."""
    dataset_accession: str
    source_files: List[str]
    selection_seed: int
    target_sample_size: int
    total_candidates_inspected: int
    total_valid_candidates: int
    total_selected: int
    rejection_counts: Dict[str, int]
    stratum_counts_available: Dict[str, int]
    stratum_counts_selected: Dict[str, int]
    selected_trace_ids: List[str]
    selection_criteria: Dict[str, Any]

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def estimate_mad_noise(y: np.ndarray) -> float:
    """
    Estimate baseline noise using Median Absolute Deviation (MAD) of first differences.
    sigma = median(|dy - median(dy)|) / (0.6745 * sqrt(2))

    Robust against monotonic peak slopes, steps, and baseline drift.
    """
    if len(y) < 4:
        return 0.0
    dy = np.diff(y)
    med_dy = np.median(dy)
    mad = np.median(np.abs(dy - med_dy))
    scale = 0.6745 * math.sqrt(2.0)
    return float(mad / scale) if scale > 0 else 0.0


def profile_raw_chromatogram(
    chrom: pyopenms.MSChromatogram,
    trace_id: str,
    source_file: str,
    min_points: int = 20,
    min_snr: float = 3.0,
) -> RawTraceProfile:
    """
    Profile a raw chromatogram strictly from raw signal values without smoothing.
    Assigns intensity and peak width strata for stratified sampling.
    """
    native_id = chrom.getNativeID() if hasattr(chrom, "getNativeID") else trace_id
    rt_raw, int_raw = chrom.get_peaks()

    rt_arr = np.array(rt_raw, dtype=np.float64)
    int_arr = np.array(int_raw, dtype=np.float64)
    n_pts = len(int_arr)

    if n_pts < min_points:
        return RawTraceProfile(
            trace_id=trace_id,
            native_id=native_id,
            source_file=source_file,
            num_points=n_pts,
            rt_start=float(rt_arr[0]) if n_pts > 0 else 0.0,
            rt_end=float(rt_arr[-1]) if n_pts > 0 else 0.0,
            apex_rt=0.0,
            apex_intensity=0.0,
            baseline_level=0.0,
            noise_sigma_mad=0.0,
            snr=0.0,
            approx_fwhm=0.0,
            intensity_stratum="invalid",
            width_stratum="invalid",
            is_valid=False,
            rejection_reason=f"Insufficient points ({n_pts} < {min_points})",
        )

    # Robust baseline estimation: 10th percentile
    baseline = float(np.percentile(int_arr, 10))
    noise_sigma = estimate_mad_noise(int_arr)
    apex_idx = int(np.argmax(int_arr))
    apex_int = float(int_arr[apex_idx])
    apex_rt = float(rt_arr[apex_idx])

    if noise_sigma <= 1e-9:
        snr = 0.0
    else:
        snr = max(0.0, (apex_int - baseline) / noise_sigma)

    if snr < min_snr:
        return RawTraceProfile(
            trace_id=trace_id,
            native_id=native_id,
            source_file=source_file,
            num_points=n_pts,
            rt_start=float(rt_arr[0]),
            rt_end=float(rt_arr[-1]),
            apex_rt=apex_rt,
            apex_intensity=apex_int,
            baseline_level=baseline,
            noise_sigma_mad=noise_sigma,
            snr=snr,
            approx_fwhm=0.0,
            intensity_stratum="invalid",
            width_stratum="invalid",
            is_valid=False,
            rejection_reason=f"Low SNR ({snr:.1f} < {min_snr:.1f})",
        )

    # Approximate FWHM: span of points above half-height
    half_height = baseline + (apex_int - baseline) / 2.0
    above_half_indices = np.where(int_arr >= half_height)[0]
    if len(above_half_indices) > 1:
        fwhm_rt = float(rt_arr[above_half_indices[-1]] - rt_arr[above_half_indices[0]])
        fwhm_pts = len(above_half_indices)
    else:
        dt = float(np.median(np.diff(rt_arr))) if len(rt_arr) > 1 else 1.0
        fwhm_rt = dt
        fwhm_pts = 1

    # Stratification by intensity
    if apex_int < 1e4:
        int_strat = "low"
    elif apex_int < 1e5:
        int_strat = "medium"
    else:
        int_strat = "high"

    # Stratification by width in points
    if fwhm_pts < 10:
        width_strat = "narrow"
    elif fwhm_pts < 25:
        width_strat = "medium"
    else:
        width_strat = "broad"

    return RawTraceProfile(
        trace_id=trace_id,
        native_id=native_id,
        source_file=source_file,
        num_points=n_pts,
        rt_start=float(rt_arr[0]),
        rt_end=float(rt_arr[-1]),
        apex_rt=apex_rt,
        apex_intensity=apex_int,
        baseline_level=baseline,
        noise_sigma_mad=noise_sigma,
        snr=snr,
        approx_fwhm=fwhm_rt,
        intensity_stratum=int_strat,
        width_stratum=width_strat,
        is_valid=True,
    )


class RealChromatogramSelector:
    """
    Deterministic stratified chromatogram selector.
    Guarantees no cherry-picking:
    - Filters traces purely on raw signal criteria (min points, min SNR)
    - Stratifies by raw intensity and peak width
    - Samples fixed quota per stratum using fixed RNG seed
    - Emits reproducible selection manifest
    """

    def __init__(
        self,
        target_sample_size: int = 100,
        seed: int = 42,
        min_points: int = 20,
        min_snr: float = 3.0,
    ):
        self.target_sample_size = target_sample_size
        self.seed = seed
        self.min_points = min_points
        self.min_snr = min_snr

    def select(
        self,
        chromatograms: List[Tuple[pyopenms.MSChromatogram, str, str]],  # (chrom, trace_id, source_file)
        dataset_accession: str = "REAL_DATA",
    ) -> Tuple[List[pyopenms.MSChromatogram], SelectionManifest, List[RawTraceProfile]]:
        """
        Perform stratified deterministic selection across candidates.
        Returns: (selected_chroms, manifest, selected_profiles)
        """
        rng = np.random.default_rng(self.seed)

        profiles: List[RawTraceProfile] = []
        chrom_map: Dict[str, pyopenms.MSChromatogram] = {}
        rejection_counts: Dict[str, int] = {
            "insufficient_points": 0,
            "low_snr": 0,
            "constant_or_zero": 0,
        }

        for chrom, t_id, src in chromatograms:
            chrom_map[t_id] = chrom
            p = profile_raw_chromatogram(
                chrom, t_id, src, min_points=self.min_points, min_snr=self.min_snr
            )
            profiles.append(p)
            if not p.is_valid:
                reason = p.rejection_reason or ""
                if "Insufficient points" in reason:
                    rejection_counts["insufficient_points"] += 1
                elif "Low SNR" in reason:
                    rejection_counts["low_snr"] += 1
                else:
                    rejection_counts["constant_or_zero"] += 1

        valid_profiles = [p for p in profiles if p.is_valid]

        # Group into 3x3 strata: (intensity_stratum, width_stratum)
        strata: Dict[str, List[RawTraceProfile]] = {}
        for p in valid_profiles:
            key = f"{p.intensity_stratum}_{p.width_stratum}"
            strata.setdefault(key, []).append(p)

        strata_available = {k: len(v) for k, v in strata.items()}

        # Allocate quota evenly across non-empty strata
        n_strata = max(1, len(strata))
        per_stratum_quota = max(1, self.target_sample_size // n_strata)

        selected_profiles: List[RawTraceProfile] = []
        strata_selected: Dict[str, int] = {}

        # Deterministic stratified sampling
        for k in sorted(strata.keys()):
            stratum_list = strata[k]
            # Sort by trace_id before shuffle to ensure byte-level determinism across platforms
            stratum_list.sort(key=lambda x: x.trace_id)
            n_take = min(len(stratum_list), per_stratum_quota)
            if n_take > 0:
                indices = rng.choice(len(stratum_list), size=n_take, replace=False)
                chosen = [stratum_list[i] for i in sorted(indices)]
                selected_profiles.extend(chosen)
                strata_selected[k] = len(chosen)
            else:
                strata_selected[k] = 0

        # If quota is under target and some strata have surplus, fill deterministically
        if len(selected_profiles) < self.target_sample_size:
            remaining_needed = self.target_sample_size - len(selected_profiles)
            pool = []
            selected_ids = {p.trace_id for p in selected_profiles}
            for p in valid_profiles:
                if p.trace_id not in selected_ids:
                    pool.append(p)
            pool.sort(key=lambda x: x.trace_id)
            if pool:
                n_extra = min(len(pool), remaining_needed)
                extra_indices = rng.choice(len(pool), size=n_extra, replace=False)
                extra_chosen = [pool[i] for i in sorted(extra_indices)]
                selected_profiles.extend(extra_chosen)
                for ep in extra_chosen:
                    k = f"{ep.intensity_stratum}_{ep.width_stratum}"
                    strata_selected[k] = strata_selected.get(k, 0) + 1

        selected_chroms = [chrom_map[p.trace_id] for p in selected_profiles]
        source_files = sorted(list({p.source_file for p in profiles}))

        manifest = SelectionManifest(
            dataset_accession=dataset_accession,
            source_files=source_files,
            selection_seed=self.seed,
            target_sample_size=self.target_sample_size,
            total_candidates_inspected=len(profiles),
            total_valid_candidates=len(valid_profiles),
            total_selected=len(selected_profiles),
            rejection_counts=rejection_counts,
            stratum_counts_available=strata_available,
            stratum_counts_selected=strata_selected,
            selected_trace_ids=[p.trace_id for p in selected_profiles],
            selection_criteria={
                "min_points": self.min_points,
                "min_snr": self.min_snr,
                "intensity_strata": ["low (<1e4)", "medium (1e4-1e5)", "high (>=1e5)"],
                "width_strata": ["narrow (<10 pts)", "medium (10-25 pts)", "broad (>=25 pts)"],
            },
        )

        return selected_chroms, manifest, selected_profiles


def load_chromatograms_from_file(
    file_path: str | Path,
    max_chroms: Optional[int] = None,
) -> List[Tuple[pyopenms.MSChromatogram, str, str]]:
    """
    Load native chromatograms from an mzML file using pyopenms.
    Returns: list of (MSChromatogram, trace_id, source_file)
    """
    path = Path(file_path)
    if not path.exists():
        return []

    exp = pyopenms.MSExperiment()
    try:
        pyopenms.MzMLFile().load(str(path), exp)
    except Exception:
        return []

    chroms = exp.getChromatograms()
    n_load = len(chroms) if max_chroms is None else min(len(chroms), max_chroms)

    res = []
    for i in range(n_load):
        c = chroms[i]
        native_id = c.getNativeID() if hasattr(c, "getNativeID") and c.getNativeID() else f"chrom_{i}"
        trace_id = f"{path.stem}_{native_id}"
        res.append((c, trace_id, path.name))
    return res


def extract_ms1_xic(
    exp: pyopenms.MSExperiment,
    target_mz: float,
    ppm: float = 15.0,
    trace_id: str = "xic",
    source_name: str = "spectra",
) -> Optional[Tuple[pyopenms.MSChromatogram, str, str]]:
    """
    Extract a single MS1 Extracted Ion Chromatogram (XIC) from an MSExperiment.
    Uses fixed ppm mass tolerance window.
    """
    tol = target_mz * (ppm * 1e-6)
    rts: List[float] = []
    ints: List[float] = []

    for spec in exp:
        if spec.getMSLevel() == 1:
            rt = spec.getRT()
            mz_arr, int_arr = spec.get_peaks()
            mask = (mz_arr >= target_mz - tol) & (mz_arr <= target_mz + tol)
            intensity = float(np.sum(int_arr[mask])) if np.any(mask) else 0.0
            rts.append(rt)
            ints.append(intensity)

    if not rts:
        return None

    chrom = pyopenms.MSChromatogram()
    chrom.setName(f"{trace_id}_mz_{target_mz:.4f}")
    chrom.set_peaks((rts, ints))
    return chrom, f"{source_name}_{trace_id}", source_name


def load_swath_windows(windows_file: str | Path) -> List[Tuple[float, float]]:
    """Load SWATH precursor isolation windows from TSV."""
    import csv
    windows: List[Tuple[float, float]] = []
    with open(windows_file, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            start = float(row.get("start") or row.get("Start") or row.get("lower"))
            end = float(row.get("end") or row.get("End") or row.get("upper"))
            windows.append((start, end))
    return windows


def extract_swath_chromatograms(
    mzml_path: str | Path,
    windows_file: str | Path,
    assay_tsv: str | Path,
    output_chrom_mzml: str | Path,
    mz_tolerance_da: float = 0.05,
    max_spectra: Optional[int] = None,
) -> int:
    """
    Extract transition chromatograms from a SWATH mzML file against an assay library.
    Guarantees no circular peak-scoring bias: extracts raw chromatogram traces directly from MS2 scans.
    """
    import csv
    windows = load_swath_windows(windows_file)

    # Load transitions binned by swath isolation window
    by_window: Dict[int, List[Dict[str, Any]]] = {i: [] for i in range(len(windows))}
    all_transitions: List[Dict[str, Any]] = []

    with open(assay_tsv, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            prec_mz = float(row["PrecursorMz"])
            prod_mz = float(row["ProductMz"])
            tr_name = row.get("transition_name", f"tr_{prec_mz}_{prod_mz}")
            exp_rt = float(row.get("Tr_recalibrated", row.get("RetentionTime", 0.0)))
            pep = row.get("PeptideSequence", "")
            prot = row.get("ProteinName", "")

            matched_win = -1
            for widx, (w_start, w_end) in enumerate(windows):
                if w_start <= prec_mz <= w_end:
                    matched_win = widx
                    break

            if matched_win != -1:
                item = {
                    "transition_id": tr_name,
                    "precursor_mz": prec_mz,
                    "product_mz": prod_mz,
                    "expected_rt": exp_rt,
                    "peptide_sequence": pep,
                    "protein_name": prot,
                    "window_idx": matched_win,
                }
                by_window[matched_win].append(item)
                all_transitions.append(item)

    trans_rts: Dict[str, List[float]] = {t["transition_id"]: [] for t in all_transitions}
    trans_ints: Dict[str, List[float]] = {t["transition_id"]: [] for t in all_transitions}

    # Precompute numpy sorted product m/zs per window for fast binary search
    window_fast_lookup = {}
    for widx, tr_list in by_window.items():
        if tr_list:
            prod_mzs = np.array([t["product_mz"] for t in tr_list], dtype=np.float64)
            tr_ids = [t["transition_id"] for t in tr_list]
            window_fast_lookup[widx] = (prod_mzs, tr_ids)

    input_path = Path(mzml_path)
    is_gz = str(input_path).lower().endswith(".gz")
    target_mzml = input_path
    temp_uncompressed: Optional[Path] = None

    if is_gz:
        # For .gz inputs (e.g. .mzML.gz), decompress to uncompressed indexed mzML before openFile
        uncompressed_candidate = input_path.with_suffix("")
        if uncompressed_candidate.exists():
            target_mzml = uncompressed_candidate
        else:
            partial = uncompressed_candidate.with_name(uncompressed_candidate.name + ".partial")
            try:
                try:
                    with gzip.open(input_path, "rb") as f_in, open(partial, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out, length=64 * 1024 * 1024)
                    os.replace(partial, uncompressed_candidate)
                    target_mzml = uncompressed_candidate
                finally:
                    if partial.exists():
                        try:
                            partial.unlink()
                        except OSError:
                            pass
            except (OSError, PermissionError):
                tf = tempfile.NamedTemporaryFile(suffix=".mzML", delete=False)
                temp_uncompressed = Path(tf.name)
                tf.close()
                try:
                    with gzip.open(input_path, "rb") as f_in, open(temp_uncompressed, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out, length=64 * 1024 * 1024)
                    target_mzml = temp_uncompressed
                except Exception:
                    if temp_uncompressed.exists():
                        try:
                            temp_uncompressed.unlink()
                        except OSError:
                            pass
                    raise

    try:
        # Open with OnDiscMSExperiment or fallback to MSExperiment (only for uncompressed inputs)
        use_ondisc = False
        od = pyopenms.OnDiscMSExperiment()
        try:
            if od.openFile(str(target_mzml)):
                use_ondisc = True
                n_spec = od.getNrSpectra()
        except Exception:
            use_ondisc = False

        if not use_ondisc:
            if is_gz:
                raise RuntimeError(
                    f"Failed to open decompressed mzML '{target_mzml}' with OnDiscMSExperiment. "
                    "In-memory MzMLFile().load fallback is rejected for compressed inputs to prevent memory exhaustion."
                )
            exp = pyopenms.MSExperiment()
            pyopenms.MzMLFile().load(str(target_mzml), exp)
            n_spec = exp.getNrSpectra()
            get_spec_fn = lambda idx: exp.getSpectrum(idx)
        else:
            get_spec_fn = lambda idx: od.getSpectrum(idx)

        limit = n_spec if max_spectra is None else min(n_spec, max_spectra)

        for i in range(limit):
            spec = get_spec_fn(i)
            if spec.getMSLevel() == 2:
                rt = spec.getRT()
                precursors = spec.getPrecursors()
                if not precursors:
                    continue
                prec_mz = precursors[0].getMZ()

                matched_win = -1
                for widx, (w_start, w_end) in enumerate(windows):
                    if w_start <= prec_mz <= w_end:
                        matched_win = widx
                        break

                if matched_win in window_fast_lookup:
                    prod_mzs, tr_ids = window_fast_lookup[matched_win]
                    s_mz, s_int = spec.get_peaks()
                    if len(s_mz) > 0:
                        for target_mz, tr_id in zip(prod_mzs, tr_ids):
                            idx_left = int(np.searchsorted(s_mz, target_mz - mz_tolerance_da))
                            idx_right = int(np.searchsorted(s_mz, target_mz + mz_tolerance_da, side="right"))
                            val = float(np.sum(s_int[idx_left:idx_right])) if idx_right > idx_left else 0.0
                            trans_rts[tr_id].append(rt)
                            trans_ints[tr_id].append(val)

        out_exp = pyopenms.MSExperiment()
        for t in all_transitions:
            rts = trans_rts[t["transition_id"]]
            ints = trans_ints[t["transition_id"]]
            if len(rts) > 0:
                c = pyopenms.MSChromatogram()
                c.setNativeID(t["transition_id"])
                c.setName(f"{t['peptide_sequence']}_{t['transition_id']}")
                c.set_peaks((rts, ints))
                c_prec = pyopenms.Precursor()
                c_prec.setMZ(t["precursor_mz"])
                c.setPrecursor(c_prec)
                c_prod = pyopenms.Product()
                c_prod.setMZ(t["product_mz"])
                c.setProduct(c_prod)
                out_exp.addChromatogram(c)

        Path(output_chrom_mzml).parent.mkdir(parents=True, exist_ok=True)
        pyopenms.MzMLFile().store(str(output_chrom_mzml), out_exp)
        return out_exp.getNrChromatograms()
    finally:
        if temp_uncompressed and temp_uncompressed.exists():
            try:
                temp_uncompressed.unlink()
            except OSError:
                pass
