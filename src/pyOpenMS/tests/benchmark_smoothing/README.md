# OpenMS Chromatogram Smoothing Benchmark (Remediated)

> **Benchmarking `ModifiedSincSmoother` against `SavitzkyGolayFilter` on Chromatograms**  
> **OpenMS GitHub Issue #10425**

---

## 1. Motivation and Background

Smoothing is an essential preprocessing step in LC-MS chromatogram analysis (extracted ion chromatograms, MRM/SRM traces, and SWATH-MS ion traces). Traditionally, the Savitzky–Golay filter has been the default standard across mass spectrometry software.

OpenMS recently incorporated `ModifiedSincSmoother` based on:
> *Schmid, Rath & Diebold, "Why and How Savitzky-Golay Filters Should Be Replaced", ACS Measurement Science Au, 2022, 2(2), 185-196.*

This benchmark framework rigorously evaluates whether `ModifiedSincSmoother` improves chromatogram quality compared to `SavitzkyGolayFilter`, addressing OpenMS Issue #10425.

---

## 2. Benchmark Architecture

```
src/pyOpenMS/tests/benchmark_smoothing/
├── __init__.py                     # Package initialization
├── synthetic_data.py               # Deterministic synthetic chromatograms (Gaussian, EMG, Physical Noise)
├── metrics.py                      # Sub-scan parabolic apex & flank FWHM metrics, noise, doublet area
├── real_data.py                    # Verified metadata, raw profiling, deterministic stratified selector, manifests
├── real_metrics.py                 # Real-data empirical metrics (apex shift, area/width ratio, flank MAD noise, replicate CV)
├── benchmark_engine.py             # Execution engine, parameter tuning, phase sweeps, real-data benchmark
├── plotter.py                      # Publication-grade figure generation (Matplotlib Agg)
├── run_benchmark.py                # Standalone CLI entry point and markdown reporter
├── test_benchmark.py               # Comprehensive unit & integration regression test suite (23 tests)
└── README.md                       # Comprehensive methodology documentation
```

---

## 3. Scientific Methodology & Audit Remediations

1. **Parameter Tuning without Test Leakage**:
   - Optimal default parameters are tuned strictly on a separate calibration trace (`tuning_composite`, generated with seed `seed + 10000`).
   - **Tuning Objective**: Minimizes global Root Mean Squared Error (RMSE) against the ground-truth noiseless signal. Global RMSE was chosen as an objective, unbiased $L_2$ metric that simultaneously balances noise suppression across flat baseline regions and peak preservation across diverse shapes without arbitrary heuristic weightings.
   - The selected parameters are frozen and evaluated on completely held-out test chromatograms.
2. **Sub-Scan Peak Fidelity**:
   - Peak apex RT and intensity are determined via 3-point parabolic interpolation around the maximum, eliminating discrete grid quantization.
   - Full Width at Half Maximum (FWHM) is computed continuously using linear flank interpolation at half-maximum ($I_{baseline} + (I_{apex} - I_{baseline}) / 2$).
3. **Narrow Peak Undersampling**:
   - The ~3-point FWHM peak ($FWHM = 3.0$ s with sampling interval $\Delta t = 1.0$ s) is severely undersampled, containing only ~3 data points across its half-maximum width. It deliberately serves as an extreme stress test of peak shape preservation.
   - Under these conditions, convolution smoothing inevitably attenuates peak height and broadens width. The observed ~58% FWHM broadening for Modified Sinc and ~73% for Savitzky–Golay reflect actual physical filter dispersion rather than measurement artifacts.
   - Evaluated across 4 sub-scan sampling phases ($0.0, 0.25, 0.50, 0.75 \times \Delta t$) to assess grid-phase jitter effects.
4. **Overlapping Doublet Evaluation**:
   - For co-eluting doublets, single-peak FWHM is not physically defined due to blended flanks.
   - Evaluated via total doublet area preservation compared against total ground truth area ($\text{Area}_1 + \text{Area}_2$), valley-to-peak ratio (VPR), and individual apex resolvability.
5. **Physical Non-Rectified Noise**:
   - Synthetic chromatograms incorporate detector baseline offsets ($B_0 \ge 80\text{--}150$) so zero-mean Gaussian noise $\mathcal{N}(0, \sigma^2)$ is not artificially clipped at zero, preserving Gaussian noise statistics.
6. **Realistic Downstream Peak Picking**:
   - Evaluated with OpenMS `PeakPickerHiRes` directly on smoothed traces across $S/N \in [1.0, 2.0, 3.0, 5.0]$.
   - $S/N = 2.0$ serves as the primary operating point (standard threshold for distinguishing real chromatographic peaks from baseline noise).
7. **Frequency-Matched Bandwidth Comparison**:
   - Uses OpenMS's built-in `ModifiedSincSmoother.savitzkyGolayBandwidth(p, m_SG)` to determine the exact 3-dB cutoff frequency of a Savitzky–Golay filter, and `ModifiedSincSmoother.bandwidthToM(is_ms1, degree, bandwidth)` to compute the matching Modified Sinc parameter.
   - Serves as a controlled secondary comparison at identical frequency cutoffs, alongside practical independent parameter tuning.
8. **Runtime Microbenchmarking & Algorithmic Buffer Analysis**:
   - Uses warmup iterations, repeated loops of 50–500 iterations, and `gc.disable()`, matching `benchmark_pyopenms.py` standards.
   - Direct OS RSS is not reported because page granularity (4 KB) cannot resolve small temporary C++ heap allocations. Instead, algorithmic buffer allocation from `ModifiedSincSmoother.cpp` is analyzed.
9. **Real-Data Checks**:
   - Sanity checks on representative OpenMS test chromatograms (`src/tests/topp/`) confirm stable numerical execution without divergence. Larger external LC-MS datasets with orthogonal gold standards would be needed for broader real-world validation.

---

## 4. Real-Data Empirical Characterization Layer (Sanity Check)

In response to maintainer feedback (*"Would be great if this is tested on real chroms"*), an empirical characterization layer is integrated to verify whether the relative behaviors observed under controlled synthetic conditions translate to experimental LC-MS chromatograms.

### A. Scientific Guardrails & Methodological Principles
- **No Ground Truth on Real Data**: Unlike synthetic signals with an analytical ground truth, real LC-MS chromatograms possess no known "true" signal. Consequently, RMSE and ground-truth error metrics are strictly avoided to eliminate circular bias.
- **Role as Empirical Sanity Check**: Real data serves exclusively as an empirical sanity check, not as "proof of superiority" or "statistically comprehensive validation."
- **Fair Comparison**: Both smoothers receive the exact same raw chromatogram. Default and bandwidth-matched parameters are evaluated without test-set tuning.
- **Zero Cherry-Picking**: Chromatogram selection is completely automated and deterministic based purely on **raw signal** characteristics before smoothing.

### B. Verified Public Reference Datasets

| Dataset | Accession / Repository | Platform & Acquisition | Sample / Matrix | Downloadable Files & Status |
|---|---|---|---|---|
| **Dataset A (SWATH/DIA)** | [`PASS00779`](http://www.peptideatlas.org/PASS/PASS00779)<br>PeptideAtlas / PASSEL | AB SCIEX TripleTOF 5600<br>SWATH-MS DIA (32 × 25 Da windows) | *Mycobacterium tuberculosis*<br>(Wayne dormancy model) | 3 mzML.gz files (R1, R2, R3; 3.3 GB compressed, ~10 GB uncompressed each); Assay library (`Mtb_TubercuList-R27_iRT_UPS.tsv`, 48 MB); SWATH windows file; Open access via FTP (`ftp://PASS00779:SWATH@ftp.peptideatlas.org/`). |
| **Dataset B (Metabolomics)** | [`MTBLS404`](https://www.ebi.ac.uk/metabolights/MTBLS404)<br>MetaboLights | Thermo LTQ-Orbitrap Discovery<br>LC-HRMS Full-Scan MS1 (negative ESI) | Human urine cohort<br>(184 adults + 26 pooled-QCs) | 234 mzML files (18 GB total) on EBI FTP (`ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS404/`); Manageable subset: 4–6 pooled-QC mzML files (`QC01.mzML`, `QC02.mzML`, ...); EMBL-EBI Terms of Use. |
| **In-Tree Fixtures** | `LOCAL_FIXTURES`<br>OpenMS In-Tree | Diverse (SCIEX TripleTOF, QTRAP, Thermo)<br>MRM / SWATH / DDA chromatograms | Standard OpenMS test traces | Packaged in `src/tests/topp/` (`NoiseFilterSGolay_2_input.chrom.mzML`, `MRMTransitionGroupPicker_1_input.mzML`, `OpenSwathWorkflow_1_output.chrom.mzML`); Hermetic offline test execution. |

> **Note on Data Volume**: Full public DIA and metabolomics datasets span tens of gigabytes. CI test runs must NEVER download multi-gigabyte datasets over the network. The benchmark executes hermetically on in-tree fixtures by default, with opt-in support for downloaded dataset directories via `--real-data-dir`.

### C. Extraction Pipeline
- **SWATH-MS / DIA (PASS00779)**: Extracted using OpenMS `OpenSwathWorkflow` with `-out_chrom` or native `MSExperiment` chromatogram parsing. Precursors can be deterministically subsampled (500–1000 targets) from `Mtb_TubercuList-R27_iRT_UPS.tsv` using a fixed seed.
- **Metabolomics MS1 XICs (MTBLS404)**: Extracted via `extract_ms1_xic` with fixed 15 ppm mass tolerance from centroided MS1 spectra across pooled-QC replicate runs.

### D. Deterministic Stratified Selection (No Cherry-Picking)
Chromatograms are selected without human intervention or post-smoothing filtering:
1. **Raw Profiling**: Traces are profiled strictly on raw signal (minimum points $\ge 20$, raw $S/N \ge 3.0$).
2. **3×3 Stratification**: Traces are binned by:
   - **Intensity**: Low ($< 10^4$), Medium ($10^4\text{--}10^5$), High ($\ge 10^5$)
   - **Width**: Narrow ($< 10$ points), Medium ($10\text{--}25$ points), Broad ($\ge 25$ points)
3. **Deterministic Sampling**: A fixed quota is drawn from each stratum using a fixed seed (`seed=42`).
4. **Manifest Generation**: An immutable `real_data_manifest.json` is exported, recording dataset accession, file sources, seed, criteria, rejection counts, and trace IDs.

### E. Real-Data Empirical Metrics
- **Signal Preservation**:
  - *Apex RT Shift*: Sub-scan parabolic apex shift ($\Delta RT = RT_{\text{smooth}} - RT_{\text{raw}}$).
  - *Apex Intensity Ratio*: $I_{\text{smooth}} / I_{\text{raw}}$.
  - *Peak Area Ratio*: Numerical trapezoidal area preservation ($\text{Area}_{\text{smooth}} / \text{Area}_{\text{raw}}$).
  - *FWHM Width Ratio*: Continuous flank-interpolated width ratio ($\text{FWHM}_{\text{smooth}} / \text{FWHM}_{\text{raw}}$).
- **Baseline Flank Noise**:
  - Confined strictly to peak-free outer flank regions (first and last 20% of points).
  - Robust Median Absolute Deviation (MAD) of first differences divided by $\sqrt{2}$:
    $$\sigma_{\text{MAD}} = \frac{\text{median}(|\Delta y - \text{median}(\Delta y)|)}{0.6745 \cdot \sqrt{2}}$$
    This formulation is immune to monotonic chromatographic slopes, steps, and baseline drift.
- **Overlapping Resolvability**:
  - Number of local maxima and Valley-to-Peak Ratio (VPR).
- **Downstream Peak Detection**:
  - OpenMS `PeakPickerHiRes` at $S/N = 2.0$ to verify peak count concordance.
- **QC Replicate Consistency (MTBLS404)**:
  - Coefficient of Variation (CV% = $\frac{\sigma}{\mu} \times 100\%$) of apex RT, apex intensity, and area across repeated injections.

---

## 5. How to Reproduce

### Running Unit Tests (Offline, Fast, Hermetic)
```bash
python -m unittest discover -s src/pyOpenMS/tests/benchmark_smoothing -p "test_*.py" -v
```

### Running Quick Synthetic Benchmark
```bash
python src/pyOpenMS/tests/benchmark_chromatogram_smoothing.py --quick
```

### Running Full Synthetic Benchmark with Markdown Report
```bash
python src/pyOpenMS/tests/benchmark_chromatogram_smoothing.py --seed 42 --report
```

### Running Real-Data Empirical Characterization (In-Tree Fixtures, Offline)
```bash
python src/pyOpenMS/tests/benchmark_chromatogram_smoothing.py --real-data --report
```

### Running Real-Data Only (Skip Synthetic Benchmark)
```bash
python src/pyOpenMS/tests/benchmark_chromatogram_smoothing.py --real-data-only --output-dir results/
```

### Running on Downloaded Public Datasets (Opt-In)
```bash
# PASS00779 (M. tuberculosis SWATH-MS)
python src/pyOpenMS/tests/benchmark_chromatogram_smoothing.py --real-data --real-dataset pass00779 --real-data-dir /path/to/PASS00779/

# MTBLS404 (Sacurine Urine Metabolomics)
python src/pyOpenMS/tests/benchmark_chromatogram_smoothing.py --real-data --real-dataset mtbls404 --real-data-dir /path/to/MTBLS404/
```


