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
├── benchmark_engine.py             # Execution engine, parameter tuning, phase sweeps, microbenchmarks
├── plotter.py                      # Publication-grade figure generation (Matplotlib Agg)
├── run_benchmark.py                # Standalone CLI entry point and markdown reporter
├── test_benchmark.py               # Unit & integration regression test suite (14 tests)
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

## 4. How to Reproduce

### Running Unit Tests
```bash
python -m unittest discover -s src/pyOpenMS/tests/benchmark_smoothing -p "test_*.py" -v
```

### Running Quick Benchmark
```bash
python src/pyOpenMS/tests/benchmark_chromatogram_smoothing.py --quick
```

### Running Full Comprehensive Benchmark
```bash
python src/pyOpenMS/tests/benchmark_chromatogram_smoothing.py --seed 42 --report
```

