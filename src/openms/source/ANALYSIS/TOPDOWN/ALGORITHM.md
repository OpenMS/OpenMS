# FLASHDeconv Algorithm Description

FLASHDeconv is an ultrafast mass spectral deconvolution algorithm for top-down
proteomics. It converts multiply-charged m/z spectra into zero-charge (neutral)
mass spectra. Its key innovation is a **log-space binning trick** that turns
charge-state determination into efficient integer bin-offset lookups.

Reference: Jeong et al., "FLASHDeconv: Ultrafast, High-Quality Feature
Deconvolution for Top-Down Proteomics," *Cell Systems*, 2020.

## Pipeline Overview

Entry point: `FLASHDeconvAlgorithm::run()` (`FLASHDeconvAlgorithm.cpp`)

```
Input MS data
    │
    ▼
1. Initialization
   ├── Pre-calculate averagine isotope distributions
   ├── Auto-estimate mass tolerances (if not specified)
   └── Set up target/decoy generators
    │
    ▼
2. Spectral Deconvolution (per spectrum)
   ├── 2a. Log-m/z binning
   ├── 2b. Candidate mass bin selection via universal pattern
   ├── 2c. Peak group construction (isotope clustering)
   └── 2d. Scoring and filtering
    │
    ▼
3. Feature Tracing (across retention time)
   └── Link deconvolved masses into chromatographic features
    │
    ▼
4. Quality Assessment
   ├── Target-decoy FDR estimation
   └── Optional isobaric quantification
    │
    ▼
Output: deconvolved spectra, mass features, quality scores
```

## Stage 1: Initialization

### Pre-calculated Averagine

Class: `FLASHHelperClasses::PrecalculatedAveragine` (`FLASHHelperClasses.h`)

Before deconvolution, theoretical isotope distributions are pre-computed for
mass bins spanning the full mass range (default 50–100,000 Da) at ~25 Da
intervals. For each bin the algorithm stores:

- The isotope distribution
- L2 norms (for fast cosine similarity)
- Apex index (most abundant isotope position)
- Left/right isotope count bounds from the apex
- SNR multiplication factors

This avoids recalculating isotope patterns per candidate — a major speedup.

### Tolerance Estimation

If tolerances are not user-specified, FLASHDeconv fits a Gaussian to a
histogram of mass errors from an initial pass to auto-estimate them
(`FLASHDeconvAlgorithm::determineTolerance_()`).

## Stage 2: Spectral Deconvolution

Core implementation: `SpectralDeconvolution` (`SpectralDeconvolution.cpp`)

### 2a. Log-m/z Binning

Function: `updateLogMzPeaks_()`, `binLogMzPeaks_()`

Raw peaks are transformed: `logMz = log(mz - charge_carrier_mass)`.
The log-m/z axis is discretized into bins of width `tolerance / 2.5`.
A `boost::dynamic_bitset` records which bins contain peaks for O(1)
occupancy checks.

**Why log-space?** For a mass M at charge z: `log(M/z) = log(M) - log(z)`.
All charge states of the same mass map to the same mass bin when offset by
`-log(z)`. This converts multiplication into addition.

### 2b. Candidate Mass Bin Selection

Function: `updateCandidateMassBins_()` — determines FLASHDeconv's runtime.

The **universal pattern** is: `universal_pattern[j] = -log(j + 1)` for
charge j+1. For each occupied m/z bin, the algorithm checks mass bins at
the universal pattern offsets:

```
mass_bin_index = mz_bin_index + binned_universal_pattern_[j]
```

Evidence accumulates for each candidate mass by counting consecutive charge
states. Three heuristic filters prune candidates early:

1. **Charge continuity**: Consecutive charge states must be present.
2. **Intensity ratio**: Ratio between consecutive charges must be within
   10× (low charge) or 5× (high charge).
3. **Isotope peak verification**: For low charges (z ≤ 10), adjacent
   isotope peaks must exist.

#### Harmonic Artifact Reduction

A harmonic pattern matrix for small primes {2, 3, 5, 7, 11} identifies
masses that are harmonic artifacts (e.g., M/2 falsely detected as M due to
charge-state misassignment). Candidate mass bins scoring better as harmonics
of a higher-quality mass are eliminated.

### 2c. Peak Group Construction

Function: `getCandidatePeakGroups_()`

For each surviving candidate mass bin with its determined charge range:

1. For each charge state, find the most intense peak matching the expected
   m/z.
2. From that anchor, search both directions for **isotope peaks** at
   spacing `iso_delta = ~1.00235 Da / charge`.
3. The search window is bounded by pre-calculated `left_count_from_apex`
   and `right_count_from_apex` from the averagine model.
4. Assign isotope indices and compute the monoisotopic mass.

Output: a `PeakGroup` — the collection of raw peaks assigned to one
monoisotopic mass.

### 2d. Scoring and Filtering

Function: `scoreAndFilterPeakGroups_()` (OpenMP-parallelized)

**Iterative isotope pattern matching**:
1. Recruit all peaks matching the candidate mass.
2. Compute **cosine similarity** between observed and theoretical
   (averagine) isotope envelopes.
3. Check if the monoisotopic assignment is off by ±N isotopes.
4. Iterate (up to 10×, early-stop when cosine decreases).

**Quality score (Qscore)** incorporates:
- Isotope cosine similarity (default minimum: 0.85)
- Per-charge signal-to-noise ratio (default minimum: 0.25)
- Overall SNR

**Filtering**: Peak groups are rejected if Qscore ≤ 0, mass is out of
range, or high-charge groups have insufficient consecutive charge states.

**Final harmonic removal**: A second pass checks each peak group against its
harmonic multiples. If a harmonic (e.g., 2× mass) has better SNR, the
candidate is discarded.

## Stage 3: Feature Tracing

Class: `MassFeatureTrace`, called by `FLASHDeconvAlgorithm::runFeatureFinding_()`

1. Deconvolved masses from individual spectra are traced across retention
   time using `MassTraceDetection`.
2. A **2D quality score** combines per-spectrum isotope pattern quality (1D)
   with chromatographic peak shape (temporal continuity).
3. For MSn data, features are grouped by precursor mass, maintaining the
   MS1 → MS2 → MS3 hierarchy.

## Stage 4: Quality Assessment

### Target-Decoy FDR

Class: `Qvalue`

Two decoy types:
- **Noise decoys**: Use perturbed isotope spacing (0.9444 Da instead of
  ~1.00235 Da). These nonsensical patterns should not match real signals.
- **Signal decoys**: Target masses are excluded from the decoy run; any
  mass found represents a false positive among remaining signals.

FDR is computed by `Qvalue::updatePeakGroupQvalues()` and propagated to
feature-level quality scores.

### Isobaric Quantification

Class: `TopDownIsobaricQuantification`

Optional quantification of isobaric species within the deconvolved data.

## MSn Handling

For MS2+ spectra (`findPrecursorPeakGroupsForMSnSpectra_()`):
- Mass range is constrained by the precursor ± isolation window.
- Charge range is constrained by the precursor charge.
- The algorithm searches MS(n-1) spectra to find the best precursor peak
  group, selecting the one with highest charge SNR.
- Optional spectrum merging combines spectra by precursor mass.

## Key Data Structures

| Structure | File | Purpose |
|---|---|---|
| `PrecalculatedAveragine` | `FLASHHelperClasses.h` | Pre-computed isotope distributions for fast lookup |
| `LogMzPeak` | `FLASHHelperClasses.h` | Log-transformed peak with charge/isotope annotations |
| `PeakGroup` | `PeakGroup.h` | Collection of peaks assigned to one monoisotopic mass |
| `DeconvolvedSpectrum` | `DeconvolvedSpectrum.h` | Container of PeakGroups from one spectrum |
| `MassFeature` | `FLASHHelperClasses.h` | Traced mass across retention time |

## What Makes It Fast

| Technique | Effect |
|---|---|
| Log-space binning | Converts charge detection from O(peaks × charges) to O(bins) integer addition |
| Bitset occupancy tracking | O(1) membership tests via `boost::dynamic_bitset` |
| Pre-computed averagine | Eliminates redundant isotope distribution calculations |
| Charge continuity + intensity ratio heuristics | Prunes implausible candidates before scoring |
| Harmonic pre-filtering | Removes artifacts without full re-scoring |
| Iterative scoring with early stopping | Converges in 2–3 iterations typically |
| OpenMP parallelization | Parallel scoring of independent peak groups |

The fundamental insight: in log-space, `m/z = M/z` becomes
`log(m/z) = log(M) - log(z)`, so charge assignment reduces to checking
fixed bin offsets — far cheaper than multiplicative mass calculation in
linear space.
