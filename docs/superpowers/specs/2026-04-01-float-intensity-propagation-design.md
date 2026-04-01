# Propagate float IntensityType Through Algorithms and I/O

**Issue:** [#8872](https://github.com/OpenMS/OpenMS/issues/8872)
**Date:** 2026-04-01
**Status:** Design

## Context

PR #8857 changed `ChromatogramPeak::IntensityType` from `double` to `float`, aligning it
with `Peak1D`, `Peak2D`, and `MobilityPeak2D`. All kernel peak types now store intensity as
`float`. However, ~30 locations in algorithms and I/O still use `double` for intensity data,
negating the memory and vectorization benefits.

## Goals

- **Maximum memory reduction**: Propagate `float` for intensity storage everywhere possible.
- **SIMD-friendly**: Uniform `float` types (including accumulators) enable clean
  auto-vectorization. AVX2 processes 8 floats vs 4 doubles per register.
- **No file format changes**: Binary wire formats stay `double` for backwards compatibility.
- **No OpenSWATH interface boundary change**: `BinaryDataArray::data` stays `vector<double>`.

## Design Decisions

### Type policy

- Use `Peak1D::IntensityType` (or the relevant peak class typedef) when the variable
  semantically represents a peak intensity and the enclosing code already includes the peak type.
- Use `float` directly for local temporaries or contexts where pulling in a peak header
  would add unnecessary coupling.
- Coordinate types (m/z, RT, IM) remain `double`.

### Accumulators

All-float, including loop accumulators and running sums. Uniform types enable clean SIMD
vectorization. Mixed float-input/double-accumulator patterns force `vcvtps2pd` conversion
instructions that halve throughput and can prevent auto-vectorization entirely.

For typical OpenSWATH data sizes (hundreds to low thousands of points per chromatogram),
float accumulation precision is adequate. If a specific algorithm shows precision issues,
pairwise summation can be added as a targeted fix.

### Cross-correlation return values

`XCorrArrayType` (`vector<pair<int, double>>`) stays as-is. The correlation coefficients
are computed results, not intensity storage. They flow into score scalars written to
TSV/OSW output.

### Boundaries that stay `double`

| Boundary | Reason |
|----------|--------|
| `BinaryDataArray::data` (`vector<double>`) | OpenSWATH interface; cascading risk too high |
| `DataAccessHelper.cpp` conversions | Bridge the float/double boundary at BinaryDataArray |
| `MzMLSqliteHandler` wire format | Backwards compatibility with existing `.sqMass` files |
| `CachedMzMLHandler` wire format | Backwards compatibility with cached `.mzML` files |
| `NonNegativeLeastSquaresSolver` | Fortran-derived NNLS solver requires `double*`; channel arrays are small (~10 elements), savings negligible |
| `CrawdadWrapper::SetChromatogram()` | External optional dependency; convert at call site |
| `OpenSwath_Scores` scalar fields (~42 doubles) | Computed score results, not intensity storage arrays; no memory benefit from changing scalars |

## PR Structure

Five PRs merged in sequence, ordered by risk tier.

### PR 0 -- Template OpenSwathAlgo for float

Template the scoring and stats helper functions in `openswathalgo` to accept both `float`
and `double`. This is the foundation PR that unblocks PR 2.

**Scoring.h / Scoring.cpp** -- 9 functions to template on scalar type `T`:

| Function | Notes |
|----------|-------|
| `standardize_data(vector<T>&)` | Uses `Eigen::Matrix<T,Dynamic,1>` instead of `VectorXd` |
| `normalizedCrossCorrelation(vector<T>&, vector<T>&, ...)` | Calls standardize_data + Post |
| `normalizedCrossCorrelationPost(vector<T>&, vector<T>&, ...)` | Calls calculateCrossCorrelation |
| `calculateCrossCorrelation(const vector<T>&, const vector<T>&, ...)` | `Eigen::Map<const Matrix<T,Dynamic,1>>`, dot products |
| `calcxcorr_legacy_mquest_(vector<T>&, vector<T>&, bool)` | Legacy cross-correlation |
| `normalize_sum(vector<T>&)` | Simple arithmetic |
| `NormalizedManhattanDist(vector<T>&, vector<T>&)` | Currently non-templated overload |
| `RootMeanSquareDeviation(const vector<T>&, const vector<T>&)` | Currently non-templated overload |
| `SpectralAngle(const vector<T>&, const vector<T>&)` | No Eigen dependency |

Return type stays `XCorrArrayType` (`vector<pair<int, double>>`) -- correlation coefficients
remain double precision.

**StatsHelpers.h / StatsHelpers.cpp** -- 3 functions:

| Function | Notes |
|----------|-------|
| `normalize(const vector<T>&, double, vector<T>&)` | Normalization; normalizer stays `double` (it may be a sum/norm computed at higher precision) |
| `dotprodScoring(vector<T>, vector<T>)` | Dot product scoring |
| `manhattanScoring(vector<T>, vector<T>)` | Manhattan scoring |

**Implementation approach**: Template on scalar type `T`, use `Eigen::Matrix<T, Dynamic, 1>`
internally. Add explicit template instantiations for `float` and `double` in the `.cpp` files.
This matches the existing pattern -- `StatsHelpers.cpp` already does explicit instantiation
for `dotProd<float>`, `dotProd<double>`.

**Testing**: Existing openswathalgo unit tests must pass for both float and double
instantiations. Add float-specific test cases.

### PR 1 -- Cleanup (zero behavioral change)

Remove redundant casts and unnecessary `double` promotions where `getIntensity()` already
returns `float`.

| File | Change |
|------|--------|
| `PROCESSING/RESAMPLING/LinearResampler.h:106` | Remove `static_cast<double>(getIntensity())` |
| `ANALYSIS/OPENSWATH/MRMScoring.cpp:556` | Remove `static_cast<double>(getIntensity())` |
| `ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h:491` | `const double` -> `const auto` / `IntensityType` |
| `QUANTITATION/IsotopeLabelingMDVs.cpp:204` | Remove `(double)(it->getIntensity())` C-style cast |
| `FORMAT/QcMLFile.cpp:1274` | `double sum` -> `float` / `IntensityType` |
| `PROCESSING/SMOOTHING/GaussFilter.cpp:62` | Remove `static_cast<double>` (prepares for PR 2) |

**Testing**: All existing unit tests pass unchanged.

### PR 2 -- Algorithm internals (main memory savings)

Change intensity storage vectors from `vector<double>` to `vector<float>` (via typedef where
appropriate) in scoring, peak picking, and signal processing.

| File | Change |
|------|--------|
| `ANALYSIS/OPENSWATH/MRMScoring.cpp:53-111` | `vector<vector<double>>` intensity matrices -> `vector<vector<float>>` |
| `ANALYSIS/OPENSWATH/OpenSwathScores.h:190-231` | `OpenSwath_Ind_Scores`: ~31 `vector<double>` fields -> `vector<float>` |
| `ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.cpp:457-458` | Per-transition score vectors -> `vector<float>` |
| `ANALYSIS/OPENSWATH/OpenSwathScoring.cpp:454-487` | Intensity arrays feeding into MRMScoring -> `vector<float>` |
| `ANALYSIS/OPENSWATH/PeakPickerChromatogram.cpp:216` | `vector<double> intensity` -> `vector<float>` |
| `ANALYSIS/OPENSWATH/IonMobilityScoring.cpp:116-135` | `extractIntensities()` return/output -> `vector<float>` |
| `ANALYSIS/OPENSWATH/DIAHelper.h:63` | `integrated_windows_intensity` parameter -> `vector<float>` (m/z and IM stay double) |
| `PROCESSING/SMOOTHING/GaussFilter.cpp:50-142` | `int_in`, `int_out` -> `vector<float>` (coordinate vectors stay double) |

**Crawdad boundary**: `PeakPickerChromatogram` converts `vector<float>` to `vector<double>`
at the `CrawdadWrapper::SetChromatogram()` call site only. This is a one-time cost per
chromatogram and Crawdad is optional (`#ifdef WITH_CRAWDAD`).

**Cascading changes**: Functions calling these (e.g., `MRMFeatureFinderScoring` calling
`MRMScoring`, `OpenSwathScoring` populating `OpenSwath_Ind_Scores`) need signature alignment.
The `OpenSwathOSWWriter` and `OpenSwathOSWParquetWriter` read `Ind_Scores` fields via
`getSeparateScore()` -- these need to accept `vector<float>`.

**Testing**: OpenSWATH tests, GaussFilter tests, PeakPicker tests. May need minor tolerance
adjustments if reference values were computed at double precision.

### PR 3 -- Quantitation and metadata

| File | Change |
|------|--------|
| `ANALYSIS/QUANTITATION/IsobaricChannelExtractor.cpp:688` | `extractSingleSpec()` return type -> `vector<float>` |
| `ANALYSIS/QUANTITATION/IsobaricIsotopeCorrector.cpp:29-68` | Buffer vectors -> `vector<float>`; convert to `double` at NNLS solver call boundary only |
| `ANALYSIS/MRM/ReactionMonitoringTransition.h:292` | `double library_intensity_` -> `Peak1D::IntensityType` |
| `ANALYSIS/TARGETED/MetaboTargetedAssay.h:33` | `double precursor_int` -> `Peak1D::IntensityType` |

**NNLS boundary**: The `NonNegativeLeastSquaresSolver::solve()` requires `Matrix<double>`
and `vector<double>`. Convert `vector<float>` to `vector<double>` at the solver call site.
Channel arrays are small (~10 elements), so the conversion cost is negligible.

**Testing**: Isobaric quantitation tests, MRM transition tests. Verify serialization of
`library_intensity_` and `precursor_int` to files is unaffected.

### PR 4 -- I/O in-memory buffers

Change in-memory buffer types during file I/O. Wire format stays `double` on disk.

| File | Change |
|------|--------|
| `FORMAT/DATAACCESS/MSChromatogramParquetConsumer.cpp:419-430` | `vector<double> intensity_data` -> `vector<float>`; use `encodeNP()` (float overload) instead of `encodeNPRaw()` |
| `FORMAT/XICParquetFile.h:91-92` | `vector<double> intensity` -> `vector<float>` in output struct |
| `FORMAT/XICParquetFile.cpp:330-351` | Decode as double from disk, narrow to float for output |
| `FORMAT/HANDLERS/CachedMzMLHandler.h:52-54` | `DatumSingleton` stays `double` (wire format type); convert to float after read, widen to double before write |

**MSNumpressCoder**: Has a `vector<float>` overload for `encodeNP()` that converts internally
to double before encoding. The `encodeNPRaw()` variant only accepts `vector<double>`. For
`MSChromatogramParquetConsumer`, switch from `encodeNPRaw()` to `encodeNP()` to use the float
overload, or convert at the call site.

**Testing**: Roundtrip tests -- write chromatograms, read back, verify intensities match at
float precision.

## Files Changed (complete list)

### PR 0 (openswathalgo)
- `src/openswathalgo/include/OpenMS/OPENSWATHALGO/ALGO/Scoring.h`
- `src/openswathalgo/source/OPENSWATHALGO/ALGO/Scoring.cpp`
- `src/openswathalgo/include/OpenMS/OPENSWATHALGO/ALGO/StatsHelpers.h`
- `src/openswathalgo/source/OPENSWATHALGO/ALGO/StatsHelpers.cpp`

### PR 1 (cleanup)
- `src/openms/include/OpenMS/PROCESSING/RESAMPLING/LinearResampler.h`
- `src/openms/source/ANALYSIS/OPENSWATH/MRMScoring.cpp`
- `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h`
- `src/openms/source/ANALYSIS/QUANTITATION/IsotopeLabelingMDVs.cpp`
- `src/openms/source/FORMAT/QcMLFile.cpp`
- `src/openms/source/PROCESSING/SMOOTHING/GaussFilter.cpp`

### PR 2 (algorithm internals)
- `src/openms/source/ANALYSIS/OPENSWATH/MRMScoring.cpp`
- `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/OpenSwathScores.h`
- `src/openms/source/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.cpp`
- `src/openms/source/ANALYSIS/OPENSWATH/OpenSwathScoring.cpp`
- `src/openms/source/ANALYSIS/OPENSWATH/PeakPickerChromatogram.cpp`
- `src/openms/source/ANALYSIS/OPENSWATH/IonMobilityScoring.cpp`
- `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/DIAHelper.h`
- `src/openms/source/ANALYSIS/OPENSWATH/DIAHelper.cpp`
- `src/openms/source/ANALYSIS/OPENSWATH/DIAPrescoring.cpp`
- `src/openms/source/PROCESSING/SMOOTHING/GaussFilter.cpp`

### PR 3 (quantitation and metadata)
- `src/openms/source/ANALYSIS/QUANTITATION/IsobaricChannelExtractor.cpp`
- `src/openms/include/OpenMS/ANALYSIS/QUANTITATION/IsobaricChannelExtractor.h`
- `src/openms/source/ANALYSIS/QUANTITATION/IsobaricIsotopeCorrector.cpp`
- `src/openms/include/OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h`
- `src/openms/source/ANALYSIS/MRM/ReactionMonitoringTransition.cpp`
- `src/openms/include/OpenMS/ANALYSIS/TARGETED/MetaboTargetedAssay.h`

### PR 4 (I/O in-memory buffers)
- `src/openms/source/FORMAT/DATAACCESS/MSChromatogramParquetConsumer.cpp`
- `src/openms/include/OpenMS/FORMAT/XICParquetFile.h`
- `src/openms/source/FORMAT/XICParquetFile.cpp`
- `src/openms/include/OpenMS/FORMAT/HANDLERS/CachedMzMLHandler.h`
- `src/openms/source/FORMAT/HANDLERS/CachedMzMLHandler.cpp`

## Risks and Mitigations

| Risk | Mitigation |
|------|------------|
| Float accumulation precision loss in large chromatograms | Monitor test results; add pairwise summation as targeted fix if needed |
| Test reference value mismatches | Expected for values computed at double precision; update tolerances where the float result is equally correct |
| OpenSwath_Ind_Scores field type change breaks downstream consumers | Trace all readers (OSWWriter, OSWParquetWriter, pyprophet); update in same PR |
| Eigen float path performance differs from double | Benchmark cross-correlation on representative data in PR 0 |
| CachedMzMLHandler read/write asymmetry during rollout | Old files (double on disk) must still read correctly; test with existing cached files |

## Out of Scope

- Changing `BinaryDataArray::data` to `vector<float>` (deferred, tracked in issue #8872)
- Changing binary wire formats on disk
- Templating the NNLS solver
- Templating CrawdadWrapper
- Changing `OpenSwath_Scores` scalar fields (42 doubles -- computed results, not arrays)
- pyOpenMS binding changes (float IntensityType is already handled by PR #8857)
