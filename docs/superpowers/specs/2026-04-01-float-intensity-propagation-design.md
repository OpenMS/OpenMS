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
- **SIMD-friendly**: Where Eigen is used (cross-correlation, standardization), `VectorXf`
  processes 2x elements per SIMD instruction vs `VectorXd`. Note: many non-Eigen loops
  (push_back, branched filter kernels) won't auto-vectorize regardless of type — the primary
  win in those paths is memory reduction, not throughput.
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

All-float, including loop accumulators and running sums. For Eigen-based scoring paths,
uniform float avoids `vcvtps2pd` conversion instructions and enables 2x SIMD throughput.
For non-Eigen loops, the benefit is primarily memory reduction (50% for intensity arrays).

For typical OpenSWATH data sizes (hundreds to low thousands of points per chromatogram),
float accumulation precision is adequate. If a specific algorithm shows precision issues,
pairwise summation can be added as a targeted fix.

Note: `GaussFilterAlgorithm::integrate_()` returns `double` internally. The `int_out`
assignment requires an explicit `static_cast<float>` to avoid narrowing warnings.
`mean_and_stddev` in StatsHelpers.h also uses hardcoded `double` internally — this stays
double as it's a statistical computation functor, not intensity storage.

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
| `IFeature::getIntensity(vector<double>&)` | Virtual interface in `openswathalgo`; `FeatureOpenMS` materializes `ConvexHull2D` points as `double`. Add `getIntensity(vector<float>&)` as non-pure-virtual with default impl delegating to double overload + conversion. `MRMScoring::fillIntensityFromFeature()` then calls the float overload. This avoids an ABI break while enabling float intensity matrices |
| `PeakIntegrator` hull points | `PeakIntegrator.h` materializes convex hull points as `DPosition<2>` (double). Stays double — these feed into `FeatureOpenMS::getIntensity()` which handles conversion |
| `OpenSwath_Scores` scalar fields (~42 doubles) | Computed score results, not intensity storage arrays; no memory benefit from changing scalars |
| `DataValue::DOUBLE_LIST` | `MRMFeature.cpp` stores `OpenSwath_Ind_Scores` vectors as `DataValue` meta values typed `DOUBLE_LIST`. Convert `vector<float>` to `DoubleList` at the `MRMFeature` boundary. OSW writers read these via `Feature::getMetaValue()` as `DOUBLE_LIST` — no writer changes needed |

## PR Structure

Five PRs merged in sequence, ordered by risk tier.

### PR 0 -- Template OpenSwathAlgo for float

Template the scoring and stats helper functions in `openswathalgo` to accept both `float`
and `double`. This is the foundation PR that unblocks PR 2.

**Scoring.h / Scoring.cpp** -- 11 functions to template on scalar type `T`:

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
| `computeAndAppendRank(const vector<T>&, vector<unsigned int>&)` | Called by computeRankVector; used by MI scoring |
| `computeRankVector(const vector<vector<T>>&, vector<vector<unsigned int>>&)` | Calls computeAndAppendRank; used in MRMScoring and MRMTransitionGroupPicker |

Return type stays `XCorrArrayType` (`vector<pair<int, double>>`) -- correlation coefficients
remain double precision.

**IFeature interface** -- Add `getIntensity(vector<float>&)` overload:

Add a non-pure-virtual method `getIntensity(std::vector<float>& intens)` to `IFeature`
with a default implementation that calls the existing `getIntensity(vector<double>&)` and
narrows. Override in `FeatureOpenMS` to provide a direct float extraction. This allows
`MRMScoring::fillIntensityFromFeature()` to populate `vector<vector<float>>` matrices
without an intermediate double copy.

Files: `ITransition.h` (interface), `MRMFeatureAccessOpenMS.h/.cpp` (FeatureOpenMS impl),
`MockObjects.h` (test mock — add float overload).

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
instantiations. Add float-specific test cases to:

- `Scoring_test.cpp` — covers cross-correlation and rank helpers; add float variants
- `DiaHelpers_test.cpp` — constructs `vector<double>` for `dotprodScoring`/`manhattanScoring`; add float coverage
- `StatsHelpers_test.cpp` (in openms tests) — constructs `vector<double>` intensity vectors; add float coverage

### PR 1 -- Cleanup (zero behavioral change)

Remove redundant casts and unnecessary `double` promotions where `getIntensity()` already
returns `float`.

| File | Change |
|------|--------|
| `PROCESSING/RESAMPLING/LinearResampler.h:106-111` | Remove four `static_cast<double>(getIntensity())` calls in the resampling block |
| `ANALYSIS/OPENSWATH/MRMScoring.cpp:556` | Remove `static_cast<double>(getIntensity())` |
| `ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h:491` | `const double` -> `const auto` / `IntensityType` |
| `QUANTITATION/IsotopeLabelingMDVs.cpp:204` | Remove `(double)(it->getIntensity())` C-style cast |
| `FORMAT/QcMLFile.cpp:1275` | `double sum` -> `float` / `IntensityType` |
| `PROCESSING/SMOOTHING/GaussFilter.cpp:62` | Remove `static_cast<double>` (prepares for PR 2) |

**Testing**: All existing unit tests pass unchanged.

### PR 2 -- Algorithm internals (main memory savings)

Change intensity storage vectors from `vector<double>` to `vector<float>` (via typedef where
appropriate) in scoring, peak picking, and signal processing.

| File | Change |
|------|--------|
| `ANALYSIS/OPENSWATH/MRMScoring.cpp:53-111` | `vector<vector<double>>` intensity matrices -> `vector<vector<float>>` |
| `ANALYSIS/OPENSWATH/OpenSwathScores.h:190-231` | `OpenSwath_Ind_Scores`: change only intensity-derived fields to `vector<float>` (see field audit below). Coordinate/score fields stay `double`. |
| `ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.cpp:457-458` | Per-transition score vectors -> `vector<float>` |
| `ANALYSIS/OPENSWATH/OpenSwathScoring.cpp:109,364,519` | Intensity arrays feeding into MRMScoring -> `vector<float>` |
| `ANALYSIS/OPENSWATH/PeakPickerChromatogram.cpp:216` | `vector<double> intensity` -> `vector<float>` |
| `ANALYSIS/OPENSWATH/IonMobilityScoring.cpp:116-135` | `extractIntensities()` return/output -> `vector<float>` |
| `ANALYSIS/OPENSWATH/DIAHelper.h:63` | `integrated_windows_intensity` parameter -> `vector<float>` (m/z and IM stay double) |
| `ANALYSIS/OPENSWATH/DIAPrescoring.cpp:163` | Update `intExp` variable to `vector<float>` to match changed `integrateWindows` signature |
| `PROCESSING/SMOOTHING/GaussFilter.cpp:50-142` | `int_in`, `int_out` -> `vector<float>` (coordinate vectors stay double) |

**Crawdad boundary**: `PeakPickerChromatogram` converts `vector<float>` to `vector<double>`
at the `CrawdadWrapper::SetChromatogram()` call site only. This is a one-time cost per
chromatogram and Crawdad is optional (`#ifdef WITH_CRAWDAD`).

**Cascading changes**: Functions calling these (e.g., `MRMFeatureFinderScoring` calling
`MRMScoring`, `OpenSwathScoring` populating `OpenSwath_Ind_Scores`) need signature alignment.
Header files with changed function signatures must also be updated (see file list below).

**DataValue boundary**: `MRMFeature.cpp` stores `OpenSwath_Ind_Scores` vectors as `DataValue`
meta values typed `DOUBLE_LIST`. Add explicit `vector<float>` -> `DoubleList` conversion in
`MRMFeature.cpp` at the storage boundary. The OSW writers (`OpenSwathOSWWriter`,
`OpenSwathOSWParquetWriter`) read these via `Feature::getMetaValue()` as `DOUBLE_LIST` and
need no changes themselves.

**Testing**: Test files requiring type or tolerance updates:

- `DIAHelper_test.cpp` — `vector<double> intInt` must become `vector<float>` (m/z and IM stay double)
- `MRMFeatureFinderScoring_test.cpp` — key regression test for the `MRMFeature.cpp` DataValue boundary
- `PeakPickerChromatogram_test.cpp`, `GaussFilter_test.cpp`, `IonMobilityScoring_test.cpp` — tolerance review
- `MRMScoring_test.cpp` — intensity fixtures flow through `MockObjects.h` which stays double (BinaryDataArray boundary); likely no type changes but verify tolerances
- `OpenSwathScoring_test.cpp` — facade tests over BinaryDataArray; no type edits expected

**OpenSwath_Ind_Scores field audit:**

| Category | Fields | Type |
|----------|--------|------|
| **Intensity** (change to `float`) | `ind_area_intensity`, `ind_total_area_intensity`, `ind_apex_intensity`, `ind_log_intensity`, `ind_intensity_score`, `ind_intensity_ratio`, `ind_im_log_intensity` | `vector<float>` |
| **Scores/coefficients** (stay `double`) | `ind_xcorr_coelution_score`, `ind_xcorr_shape_score`, `ind_log_sn_score`, `ind_isotope_correlation`, `ind_isotope_overlap`, `ind_massdev_score`, `ind_mi_score`, `ind_mi_ratio`, `ind_total_mi`, `ind_im_delta_score`, `ind_im_contrast_*`, `ind_im_sum_contrast_*` | `vector<double>` |
| **Coordinates/positions** (stay `double`) | `ind_apex_position`, `ind_start_position_at_5..50`, `ind_end_position_at_5..50`, `ind_im_drift`, `ind_im_drift_left`, `ind_im_drift_right`, `ind_im_delta` | `vector<double>` |
| **Shape metrics** (stay `double`) | `ind_tailing_factor`, `ind_asymmetry_factor`, `ind_slope_of_baseline`, `ind_baseline_delta_2_height`, `ind_fwhm` | `vector<double>` |

Only ~7 of ~40 fields are genuine intensity values. The rest are computed scores, RT/IM
coordinates, or shape metrics that should remain `double` per the type policy.

**`normalized_library_intensity` handling:** This parameter in `OpenSwathScoring` and
`MRMScoring` stays `vector<double>` (deferred). The templated scoring functions will be
instantiated for both float (intensity data) and double (library intensity) parameter types
at the respective call sites.

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
| `FORMAT/DATAACCESS/MSChromatogramParquetConsumer.cpp:419-430` | Keep `intensity_data` as `vector<double>` — the `encodeNP(vector<float>)` overload allocates a temporary `vector<double>` internally, so switching to float would add overhead vs the current `encodeNPRaw()` path. No change needed here. |
| `FORMAT/XICParquetFile.h:91-92` | `vector<double> intensity` -> `vector<float>` in output struct |
| `FORMAT/XICParquetFile.cpp:330-351` | Decode as double from disk, narrow to float for output |
| `FORMAT/HANDLERS/CachedMzMLHandler.h:52-54` | `DatumSingleton` stays `double` (wire format type); convert to float after read, widen to double before write |

**MSNumpressCoder**: The `encodeNP(vector<float>)` overload allocates a temporary
`vector<double>` internally, making it strictly slower than the current `encodeNPRaw()`
path. `MSChromatogramParquetConsumer` therefore keeps its intensity buffer as `vector<double>`
— the buffer is transient (per-chromatogram) and the Numpress encoding requires double
precision anyway.

**Testing**: Roundtrip tests -- write chromatograms, read back, verify intensities match at
float precision.

## Files Changed (complete list)

### PR 0 (openswathalgo + IFeature interface)
- `src/openswathalgo/include/OpenMS/OPENSWATHALGO/ALGO/Scoring.h`
- `src/openswathalgo/source/OPENSWATHALGO/ALGO/Scoring.cpp`
- `src/openswathalgo/include/OpenMS/OPENSWATHALGO/ALGO/StatsHelpers.h`
- `src/openswathalgo/source/OPENSWATHALGO/ALGO/StatsHelpers.cpp`
- `src/openswathalgo/include/OpenMS/OPENSWATHALGO/DATAACCESS/ITransition.h`
- `src/openswathalgo/include/OpenMS/OPENSWATHALGO/DATAACCESS/MockObjects.h`
- `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/MRMFeatureAccessOpenMS.h`
- `src/openms/source/ANALYSIS/OPENSWATH/DATAACCESS/MRMFeatureAccessOpenMS.cpp`

### PR 1 (cleanup)
- `src/openms/include/OpenMS/PROCESSING/RESAMPLING/LinearResampler.h`
- `src/openms/source/ANALYSIS/OPENSWATH/MRMScoring.cpp`
- `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h`
- `src/openms/source/ANALYSIS/QUANTITATION/IsotopeLabelingMDVs.cpp`
- `src/openms/source/FORMAT/QcMLFile.cpp`
- `src/openms/source/PROCESSING/SMOOTHING/GaussFilter.cpp`

### PR 2 (algorithm internals)
- `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MRMScoring.h`
- `src/openms/source/ANALYSIS/OPENSWATH/MRMScoring.cpp`
- `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/OpenSwathScores.h`
- `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/OpenSwathScoring.h`
- `src/openms/source/ANALYSIS/OPENSWATH/OpenSwathScoring.cpp`
- `src/openms/source/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.cpp`
- `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/PeakPickerChromatogram.h`
- `src/openms/source/ANALYSIS/OPENSWATH/PeakPickerChromatogram.cpp`
- `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/IonMobilityScoring.h`
- `src/openms/source/ANALYSIS/OPENSWATH/IonMobilityScoring.cpp`
- `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/DIAHelper.h`
- `src/openms/source/ANALYSIS/OPENSWATH/DIAHelper.cpp`
- `src/openms/source/ANALYSIS/OPENSWATH/DIAPrescoring.cpp`
- `src/openms/source/KERNEL/MRMFeature.cpp`
- `src/openms/source/PROCESSING/SMOOTHING/GaussFilter.cpp`

### PR 3 (quantitation and metadata)
- `src/openms/source/ANALYSIS/QUANTITATION/IsobaricChannelExtractor.cpp`
- `src/openms/include/OpenMS/ANALYSIS/QUANTITATION/IsobaricChannelExtractor.h`
- `src/openms/source/ANALYSIS/QUANTITATION/IsobaricIsotopeCorrector.cpp`
- `src/openms/include/OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h`
- `src/openms/source/ANALYSIS/MRM/ReactionMonitoringTransition.cpp`
- `src/openms/include/OpenMS/ANALYSIS/TARGETED/MetaboTargetedAssay.h`

### PR 4 (I/O in-memory buffers)
- `src/openms/include/OpenMS/FORMAT/XICParquetFile.h`
- `src/openms/source/FORMAT/XICParquetFile.cpp`
- `src/openms/include/OpenMS/FORMAT/HANDLERS/CachedMzMLHandler.h`
- `src/openms/source/FORMAT/HANDLERS/CachedMzMLHandler.cpp`

## Risks and Mitigations

| Risk | Mitigation |
|------|------------|
| Float accumulation precision loss in large chromatograms | Monitor test results; add pairwise summation as targeted fix if needed |
| Test reference value mismatches | Expected for values computed at double precision; update tolerances where the float result is equally correct |
| OpenSwath_Ind_Scores field type change breaks downstream consumers | Convert `vector<float>` to `DoubleList` at `MRMFeature.cpp` boundary; OSW writers read via `DataValue` and need no changes |
| Eigen float path performance differs from double | Benchmark cross-correlation on representative data in PR 0 |
| CachedMzMLHandler read/write asymmetry during rollout | Old files (double on disk) must still read correctly; test with existing cached files |
| New double->float conversion at `IFeature::getIntensity()` boundary | Default impl converts from double; `FeatureOpenMS` override can extract float directly from hull points to avoid copy. Benchmark to confirm no regression |
| `GaussFilterAlgorithm::integrate_()` double->float narrowing | Add explicit `static_cast<float>` at assignment to suppress compiler warnings |

## Conversion Map

Data flow through the main OpenSWATH chromatogram scoring path, showing where float/double
boundaries exist after this work:

```
MSChromatogram<float>
  -> BinaryDataArray<double>         [stays double — interface boundary]
    -> DataAccessHelper               [double -> float conversion]
      -> MSChromatogram<float>
        -> GaussFilter<float>          [PR 2: was double, now float internally]
          -> MSChromatogram<float>
            -> PeakIntegrator          [stays double — hull points]
              -> FeatureOpenMS
                -> getIntensity(vector<float>&)  [PR 0: new float overload]
                  -> MRMScoring<float>             [PR 2: was double]
                    -> Scoring::normalizedCrossCorrelationPost<float>  [PR 0: templated]
                      -> XCorrArrayType<double>    [stays double — correlation coefficients]
                        -> OpenSwath_Scores<double>  [stays double — scalar results]

OpenSwath_Ind_Scores (7 float fields, ~33 double fields)
  -> MRMFeature::IDScoresAsMetaValue()
    -> DataValue::DOUBLE_LIST          [stays double — conversion at boundary]
      -> OSW writers                     [unchanged]
```

New conversions introduced by this work:
- `FeatureOpenMS::getIntensity(vector<float>&)`: double hull points -> float (one-time per feature)
- `MRMFeature.cpp`: `vector<float>` -> `DoubleList` for Ind_Scores storage
- `CrawdadWrapper`: `vector<float>` -> `vector<double>` (optional path only)

Conversions removed:
- `GaussFilter.cpp`: eliminated float->double->float round-trip
- Scoring functions: eliminated double intermediates for intensity computation

## Out of Scope

- Changing `BinaryDataArray::data` to `vector<float>` (deferred, tracked in issue #8872)
- Changing binary wire formats on disk
- Templating the NNLS solver
- Templating CrawdadWrapper
- Changing `OpenSwath_Scores` scalar fields (42 doubles -- computed results, not arrays)
- Mobilogram/XIM analogous paths (`PeakPickerMobilogram.h`, `MobilogramParquetConsumer.cpp`, `XIMParquetFile.h`, `PeakPickerIM.cpp`) -- similar patterns exist but are deferred to a follow-up
- Additional intensity-double sites in `ConfidenceScoring.cpp`, `ModifiedSincSmoother.cpp`, `OpenSwathOSWParquetReader.h` -- deferred to follow-up
- `normalized_library_intensity` parameter in `MRMScoring.h`, `OpenSwathScoring.h`, `DIAScoring.h` -- deferred
- pyOpenMS binding changes (float IntensityType is already handled by PR #8857)
