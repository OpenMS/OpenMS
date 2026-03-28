# Bruker .d (IM_PEAK) Support in ProteomicsLFQ

## Summary

Add Bruker `.d` file support to ProteomicsLFQ by introducing an IM_PEAK-aware code path that bypasses incompatible preprocessing steps (PeakPickerHiRes, PrecursorCorrection, MassTraceDetection-based FWHM, FeatureFinderMultiplex seeding) while preserving the existing FFId-based targeted extraction pipeline. Biosaur2Algorithm replaces FeatureFinderMultiplexAlgorithm for seed generation (user-selectable for mzML, mandatory for .d).

## Scope

**In scope:**
- Accept `.d` directories as input to ProteomicsLFQ
- IM_PEAK-aware preprocessing branch (skip PeakPickerHiRes, PrecursorCorrection, clearMetaDataArrays)
- Biosaur2Algorithm as user-selectable seed generator (`Seeding:algorithm` parameter)
- FWHM estimation from Biosaur2 features for .d input
- Preserve IM meta values on features through the pipeline

**Out of scope:**
- Replacing FeatureFinderIdentificationAlgorithm (it already handles IM_PEAK correctly)
- IM-aware feature grouping or alignment (existing RT x m/z matching unchanged)
- IDMapper-based workflows
- Changes to DDAWorkflowCommons
- Other TOPP tools

## Prerequisites

- SageAdapter produces idXML with `"IM"` meta value on PeptideIdentification (already working)
- `FileHandler` can load `.d` via `BrukerTimsFile` (already working)
- `FeatureFinderIdentificationAlgorithm` supports IM_PEAK format with 2D chromatogram extraction (already working)
- `Biosaur2Algorithm` handles IM_PEAK/PASEF data natively (already working)

## Design

### 1. Input Format Registration

In `registerOptionsAndFlags_()`:

```cpp
setValidFormats_("in", ListUtils::create<String>("mzML,d"));
```

No new user-facing parameters for .d detection — file type is determined automatically by `FileHandler`.

### 2. Format-Aware Loading and Preprocessing

Replace the direct call to `centroidAndCorrectPrecursors_()` in `quantifyFraction_()` with a new method `loadAndPreprocess_()` that branches on file type.

**Method signature:**
```cpp
ExitCodes loadAndPreprocess_(const String& mz_file, MSExperiment& ms_out, bool& is_im_peak_data)
```

**mzML path (unchanged behavior):**
1. `FileHandler().loadExperiment(mz_file, ms_raw, {FileTypes::MZML})`
2. `ms_raw.clearMetaDataArrays()` — safe for non-IM data, frees memory
3. `PeakPickerHiRes::pickExperiment()` — centroid MS1
4. `PrecursorCorrection::correctToHighestIntensityMS1Peak()` — fix precursor m/z
5. `is_im_peak_data = false`

**.d path (new):**
1. `FileHandler().loadExperiment(mz_file, ms_out, {FileTypes::BRUKER_TDF})` — loads with IM float arrays preserved
2. **No** `clearMetaDataArrays()` — IM per-peak arrays must survive
3. Remove MS2 peak data, sort by position (same housekeeping as mzML path)
4. **No** PeakPickerHiRes — Bruker TOF data is sparse/centroid-like; IM arrays would be destroyed
5. **No** PrecursorCorrection — `findHighestInWindow()` is IM-unaware, would pick peaks from wrong IM slice
6. `is_im_peak_data = true`

**Rationale for skipping PrecursorCorrection on .d data:** With IM_PEAK MS1 frames, each "spectrum" contains peaks across the full ion mobility range. `findHighestInWindow()` selects the highest-intensity peak within an m/z window regardless of IM value, which can match a precursor to a completely different ion species at a different mobility. Sage already handles precursor mass assignment during search, so correction is not critical.

### 3. Seed Generation — Biosaur2 as Alternative to Multiplex

Implements the existing design spec from `2026-03-24-biosaur2-seeding-in-proteomicslfq-design.md` with one addition: for .d/IM_PEAK input, Biosaur2 is forced (Multiplex cannot handle IM_PEAK).

**Parameter registration:**

```cpp
registerStringOption_("Seeding:algorithm", "<choice>", "multiplex",
  "Algorithm for untargeted seed feature detection.\n"
  "multiplex: FeatureFinderMultiplexAlgorithm (default, current behavior).\n"
  "biosaur2: Biosaur2Algorithm (handles IM_PEAK/PASEF data natively).",
  false, false);
setValidStrings_("Seeding:algorithm", {"multiplex", "biosaur2"});
```

Biosaur2 parameter subsection added to the combined Param:
```cpp
Param bio_defaults = Biosaur2Algorithm().getDefaults();
for (auto it = bio_defaults.begin(); it != bio_defaults.end(); ++it)
{
  bio_defaults.addTag(it.getName(), "advanced");
}
combined.insert("Seeding:Biosaur2:", bio_defaults);
```

**Seed generation branch in `quantifyFraction_()`:**

```
if is_im_peak_data:
  force seeding_algorithm = "biosaur2" (warn if user set "multiplex")

if seeding_algorithm == "biosaur2":
  Biosaur2Algorithm bio;
  bio.setParameters(getParam_().copy("Seeding:Biosaur2:", true));
  bio.setMSData(ms_centroided);  // copy overload — run() consumes internal data
  bio.run(seeds);
else:
  DDAWorkflowCommons::calculateSeeds(...)  // existing path
```

**Key details:**
- `setMSData()` uses the copy overload because `run()` destructively consumes internal MSExperiment data. `ms_centroided` is still needed downstream by FFId.
- Biosaur2 auto-detects IM_PEAK format and handles PASEF centroiding via `paseftol` parameter.
- `profile_mode` defaults to `false` — appropriate for Bruker TOF data.
- FAIMS splitting is handled internally by Biosaur2 (no conflict with FFId's own FAIMS handling).

### 4. FWHM Estimation for .d Input

**Problem:** `DDAWorkflowCommons::estimateMedianChromatographicFWHM()` uses `MassTraceDetection` which doesn't understand IM_PEAK format.

**Solution:** For .d input, estimate FWHM from Biosaur2 seed features. This requires reordering: seeds are generated first, then FWHM is derived from them.

**Current order (mzML):**
1. FWHM estimation via MassTraceDetection
2. Seed generation (uses FWHM for RT parameters)
3. FFId (uses FWHM for peak_width)

**New order (.d path):**
1. Seed generation via Biosaur2 (doesn't need external FWHM — has own peak width handling)
2. FWHM estimation from Biosaur2 features: use the `run(fm, hills, peptide_features)` overload; for each `PeptideFeature`, compute `fwhm = rt_end - rt_start`; take the median across all features
3. FFId uses that FWHM for `detect:peak_width`

**Fallback:** If Biosaur2 produces fewer than 10 features, use a hardcoded default (30 seconds) with a warning log message.

### 5. Feature Meta Value Preservation

**Current behavior (line 1125-1133):** Strips all feature meta values except `"OffsetPeptide"`.

**Change:** Extend the keep set to include IM-related meta values:

```cpp
unordered_set<String> keep_meta = {"OffsetPeptide", "IM_median", "IM_min", "IM_max"};
```

This is safe for mzML input — if those keys don't exist on a feature, nothing changes.

### 6. Downstream Pipeline (Unchanged)

All steps after feature detection are unaffected:

- **FAIMS annotation** (`addMissingFAIMSToPeptideIDs`): No-op for TIMS data (no FAIMS CVs). No change needed.
- **ID meta value survival:** `loadAndCleanupIDFile_()` strips PeptideHit meta values but explicitly preserves PeptideIdentification meta values (including `"IM"` from Sage). Confirmed safe.
- **FFId:** Auto-detects IM_PEAK, reads `"IM"` from peptide IDs, computes per-peptide IM statistics, performs 2D chromatogram extraction (m/z + IM windowing). No changes needed.
- **Alignment** (`MapAlignmentAlgorithmIdentification`): RT-based, uses identified peptide RTs. Unaffected by IM data.
- **Linking** (`FeatureGroupingAlgorithmQT`): RT x m/z matching. IM not used but IM meta values survive on features.
- **Normalization, protein inference, FDR, quantification, export:** Operate on FeatureMap/ConsensusMap intensities. Unaffected.

## Files Changed

- `src/topp/ProteomicsLFQ.cpp` — all functional changes in a single file:
  - Add `#include <OpenMS/FEATUREFINDER/Biosaur2Algorithm.h>`
  - Add `"d"` to valid input formats
  - Add `Seeding:algorithm` parameter and `Seeding:Biosaur2:` subsection
  - New `loadAndPreprocess_()` method branching on file type
  - Reordered .d path in `quantifyFraction_()`: seeds → FWHM → FFId
  - Seed generation branch (multiplex vs biosaur2)
  - Extended feature meta value keep set

## Testing

- Existing ProteomicsLFQ tests pass unchanged (default `Seeding:algorithm = "multiplex"` preserves current behavior for mzML input).
- New smoke test with `Seeding:algorithm = "biosaur2"` on existing mzML test data to verify the biosaur2 seeding path completes.
- If .d test data is available: end-to-end test with `.d` input + Sage idXML containing IM annotations.

## Known Limitations

- **Feature grouping ignores IM:** `FeatureGroupingAlgorithmQT` matches on RT x m/z only. Could mis-group isobaric species with different mobilities in complex samples. This is a pre-existing limitation, not introduced by this change.
- **Map alignment ignores IM:** RT-only alignment. Same as current FAIMS workflow.
- **PrecursorCorrection skipped for .d:** Precursor m/z may be slightly off for some spectra. Sage's own precursor handling mitigates this.
- **FWHM from Biosaur2 features is approximate:** Uses `PeptideFeature::rt_end - rt_start` rather than Gaussian fitting on mass traces. Sufficient for setting FFId's `detect:peak_width` parameter.
- **IDMapper not used:** IM dimension not considered during ID-to-feature mapping. FFId handles this internally via targeted extraction.
