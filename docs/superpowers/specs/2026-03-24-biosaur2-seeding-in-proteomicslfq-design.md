# Biosaur2 as User-Selectable Seed Generator in ProteomicsLFQ

## Summary

Add a `Seeding:algorithm` parameter to ProteomicsLFQ that lets users choose between the current `FeatureFinderMultiplexAlgorithm`-based seeding (`"multiplex"`, default) and `Biosaur2Algorithm` (`"biosaur2"`) for untargeted seed generation. The rest of the pipeline (FeatureFinderIdentificationAlgorithm, alignment, linking, quantification) remains unchanged.

## Scope

**In scope:** User-selectable seed algorithm in ProteomicsLFQ only.

**Out of scope:** Replacing FeatureFinderIdentificationAlgorithm, IDMapper-based workflows, changes to DDAWorkflowCommons, or Biosaur2 integration in other workflow tools.

## Design

### Parameter Registration

In `registerOptionsAndFlags_()`:

1. Add `Seeding:algorithm` — string parameter, valid values `"multiplex"` (default) and `"biosaur2"`. Controls which algorithm generates untargeted feature seeds.

2. Add `Seeding:Biosaur2:` subsection — populated from `Biosaur2Algorithm().getDefaults()`. Add to the existing `combined` Param object (alongside Centroiding, PeptideQuantification, etc.) with prefix `"Seeding:Biosaur2:"` before the single `registerFullParam_(combined)` call. Mark all Biosaur2 parameters as advanced by iterating over the defaults and calling `addTag(key, "advanced")` on each, matching the existing pattern used for PeakPickerHiRes parameters (lines 247-251).

3. Existing `Seeding:intThreshold`, `Seeding:charge`, `Seeding:traceRTTolerance` are kept for backwards compatibility. Note: `Seeding:charge` and `Seeding:traceRTTolerance` are currently dead code — they are registered but never read. The `calculateSeeds()` call hardcodes charge range `2:5`. This is a pre-existing issue, not introduced by this change.

### Seed Generation Branch

In `quantifyFraction_()` (~line 927-934), replace the unconditional `DDAWorkflowCommons::calculateSeeds()` call with a branch:

- **`"multiplex"`** (default): Call `DDAWorkflowCommons::calculateSeeds()` as today.
- **`"biosaur2"`**: Instantiate `Biosaur2Algorithm`, set parameters from `Seeding:Biosaur2:` subsection, pass `ms_centroided` via the **copy** overload of `setMSData()`, call `run(seeds)`.

**Why copy, not move:** `Biosaur2Algorithm::run()` destructively consumes its internal `MSExperiment` — it erases non-MS1 spectra and moves the data into FAIMS groups (`std::move(ms_data_)` at line 228). After `run()`, `getMSData()` returns a moved-from, empty experiment. Since `ms_centroided` is still needed downstream by `FeatureFinderIdentificationAlgorithm`, it must be preserved. The copy overload `setMSData(ms_centroided)` is the simplest solution. This doubles peak memory for the MS1 data during seeding — acceptable since `calculateSeeds()` also creates a filtered copy internally.

**Logging:** The biosaur2 branch should log which algorithm was selected and the number of seeds produced, matching the existing log output pattern from `calculateSeeds()`.

### Caveat: Redundant Processing

Biosaur2Algorithm internally filters to MS1 spectra and supports its own centroiding via `profile_mode`. Since ProteomicsLFQ already provides centroided MS1 data, the MS1 filtering is redundant (harmless) and `profile_mode` already defaults to `false` in Biosaur2Algorithm — no explicit override needed.

### FAIMS Interaction

Biosaur2Algorithm has its own FAIMS handling (splits by compensation voltage, processes groups, merges features). When used as a seeder, this is acceptable — biosaur2 produces a merged `FeatureMap` of seeds across all CV groups. `FeatureFinderIdentificationAlgorithm` then does its own FAIMS-aware processing on the full experiment. Seeds are matched by RT/m/z/charge regardless of FAIMS origin, so there is no conflict. However, if FAIMS-specific seed filtering is ever needed, this interaction should be revisited.

### Downstream Impact

None. The seeds `FeatureMap` is passed to `FeatureFinderIdentificationAlgorithm::run()` exactly as before — it only uses RT, m/z, charge, and intensity from seeds, all of which Biosaur2Algorithm provides. Alignment, linking, quantification, and SVM-based quality filtering are unaffected.

## Testing

- Existing ProteomicsLFQ tests pass unchanged (default `algorithm = "multiplex"` preserves current behavior).
- New or extended test with `Seeding:algorithm = "biosaur2"` to verify seeds are generated and the pipeline completes. Existing test data can be reused since Biosaur2Algorithm works on any centroided MS1 input.

## Files Changed

- `src/topp/ProteomicsLFQ.cpp` — parameter registration and seed generation branch (only file with functional changes)
