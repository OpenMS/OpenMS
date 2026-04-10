# CometAdapter Bruker .d Support

**Date:** 2026-04-10
**Scope:** Add direct Bruker TimsTOF .d directory input to CometAdapter (DDA-PASEF)

## Context

CometAdapter currently only accepts mzML input. Comet itself (the external binary) only accepts mzML/mzXML — it cannot read .d natively. Several other OpenMS tools (FileConverter, PeakPickerIM, SimpleSearchEngine, SageAdapter) already support .d input via `BrukerTimsFile`. SageAdapter passes .d directly to Sage (which handles it natively), but CometAdapter must convert .d to a temporary indexed mzML before invoking Comet.

**Target use case:** DDA-PASEF data from Bruker TimsTOF instruments. DIA-PASEF is out of scope.

## Approach

In-memory load via `BrukerTimsFile` followed by writing a temporary indexed mzML for Comet. The loaded `MSExperiment` is reused for IM annotation in post-processing.

## Changes

### 1. CometAdapter.cpp — Input Registration

In `registerOptionsAndFlags_()`:

- Add `"d"` to `setValidFormats_("in", ...)` behind `#ifdef WITH_OPENTIMS`.
- Register a `bruker` TOPP subsection with parameters:
  - `bruker:export_mode` — `"auto"` default (auto-detects DDA vs DIA)
  - `bruker:calibration_tolerance` — m/z recalibration tolerance (0.0 default)
  - `bruker:calibrate` — enable m/z recalibration (false default)
- No MS1 centroiding parameters — only MS2 spectra are needed.

Add `getBrukerConfig_()` helper matching the pattern in PeakPickerIM/FileConverter.

### 2. CometAdapter.cpp — .d to Temp mzML Conversion

In `main_()`, after `getRawfileName()` and before Comet invocation:

- Call `FileHandler::getType(inputfile_name)` to detect `BRUKER_TDF`.
- If .d input:
  1. Build `BrukerTimsFile::Config` via `getBrukerConfig_()`.
  2. Load .d into `MSExperiment` via `BrukerTimsFile::load()` (full load with peak data, all MS levels — BrukerTimsFile has no MS-level filter).
  3. Remove non-MS2 spectra from the experiment (`std::erase_if` or equivalent) before writing.
  4. Write to a temporary mzML via `MzMLFile::store()` (produces indexed mzML, satisfying Comet's requirement).
  5. Set `input_file_with_index` to the temp mzML path.
  6. Skip the existing mzML index check block (not needed — we just wrote a fresh indexed file).
  7. Log debug-level native IDs from the loaded experiment for traceability.
- If mzML input: existing flow unchanged.

The `MSExperiment` from the .d load is kept for post-processing (IM annotation). Peak data can be cleared after writing the temp mzML to save memory, but spectrum metadata and native IDs are retained.

### 3. CometAdapter.cpp — Post-Processing

After Comet runs and pepXML is parsed:

- `setPrimaryMSRunPath` points to the original .d path (not the temp mzML).
- `reindex_()` — unchanged, works on native IDs regardless of source format.
- `addMissingIMToPeptideIDs()` — uses the `MSExperiment` from the .d load instead of reloading from mzML. For DDA-PASEF, MS2 spectra have `IM_SPECTRUM` format (scalar drift times), which this function handles. If `all_ids_have_im` is true, annotate the IM unit on the `ProteinIdentification`.
- `addMissingFAIMSToPeptideIDs()` — called unconditionally. Early-returns safely when no FAIMS data is present (line 238-243 in SpectrumMetaDataLookup.cpp checks `FAIMSHelper::getCompensationVoltages()` which returns empty for .d data).
- Native ID validation: log a warning if `addMissingIMToPeptideIDs` returns false for .d input (since DDA-PASEF data should always have IM).

### 4. SearchEngineBase.cpp — Centroid Warning Suppression

Add `case FileTypes::BRUKER_TDF: break;` in `getRawfileName()` switch statement. TimsTOF data is inherently centroided, so the generic "make sure spectra are centroided" warning is misleading. This is a one-liner.

Add `#include <OpenMS/FORMAT/FileTypes.h>` if needed for the `BRUKER_TDF` enum.

### 5. Compile Guards

All .d-related code is guarded by `#ifdef WITH_OPENTIMS` to match existing patterns. When built without opentims support, CometAdapter behaves exactly as before.

## Data Flow

```
.d directory
    |
    v
BrukerTimsFile::load() --> MSExperiment (MS2 only, with peak data)
    |                            |
    v                            | (retained for post-processing)
MzMLFile::store()                |
    |                            |
    v                            |
temp indexed mzML                |
    |                            |
    v                            |
Comet binary                     |
    |                            |
    v                            |
pepXML                           |
    |                            |
    v                            v
PepXMLFile::load()     addMissingIMToPeptideIDs()
    |                            |
    v                            v
PeptideIdentifications (with IM metavalues)
    |
    v
idXML output
```

## What Does NOT Change

- `reindex_()` — format-agnostic, works on native IDs
- `addMissingIMToPeptideIDs()` — works on any MSExperiment with IM_SPECTRUM spectra
- `addMissingFAIMSToPeptideIDs()` — safe no-op when no FAIMS data
- Comet invocation — always receives an indexed mzML
- pepXML parsing — format-agnostic
- All existing mzML inputs — identical behavior

## Testing

- Run CometAdapter with a DDA-PASEF .d directory; verify idXML output contains `IM` metavalues on peptide identifications.
- Verify `setPrimaryMSRunPath` points to the original .d path.
- Debug log confirms native ID round-trip consistency.
- Regression: existing mzML inputs produce identical results.
- Build without `WITH_OPENTIMS`: CometAdapter compiles and works as before (mzML only).
