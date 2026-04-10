# CometAdapter Bruker .d Support Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add direct Bruker TimsTOF .d directory input to CometAdapter, converting to temporary indexed mzML for Comet and preserving ion mobility annotations in the output idXML.

**Architecture:** Load .d via `BrukerTimsFile` into `MSExperiment`, filter to MS2, write temp indexed mzML for Comet, reuse experiment for IM annotation in post-processing. All .d code guarded by `#ifdef WITH_OPENTIMS`.

**Tech Stack:** C++, OpenMS (BrukerTimsFile, MzMLFile, FileHandler, SpectrumMetaDataLookup), CMake (WITH_OPENTIMS guard)

**Spec:** `docs/superpowers/specs/2026-04-10-cometadapter-bruker-d-support-design.md`

---

## File Map

- **Modify:** `src/topp/CometAdapter.cpp` — add .d input format, Bruker params, conversion logic, post-processing branch
- **Modify:** `src/openms/source/APPLICATIONS/SearchEngineBase.cpp` — add `BRUKER_TDF` case to `getRawfileName()`

---

### Task 1: Add BRUKER_TDF case to SearchEngineBase::getRawfileName

**Files:**
- Modify: `src/openms/source/APPLICATIONS/SearchEngineBase.cpp:36-82`

- [ ] **Step 1: Add the BRUKER_TDF case**

In `SearchEngineBase::getRawfileName()`, add a case for `BRUKER_TDF` before the `default` case. TimsTOF data is inherently centroided, so no warning is needed. The `FileHandler.h` include already provides `FileTypes`.

In `src/openms/source/APPLICATIONS/SearchEngineBase.cpp`, find the switch statement and add the new case:

```cpp
      case FileTypes::MGF:
        // no warning required. MGF files should be centroided by definition
        break;
      case FileTypes::BRUKER_TDF:
        // no warning required. TimsTOF data is inherently centroided
        break;
      default:
```

- [ ] **Step 2: Verify it compiles**

Run: `cmake --build OpenMS-build --target CometAdapter -j$(nproc) 2>&1 | tail -20`
Expected: Build succeeds (SearchEngineBase.cpp recompiles without errors)

- [ ] **Step 3: Commit**

```bash
git add src/openms/source/APPLICATIONS/SearchEngineBase.cpp
git commit -m "feat(SearchEngineBase): add BRUKER_TDF case to getRawfileName

TimsTOF data is inherently centroided, so suppress the misleading
'make sure spectra are centroided' warning for .d input."
```

---

### Task 2: Register .d input format and Bruker parameters in CometAdapter

**Files:**
- Modify: `src/topp/CometAdapter.cpp:9-38` (includes) and `src/topp/CometAdapter.cpp:110-253` (registerOptionsAndFlags_)

- [ ] **Step 1: Add BrukerTimsFile include**

At the top of `src/topp/CometAdapter.cpp`, after the existing includes (after line 34 `#include <OpenMS/ANALYSIS/ID/CometModification.h>`), add:

```cpp
#ifdef WITH_OPENTIMS
#include <OpenMS/FORMAT/BrukerTimsFile.h>
#endif
```

- [ ] **Step 2: Add .d to valid input formats**

In `registerOptionsAndFlags_()`, replace line 114:

```cpp
    setValidFormats_("in", { "mzML" } );
```

with:

```cpp
    setValidFormats_("in", { "mzML",
#ifdef WITH_OPENTIMS
      "d",
#endif
    });
```

- [ ] **Step 3: Register Bruker subsection parameters**

In `registerOptionsAndFlags_()`, just before the closing of the method (before line 253 which has `registerPeptideIndexingParameter_`), add:

```cpp
#ifdef WITH_OPENTIMS
    registerTOPPSubsection_("bruker", "Options for reading Bruker TimsTOF .d files (requires WITH_OPENTIMS)");
    registerStringOption_("bruker:export_mode", "<mode>", "auto", "Export mode: 'auto' detects DDA/DIA acquisition type, "
      "'spectrum' forces per-precursor spectra (DDA style).", false, true);
    setValidStrings_("bruker:export_mode", {"auto", "spectrum"});
    registerDoubleOption_("bruker:calibration_tolerance", "<float>", 0.0, "m/z recalibration tolerance (0 = library default)", false, true);
    setMinFloat_("bruker:calibration_tolerance", 0.0);
    registerStringOption_("bruker:calibrate", "<toggle>", "false", "Enable m/z recalibration (may fail on some datasets)", false, true);
    setValidStrings_("bruker:calibrate", {"true", "false"});
#endif
```

- [ ] **Step 4: Add getBrukerConfig_ helper**

In the `protected:` section of `TOPPCometAdapter` (after line 108 `map<string,int> num_enzyme_termini ...`), add:

```cpp
#ifdef WITH_OPENTIMS
  BrukerTimsFile::Config getBrukerConfig_()
  {
    BrukerTimsFile::Config c;
    c.calibration_tolerance = getDoubleOption_("bruker:calibration_tolerance");
    c.calibrate = (getStringOption_("bruker:calibrate") == "true");
    String mode = getStringOption_("bruker:export_mode");
    if (mode == "spectrum") c.export_mode = BrukerTimsFile::Config::SPECTRUM;
    else c.export_mode = BrukerTimsFile::Config::AUTO;
    return c;
  }
#endif
```

- [ ] **Step 5: Verify it compiles**

Run: `cmake --build OpenMS-build --target CometAdapter -j$(nproc) 2>&1 | tail -20`
Expected: Build succeeds

- [ ] **Step 6: Commit**

```bash
git add src/topp/CometAdapter.cpp
git commit -m "feat(CometAdapter): register .d input format and Bruker parameters

Add 'd' to valid input formats (guarded by WITH_OPENTIMS).
Register bruker:export_mode, bruker:calibration_tolerance, and
bruker:calibrate parameters matching the PeakPickerIM/FileConverter
pattern. Add getBrukerConfig_() helper."
```

---

### Task 3: Implement .d to temporary mzML conversion in main_()

**Files:**
- Modify: `src/topp/CometAdapter.cpp:602-684` (main_ method, input handling section)

This is the core change. The existing code (lines 660-684) does two things for mzML input:
1. Checks for mzML index and creates a temp indexed copy if missing (lines 660-676)
2. Loads metadata-only spectra for IM annotation (lines 678-684)

For .d input, we replace both with: load via BrukerTimsFile, filter to MS2, write temp mzML.

- [ ] **Step 1: Add file type detection and .d conversion branch**

In `main_()`, replace the block from line 660 (`// check for mzML index`) through line 684 (`mzml_file.load(inputfile_name, exp)`) with:

```cpp
    // Load input data — branch on file type
    MSExperiment exp;
    String input_file_with_index = inputfile_name;

#ifdef WITH_OPENTIMS
    const bool is_bruker_d = (FileHandler::getType(inputfile_name) == FileTypes::BRUKER_TDF);
    if (is_bruker_d)
    {
      // Load .d via BrukerTimsFile
      auto bruker_config = getBrukerConfig_();
      BrukerTimsFile tims_file;
      tims_file.setLogType(log_type_);
      tims_file.load(inputfile_name, exp, bruker_config);

      // Filter to target MS level only (Comet only needs MS2)
      std::erase_if(exp.getSpectra(), [&](const MSSpectrum& s) { return s.getMSLevel() != ms_level; });

      OPENMS_LOG_INFO << "Loaded " << exp.size() << " MS" << ms_level << " spectra from Bruker .d directory." << std::endl;

      if (!exp.empty())
      {
        writeDebug_("First native ID from .d: " + exp[0].getNativeID(), 2);
        writeDebug_("Last native ID from .d: " + exp.back().getNativeID(), 2);
      }

      // Write to temporary indexed mzML for Comet
      auto tmp_mzml = File::getTemporaryFile() + ".mzML";
      MzMLFile().store(tmp_mzml, exp);
      input_file_with_index = tmp_mzml;
    }
    else
#endif
    {
      // Existing mzML path
      MzMLFile mzml_file{};
      if (!mzml_file.hasIndex(inputfile_name))
      {
        OPENMS_LOG_WARN << "The mzML file provided to CometAdapter is not indexed, but comet requires one. "
                        << "We will add an index by writing a temporary file. If you run this analysis more often, consider indexing your mzML in advance!" << std::endl;
        auto tmp_file_mzml = File::getTemporaryFile() + ".mzML";
        PlainMSDataWritingConsumer consumer(tmp_file_mzml);
        consumer.getOptions().addMSLevel(ms_level);
        bool skip_full_count = true;
        mzml_file.transform(inputfile_name, &consumer, skip_full_count);
        input_file_with_index = tmp_file_mzml;
      }

      // Load spectra metadata to map to idXML
      mzml_file.getOptions().setMetadataOnly(false);
      mzml_file.getOptions().setFillData(false);
      mzml_file.getOptions().clearMSLevels();
      mzml_file.getOptions().addMSLevel(ms_level);
      mzml_file.load(inputfile_name, exp);
    }
```

Note: The variable declarations for `exp` and `input_file_with_index` that were previously at lines 661 and 663 are now at the top of this block. Remove the old declarations.

- [ ] **Step 2: Verify it compiles**

Run: `cmake --build OpenMS-build --target CometAdapter -j$(nproc) 2>&1 | tail -20`
Expected: Build succeeds

- [ ] **Step 3: Commit**

```bash
git add src/topp/CometAdapter.cpp
git commit -m "feat(CometAdapter): implement .d to temp mzML conversion

When input is a Bruker .d directory, load via BrukerTimsFile, filter
to MS2 spectra, and write a temporary indexed mzML for Comet. The
MSExperiment is retained for IM annotation in post-processing.
Existing mzML path is unchanged."
```

---

### Task 4: Update post-processing for .d input

**Files:**
- Modify: `src/topp/CometAdapter.cpp:722-749` (post-processing section)

Two changes needed:
1. `setPrimaryMSRunPath` must point to original .d path (already does — `inputfile_name` is still the .d path)
2. Add a warning if IM annotation fails for .d input (DDA-PASEF should always have IM)

- [ ] **Step 1: Add .d-specific IM validation warning**

The existing IM annotation code at lines 741-746 already works for both paths — `exp` contains the right spectra regardless of source. Add a warning for the .d case. Replace lines 741-746:

```cpp
    // Parse ion mobility information if present
    bool all_ids_have_im = SpectrumMetaDataLookup::addMissingIMToPeptideIDs(peptide_identifications, exp);
    if (all_ids_have_im)
    {
      protein_identifications[0].setMetaValue(Constants::UserParam::IM, exp.getSpectrum(0).getDriftTimeUnitAsString());
    }
```

with:

```cpp
    // Parse ion mobility information if present
    bool all_ids_have_im = SpectrumMetaDataLookup::addMissingIMToPeptideIDs(peptide_identifications, exp);
    if (all_ids_have_im)
    {
      protein_identifications[0].setMetaValue(Constants::UserParam::IM, exp.getSpectrum(0).getDriftTimeUnitAsString());
    }
#ifdef WITH_OPENTIMS
    else if (is_bruker_d)
    {
      OPENMS_LOG_WARN << "Warning: Bruker .d input but not all peptide IDs could be annotated with ion mobility values. "
                      << "This may indicate a native ID mismatch between the .d data and Comet results." << std::endl;
    }
#endif
```

- [ ] **Step 2: Verify `is_bruker_d` is in scope**

The `is_bruker_d` variable was declared inside `#ifdef WITH_OPENTIMS` in Task 3. The new code in Step 1 is also inside `#ifdef WITH_OPENTIMS`, so it compiles. However, `is_bruker_d` is currently declared inside a block scope (the if/else). We need to move it to function scope.

Check that in Task 3's code, `is_bruker_d` is declared at function level inside the `#ifdef WITH_OPENTIMS` guard (not inside a nested block). The code from Task 3 declares it as:
```cpp
#ifdef WITH_OPENTIMS
    const bool is_bruker_d = (FileHandler::getType(inputfile_name) == FileTypes::BRUKER_TDF);
```
This is at function scope within `main_()`, so it's accessible later. Good.

- [ ] **Step 3: Verify it compiles**

Run: `cmake --build OpenMS-build --target CometAdapter -j$(nproc) 2>&1 | tail -20`
Expected: Build succeeds

- [ ] **Step 4: Commit**

```bash
git add src/topp/CometAdapter.cpp
git commit -m "feat(CometAdapter): add IM validation warning for .d input

Log a warning if not all peptide IDs from a Bruker .d search could
be annotated with ion mobility values, indicating a potential native
ID mismatch."
```

---

### Task 5: Build verification and smoke test

**Files:**
- No file changes — verification only

- [ ] **Step 1: Full build of CometAdapter**

Run: `cmake --build OpenMS-build --target CometAdapter -j$(nproc) 2>&1 | tail -30`
Expected: Build succeeds with no warnings related to the changes

- [ ] **Step 2: Verify CometAdapter --help shows .d format and Bruker params**

Run: `OpenMS-build/bin/CometAdapter --help 2>&1 | grep -A2 -E "(Valid formats|bruker)"`
Expected: Output includes `d` in valid formats for `-in` and shows `bruker:export_mode`, `bruker:calibration_tolerance`, `bruker:calibrate` parameters.

If built without `WITH_OPENTIMS`, only `mzML` should appear and no `bruker:` parameters.

- [ ] **Step 3: Verify build without WITH_OPENTIMS (if feasible)**

Run: `cmake --build OpenMS-build --target CometAdapter -j$(nproc) -DWITH_OPENTIMS=OFF 2>&1 | tail -20`

Note: This may require a cmake reconfigure. If not easily testable, skip — the `#ifdef` guards follow the exact pattern used by PeakPickerIM and FileConverter, which are already tested in CI.

- [ ] **Step 4: Run existing CometAdapter tests (regression)**

Run: `ctest --test-dir OpenMS-build -R CometAdapter -V 2>&1 | tail -40`
Expected: All existing tests pass — mzML input path is unchanged.
