# ProteomicsLFQ Bruker .d Support Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add Bruker `.d` file support to ProteomicsLFQ with an IM_PEAK-aware preprocessing path and Biosaur2-based seed generation.

**Architecture:** Single-file change to `src/topp/ProteomicsLFQ.cpp`. Detect `.d` input early, branch into an IM_PEAK-aware loading path that skips PeakPickerHiRes/PrecursorCorrection/clearMetaDataArrays. Use Biosaur2Algorithm for seed generation (user-selectable for mzML, forced for `.d`). Derive FWHM from Biosaur2 features instead of MassTraceDetection. Rest of pipeline (FFId, alignment, linking, quant) unchanged.

**Tech Stack:** C++, OpenMS framework (TOPPBase, FileHandler, Biosaur2Algorithm, FeatureFinderIdentificationAlgorithm, IMTypes)

**Spec:** `docs/superpowers/specs/2026-03-28-proteomicslfq-dot-d-support-design.md`

---

### Task 1: Add Biosaur2Algorithm Include and Input Format Registration

**Files:**
- Modify: `src/topp/ProteomicsLFQ.cpp:54` (add include)
- Modify: `src/topp/ProteomicsLFQ.cpp:141` (add "d" to valid formats)

- [ ] **Step 1: Add the Biosaur2Algorithm and IMTypes includes**

After line 54 (`#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>`), add:

```cpp
#include <OpenMS/FEATUREFINDER/Biosaur2Algorithm.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
```

- [ ] **Step 2: Add "d" to valid input formats**

Change line 141 from:

```cpp
    setValidFormats_("in", ListUtils::create<String>("mzML"));
```

To:

```cpp
    setValidFormats_("in", ListUtils::create<String>("mzML,d"));
```

- [ ] **Step 3: Build to verify compilation**

Run: `cmake --build OpenMS-build --target ProteomicsLFQ -j$(nproc)`
Expected: Compiles successfully.

- [ ] **Step 4: Verify new format is accepted**

Run: `OpenMS-build/bin/ProteomicsLFQ --helphelp 2>&1 | grep -A1 '"in"'`
Expected: Shows `mzML,d` as valid formats.

- [ ] **Step 5: Commit**

```bash
git add src/topp/ProteomicsLFQ.cpp
git commit -m "feat(ProteomicsLFQ): add Biosaur2 include and .d input format"
```

---

### Task 2: Register Seeding:algorithm Parameter and Biosaur2 Subsection

**Files:**
- Modify: `src/topp/ProteomicsLFQ.cpp:239` (add Seeding:algorithm parameter)
- Modify: `src/topp/ProteomicsLFQ.cpp:293-301` (add Biosaur2 defaults to combined Param)

- [ ] **Step 1: Add Seeding:algorithm parameter**

After line 239 (`registerDoubleOption_("Seeding:traceRTTolerance", ...)`), add:

```cpp
    registerStringOption_("Seeding:algorithm", "<choice>", "multiplex",
      "Algorithm for untargeted seed feature detection.\n"
      "multiplex: FeatureFinderMultiplexAlgorithm (default, current behavior).\n"
      "biosaur2: Biosaur2Algorithm (handles IM_PEAK/PASEF data natively).",
      false, false);
    setValidStrings_("Seeding:algorithm", {"multiplex", "biosaur2"});
```

- [ ] **Step 2: Add Biosaur2 defaults to the combined Param**

Before line 294 (`Param combined;`), add:

```cpp
    Param bio_defaults = Biosaur2Algorithm().getDefaults();
    for (auto it = bio_defaults.begin(); it != bio_defaults.end(); ++it)
    {
      bio_defaults.addTag(it.getName(), "advanced");
    }
```

Then after line 299 (`combined.insert("ProteinQuantification:", pq_defaults);`), add:

```cpp
    combined.insert("Seeding:Biosaur2:", bio_defaults);
```

- [ ] **Step 3: Build and verify parameters**

Run: `cmake --build OpenMS-build --target ProteomicsLFQ -j$(nproc)`
Expected: Compiles successfully.

Run: `OpenMS-build/bin/ProteomicsLFQ --helphelp 2>&1 | grep -c "Seeding:Biosaur2"`
Expected: A count > 0, showing Biosaur2 parameters are registered.

Run: `OpenMS-build/bin/ProteomicsLFQ --helphelp 2>&1 | grep "Seeding:algorithm"`
Expected: Shows the new parameter.

- [ ] **Step 4: Commit**

```bash
git add src/topp/ProteomicsLFQ.cpp
git commit -m "feat(ProteomicsLFQ): register Seeding:algorithm and Biosaur2 param subsection"
```

---

### Task 3: Implement loadAndPreprocess_ Method

**Files:**
- Modify: `src/topp/ProteomicsLFQ.cpp:304-378` (add new method, refactor existing)

This task creates a new `loadAndPreprocess_()` method that branches on file type, and updates `quantifyFraction_()` to call it instead of `centroidAndCorrectPrecursors_()`.

- [ ] **Step 1: Add the loadAndPreprocess_ method**

Add the following method immediately before the existing `centroidAndCorrectPrecursors_()` method (before line 304):

```cpp
  ExitCodes loadAndPreprocess_(const String& mz_file, MSExperiment& ms_out, bool& is_im_peak_data)
  {
    const FileTypes::Type file_type = FileHandler::getType(mz_file);

    if (file_type == FileTypes::BRUKER_TDF)
    {
      // .d path: load with IM float arrays preserved
      FileHandler().loadExperiment(mz_file, ms_out, {FileTypes::BRUKER_TDF}, log_type_);
      ms_out.updateRanges();

      if (ms_out.empty())
      {
        OPENMS_LOG_WARN << "The given file does not contain any spectra.";
        return INCOMPATIBLE_INPUT_DATA;
      }

      // Remove MS2 peak data and sort (same housekeeping as mzML path)
      for (auto& spec : ms_out)
      {
        if (spec.getMSLevel() == 2)
        {
          spec.clear(false);
        }
        if (!spec.isSorted())
        {
          spec.sortByPosition();
        }
      }

      // No PeakPickerHiRes — Bruker TOF data is centroid-like; IM arrays must be preserved
      // No PrecursorCorrection — findHighestInWindow() is IM-unaware
      // No clearMetaDataArrays — IM per-peak float arrays must survive
      is_im_peak_data = true;
      OPENMS_LOG_INFO << "Loaded Bruker .d file with IM_PEAK data. Skipping centroiding and precursor correction.\n";
      return EXECUTION_OK;
    }
    else
    {
      // mzML path: existing centroid + precursor correction
      is_im_peak_data = false;
      return centroidAndCorrectPrecursors_(mz_file, ms_out);
    }
  }
```

- [ ] **Step 2: Update quantifyFraction_ to call loadAndPreprocess_**

In `quantifyFraction_()`, replace lines 886 and 897-900:

Replace:

```cpp
      MSExperiment ms_centroided;
```

With:

```cpp
      MSExperiment ms_centroided;
      bool is_im_peak_data = false;
```

Replace:

```cpp
      if (requires_ms_data)
      {
        ExitCodes e = centroidAndCorrectPrecursors_(mz_file, ms_centroided);
        if (e != EXECUTION_OK) { return e; }
```

With:

```cpp
      if (requires_ms_data)
      {
        ExitCodes e = loadAndPreprocess_(mz_file, ms_centroided, is_im_peak_data);
        if (e != EXECUTION_OK) { return e; }
```

- [ ] **Step 3: Build**

Run: `cmake --build OpenMS-build --target ProteomicsLFQ -j$(nproc)`
Expected: Compiles successfully.

- [ ] **Step 4: Run existing tests to verify mzML path is unchanged**

Run: `ctest --test-dir OpenMS-build -R "TOPP_ProteomicsLFQ_9" --output-on-failure`
Expected: Passes (single-file mzML test, simplest and fastest).

- [ ] **Step 5: Commit**

```bash
git add src/topp/ProteomicsLFQ.cpp
git commit -m "feat(ProteomicsLFQ): add loadAndPreprocess_ with .d/IM_PEAK branch"
```

---

### Task 4: Implement Biosaur2 Seed Generation Branch

**Files:**
- Modify: `src/topp/ProteomicsLFQ.cpp:897-930` (seed generation in quantifyFraction_)

- [ ] **Step 1: Disable mass_recalibration for .d input**

After the existing check at lines 889-893, add a check for IM_PEAK + mass_recalibration:

```cpp
      if (is_im_peak_data && mass_recalibration)
      {
        OPENMS_LOG_WARN << "Warning: mass_recalibration is not supported for .d input. Disabling.\n";
      }
```

Then wrap the existing recalibration block (lines 904-908) to also check `!is_im_peak_data`:

Replace:

```cpp
        if (mass_recalibration)
```

With:

```cpp
        if (mass_recalibration && !is_im_peak_data)
```

- [ ] **Step 2: Branch FWHM estimation and seed generation for .d**

Replace lines 910-930 (from `median_fwhm = DDAWorkflowCommons::...` through the end of the seed generation block) with:

```cpp
        if (!is_im_peak_data)
        {
          median_fwhm = DDAWorkflowCommons::estimateMedianChromatographicFWHM(ms_centroided);
          OPENMS_LOG_INFO << "Median chromatographic FWHM: " << median_fwhm << "\n";
        }
        // For .d/IM_PEAK: FWHM is estimated from Biosaur2 features below (after seed generation)
      }

      StringList id_msfile_ref;
      protein_ids[0].getPrimaryMSRunPath(id_msfile_ref);
      id_MS_run_ref.push_back(id_msfile_ref[0]);

      FeatureMap seeds;
      seeds.setPrimaryMSRunPath({mz_file});

      const bool targeted_only = getStringOption_("targeted_only") != "false";

      if (!targeted_only && in_feat_list.empty())
      {
        String seeding_algorithm = getStringOption_("Seeding:algorithm");

        // Force biosaur2 for IM_PEAK data (Multiplex cannot handle it)
        if (is_im_peak_data && seeding_algorithm != "biosaur2")
        {
          OPENMS_LOG_WARN << "Warning: IM_PEAK data detected. Forcing Seeding:algorithm to 'biosaur2' "
                          << "(FeatureFinderMultiplex does not support IM_PEAK format).\n";
          seeding_algorithm = "biosaur2";
        }

        if (seeding_algorithm == "biosaur2")
        {
          OPENMS_LOG_INFO << "Using Biosaur2Algorithm for seed detection.\n";
          Biosaur2Algorithm bio;
          Param bio_param = getParam_().copy("Seeding:Biosaur2:", true);
          bio.setParameters(bio_param);
          bio.setMSData(ms_centroided); // copy — run() destructively consumes internal data
          std::vector<Biosaur2Algorithm::Hill> hills;
          std::vector<Biosaur2Algorithm::PeptideFeature> peptide_features;
          bio.run(seeds, hills, peptide_features);
          OPENMS_LOG_INFO << "Biosaur2 produced " << seeds.size() << " seed features.\n";

          // For .d/IM_PEAK: estimate FWHM from Biosaur2 PeptideFeatures
          if (is_im_peak_data)
          {
            std::vector<double> fwhm_values;
            fwhm_values.reserve(peptide_features.size());
            for (const auto& pf : peptide_features)
            {
              double fwhm = pf.rt_end - pf.rt_start;
              if (fwhm > 0.0) { fwhm_values.push_back(fwhm); }
            }
            if (fwhm_values.size() >= 10)
            {
              median_fwhm = Math::median(fwhm_values.begin(), fwhm_values.end());
            }
            else
            {
              median_fwhm = 30.0; // fallback
              OPENMS_LOG_WARN << "Warning: Too few Biosaur2 features (" << fwhm_values.size()
                              << ") to estimate FWHM. Using default: " << median_fwhm << " seconds.\n";
            }
            OPENMS_LOG_INFO << "Median chromatographic FWHM (from Biosaur2): " << median_fwhm << "\n";
          }
        }
        else
        {
          DDAWorkflowCommons::calculateSeeds(ms_centroided, getDoubleOption_("Seeding:intThreshold"), seeds, median_fwhm, 2, 5);
        }

        if (debug_level_ > 666)
        {
          FileHandler().storeFeatures("debug_seeds_fraction_" + String(ms_files.first) + "_" + String(fraction_group) + ".featureXML", seeds, {FileTypes::FEATUREXML}, log_type_);
        }
      }
```

Note: This replaces lines 910 through 930 of the original file. The existing lines 914-921 (`StringList id_msfile_ref` through `const bool targeted_only`) are preserved unchanged within the replacement.

- [ ] **Step 3: Build**

Run: `cmake --build OpenMS-build --target ProteomicsLFQ -j$(nproc)`
Expected: Compiles successfully.

- [ ] **Step 4: Run existing tests to verify mzML path unchanged**

Run: `ctest --test-dir OpenMS-build -R "TOPP_ProteomicsLFQ_9" --output-on-failure`
Expected: Passes (default `Seeding:algorithm = "multiplex"` preserves existing behavior).

- [ ] **Step 5: Commit**

```bash
git add src/topp/ProteomicsLFQ.cpp
git commit -m "feat(ProteomicsLFQ): implement biosaur2 seed generation and .d FWHM estimation"
```

---

### Task 5: Preserve IM Meta Values on Features

**Files:**
- Modify: `src/topp/ProteomicsLFQ.cpp:1125` (feature meta value keep set)

- [ ] **Step 1: Extend the keep_meta set**

Replace line 1125:

```cpp
      unordered_set<String> keep_meta = {"OffsetPeptide"};
```

With:

```cpp
      unordered_set<String> keep_meta = {"OffsetPeptide", "IM_median", "IM_min", "IM_max"};
```

- [ ] **Step 2: Build**

Run: `cmake --build OpenMS-build --target ProteomicsLFQ -j$(nproc)`
Expected: Compiles successfully.

- [ ] **Step 3: Run existing tests**

Run: `ctest --test-dir OpenMS-build -R "TOPP_ProteomicsLFQ_9" --output-on-failure`
Expected: Passes. For mzML data without IM, these keys don't exist on features so nothing changes.

- [ ] **Step 4: Commit**

```bash
git add src/topp/ProteomicsLFQ.cpp
git commit -m "feat(ProteomicsLFQ): preserve IM meta values on features"
```

---

### Task 6: Run Full Existing Test Suite

**Files:**
- No changes — verification only

- [ ] **Step 1: Build all dependencies**

Run: `cmake --build OpenMS-build -j$(nproc)`
Expected: Full build succeeds.

- [ ] **Step 2: Run all existing ProteomicsLFQ tests**

Run: `ctest --test-dir OpenMS-build -R "TOPP_ProteomicsLFQ" --output-on-failure`
Expected: All existing tests pass unchanged. The default `Seeding:algorithm = "multiplex"` and mzML input means no new code path is exercised.

- [ ] **Step 3: If any tests fail, investigate**

Do NOT update reference files. These tests use mzML input, so only the unchanged mzML code path should run. Any failure indicates a regression in the existing path.

---

### Task 7: Add Biosaur2 Seeding Smoke Test (mzML)

**Files:**
- Modify: `src/tests/topp/CMakeLists.txt:3419` (after TOPP_ProteomicsLFQ_9 test block)

This test verifies the biosaur2 seeding path works on existing mzML test data (no .d test data needed).

- [ ] **Step 1: Add the smoke test**

After line 3419 (`set_tests_properties("TOPP_ProteomicsLFQ_9_out_4" ...)`), add:

```cmake
# biosaur2 seeding algorithm smoke test (mzML input, verifies pipeline completes)
add_test("TOPP_ProteomicsLFQ_biosaur2_seeds" ${TOPP_BIN_PATH}/ProteomicsLFQ
         -in
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F2.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F2.mzML
         -ids
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F2.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F2.idXML
         -Alignment:align_algorithm:max_rt_shift 0
         -Linking:min_nr_diffs_per_bin 10
         -fasta
         ${DATA_DIR_SHARE}/examples/TOPPAS/data/BSA_Identification/18Protein_SoCe_Tr_detergents_trace_target_decoy.fasta
         -targeted_only false
         -mass_recalibration false
         -Seeding:algorithm biosaur2
         -out_cxml ${TESTS_TEMP_DIR}/BSA_biosaur2_seeds.tmp.consensusXML
         -out ${TESTS_TEMP_DIR}/BSA_biosaur2_seeds.tmp.mzTab
         -threads 1
         -proteinFDR 0.8
         -test
         )
```

This is a smoke test — it only verifies the pipeline completes without error when using biosaur2 seeds. No reference file comparison since biosaur2 produces different seeds than multiplex.

- [ ] **Step 2: Build and run**

Run: `cmake --build OpenMS-build -j$(nproc)`
Then: `ctest --test-dir OpenMS-build -R "TOPP_ProteomicsLFQ_biosaur2_seeds" --output-on-failure`
Expected: Test passes (exit code 0).

- [ ] **Step 3: Commit**

```bash
git add src/tests/topp/CMakeLists.txt
git commit -m "test(ProteomicsLFQ): add biosaur2 seeding smoke test"
```

---

### Task 8: Add .d Integration Test (Optional, requires WITH_OPENTIMS)

**Files:**
- Modify: `src/tests/topp/CMakeLists.txt:2992` (inside the `if (WITH_OPENTIMS AND OPENTIMS_DDA_TEST_DATA)` block, after SageAdapter FDR test)

This test exercises the full .d code path: Bruker `.d` input with Sage-generated idXML containing IM annotations. It is guarded by `WITH_OPENTIMS` and `SAGE_BINARY` since it depends on the SageAdapter FDR-filtered output from the existing test chain.

- [ ] **Step 1: Add the .d ProteomicsLFQ integration test**

After line 2993 (`endif()` closing the Sage block, but still inside the outer `if (WITH_OPENTIMS AND OPENTIMS_DDA_TEST_DATA)` block — i.e. before the line 2994 `endif()`), add:

```cmake
    # ProteomicsLFQ with Bruker .d input (DDA-PASEF, smoke test — verifies .d path completes)
    add_test("TOPP_ProteomicsLFQ_DDA_PASEF" ${TOPP_BIN_PATH}/ProteomicsLFQ -test
      -in ${OPENTIMS_DDA_SYMLINK}
      -ids ${TESTS_TEMP_DIR}/SageAdapter_DDA_fdr.tmp.idXML
      -fasta ${TESTS_TEMP_DIR}/DecoyDatabase_DDA_out.tmp.fasta
      -targeted_only true
      -mass_recalibration false
      -out ${TESTS_TEMP_DIR}/ProteomicsLFQ_DDA_PASEF.tmp.mzTab
      -out_cxml ${TESTS_TEMP_DIR}/ProteomicsLFQ_DDA_PASEF.tmp.consensusXML
      -proteinFDR 1.0
      -threads 1)
    set_tests_properties("TOPP_ProteomicsLFQ_DDA_PASEF" PROPERTIES DEPENDS "TOPP_FalseDiscoveryRate_SageDDA")
```

This is a smoke test — it verifies the `.d` loading path, IM_PEAK detection, and FFId with IM annotations all complete without error. Uses `targeted_only true` to skip seed generation (simpler first test). Uses `proteinFDR 1.0` to avoid filtering out all proteins from a small test dataset.

- [ ] **Step 2: Build and run (only possible if WITH_OPENTIMS is enabled and Sage is available)**

Run: `cmake --build OpenMS-build -j$(nproc)`
Then: `ctest --test-dir OpenMS-build -R "TOPP_ProteomicsLFQ_DDA_PASEF" --output-on-failure`
Expected: Test passes (exit code 0). If `WITH_OPENTIMS` is not enabled, the test won't exist — that's fine.

- [ ] **Step 3: Commit**

```bash
git add src/tests/topp/CMakeLists.txt
git commit -m "test(ProteomicsLFQ): add optional .d DDA-PASEF integration test"
```

---

### Task 9: Update Documentation

**Files:**
- Modify: `src/topp/ProteomicsLFQ.cpp:70-123` (Doxygen documentation block)

- [ ] **Step 1: Update the tool documentation**

Replace line 75 (`  - Spectra in mzML format`) with:

```cpp
  - Spectra in mzML format or Bruker .d directories (TimsTOF PASEF)
```

After line 101 (`No special preparation of the input mzML file is required.`), add:

```cpp

@b Bruker .d (TimsTOF PASEF): @n
Bruker .d directories containing DDA-PASEF data are supported directly.
When .d input is detected, the tool automatically:
  - Skips centroiding (PeakPickerHiRes) to preserve per-peak ion mobility data
  - Skips precursor mass correction (not IM-aware)
  - Forces Biosaur2Algorithm for seed generation (FeatureFinderMultiplex does not support IM_PEAK)
  - Estimates chromatographic FWHM from Biosaur2 feature extents
Identification files should be generated with SageAdapter, which annotates
ion mobility values in the idXML output. FeatureFinderIdentificationAlgorithm
uses these IM annotations for targeted 2D chromatogram extraction (m/z + IM windowing).
```

- [ ] **Step 2: Build to verify no doc errors**

Run: `cmake --build OpenMS-build --target ProteomicsLFQ -j$(nproc)`
Expected: Compiles successfully.

- [ ] **Step 3: Commit**

```bash
git add src/topp/ProteomicsLFQ.cpp
git commit -m "docs(ProteomicsLFQ): document Bruker .d and Biosaur2 seeding support"
```
