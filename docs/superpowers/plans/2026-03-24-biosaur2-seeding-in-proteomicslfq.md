# Biosaur2 Seeding in ProteomicsLFQ Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a user-selectable `Seeding:algorithm` parameter to ProteomicsLFQ so users can choose between the existing multiplex-based seeding and Biosaur2Algorithm for untargeted feature seed generation.

**Architecture:** Single-file change in `src/topp/ProteomicsLFQ.cpp`. Add a `Seeding:algorithm` string parameter with `"multiplex"` (default) and `"biosaur2"` choices. Branch in `quantifyFraction_()` to call the chosen algorithm. Biosaur2 parameters go under `Seeding:Biosaur2:` subsection. Downstream pipeline is unchanged.

**Tech Stack:** C++, OpenMS framework (TOPPBase parameter system, DefaultParamHandler, Biosaur2Algorithm, DDAWorkflowCommons)

**Spec:** `docs/superpowers/specs/2026-03-24-biosaur2-seeding-in-proteomicslfq-design.md`

---

### Task 1: Register Seeding:algorithm Parameter

**Files:**
- Modify: `src/topp/ProteomicsLFQ.cpp:240-243` (after existing Seeding subsection registration)

- [ ] **Step 1: Add the algorithm parameter after existing Seeding params**

After line 243 (`registerDoubleOption_("Seeding:traceRTTolerance", ...)`), add:

```cpp
    registerStringOption_("Seeding:algorithm", "<choice>", "multiplex",
      "Algorithm for untargeted seed feature detection.\n"
      "multiplex: FeatureFinderMultiplexAlgorithm (default, current behavior).\n"
      "biosaur2: Biosaur2Algorithm (C++ reimplementation of biosaur2 feature finder).",
      false, false);
    setValidStrings_("Seeding:algorithm", {"multiplex", "biosaur2"});
```

- [ ] **Step 2: Build and verify parameter registration**

Run: `cmake --build OpenMS-build --target ProteomicsLFQ -j$(nproc)`
Expected: Compiles successfully.

Then verify the parameter shows up:
Run: `OpenMS-build/bin/ProteomicsLFQ --helphelp 2>&1 | grep -A2 "Seeding:algorithm"`
Expected: Shows the new parameter with valid strings.

- [ ] **Step 3: Commit**

```bash
git add src/topp/ProteomicsLFQ.cpp
git commit -m "feat(ProteomicsLFQ): add Seeding:algorithm parameter for seed algo selection"
```

---

### Task 2: Register Seeding:Biosaur2 Parameter Subsection

**Files:**
- Modify: `src/topp/ProteomicsLFQ.cpp:1-50` (add include)
- Modify: `src/topp/ProteomicsLFQ.cpp:297-305` (add to combined Param)

- [ ] **Step 1: Add the Biosaur2Algorithm include**

At the top of the file, add among the other FEATUREFINDER includes (after line 31):

```cpp
#include <OpenMS/FEATUREFINDER/Biosaur2Algorithm.h>
```

- [ ] **Step 2: Add Biosaur2 defaults to the combined Param**

Before line 297 (`// combine parameters of the individual algorithms`), add:

```cpp
    Param bio_defaults = Biosaur2Algorithm().getDefaults();
    for (auto it = bio_defaults.begin(); it != bio_defaults.end(); ++it)
    {
      bio_defaults.addTag(it.getName(), "advanced");
    }
```

Then at line 303 (in the combined Param insertions, before `registerFullParam_(combined)`), add:

```cpp
    combined.insert("Seeding:Biosaur2:", bio_defaults);
```

- [ ] **Step 3: Build and verify biosaur2 parameters are visible**

Run: `cmake --build OpenMS-build --target ProteomicsLFQ -j$(nproc)`
Expected: Compiles successfully.

Run: `OpenMS-build/bin/ProteomicsLFQ --helphelp 2>&1 | grep "Seeding:Biosaur2"`
Expected: Shows Biosaur2-prefixed parameters (e.g., `Seeding:Biosaur2:htol`, `Seeding:Biosaur2:mini`, etc.).

- [ ] **Step 4: Commit**

```bash
git add src/topp/ProteomicsLFQ.cpp
git commit -m "feat(ProteomicsLFQ): register Biosaur2Algorithm params under Seeding:Biosaur2"
```

---

### Task 3: Implement Seed Generation Branch

**Files:**
- Modify: `src/topp/ProteomicsLFQ.cpp:927-934` (seed generation in quantifyFraction_)

- [ ] **Step 1: Replace the unconditional calculateSeeds call with a branch**

Replace lines 927-934:

```cpp
      if (!targeted_only && in_feat_list.empty())
      {
        DDAWorkflowCommons::calculateSeeds(ms_centroided, getDoubleOption_("Seeding:intThreshold"), seeds, median_fwhm, 2, 5);
        if (debug_level_ > 666)
        {
          FileHandler().storeFeatures("debug_seeds_fraction_" + String(ms_files.first) + "_" + String(fraction_group) + ".featureXML", seeds, {FileTypes::FEATUREXML}, log_type_);
        }
      }
```

With:

```cpp
      if (!targeted_only && in_feat_list.empty())
      {
        const String seeding_algorithm = getStringOption_("Seeding:algorithm");
        if (seeding_algorithm == "biosaur2")
        {
          OPENMS_LOG_INFO << "Using Biosaur2Algorithm for seed detection.\n";
          Biosaur2Algorithm bio;
          Param bio_param = getParam_().copy("Seeding:Biosaur2:", true);
          bio.setParameters(bio_param);
          bio.setMSData(ms_centroided); // copy — run() destructively consumes its internal MSExperiment
          bio.run(seeds);
          OPENMS_LOG_INFO << "Using " << seeds.size() << " seeds from Biosaur2 feature detection.\n";
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

**Important notes for implementer:**
- `bio.setMSData(ms_centroided)` uses the **copy** overload (not move). This is intentional because `Biosaur2Algorithm::run()` erases non-MS1 spectra and moves internal data into FAIMS groups via `std::move(ms_data_)`. After `run()`, the internal MSExperiment is empty. `ms_centroided` is still needed downstream by `FeatureFinderIdentificationAlgorithm`.
- `getParam_().copy("Seeding:Biosaur2:", true)` strips the `Seeding:Biosaur2:` prefix, yielding parameter names that match `Biosaur2Algorithm`'s `DefaultParamHandler` keys.
- Biosaur2's `profile_mode` already defaults to `false`, so no explicit override needed (data is already centroided by ProteomicsLFQ).

- [ ] **Step 2: Build**

Run: `cmake --build OpenMS-build --target ProteomicsLFQ -j$(nproc)`
Expected: Compiles successfully.

- [ ] **Step 3: Commit**

```bash
git add src/topp/ProteomicsLFQ.cpp
git commit -m "feat(ProteomicsLFQ): implement biosaur2 seed generation branch"
```

---

### Task 4: Verify Existing Tests Still Pass (Default Path)

**Files:**
- No changes — verification only

- [ ] **Step 1: Build all dependencies**

Run: `cmake --build OpenMS-build -j$(nproc)`
Expected: Full build succeeds.

- [ ] **Step 2: Run existing ProteomicsLFQ tests**

Run: `ctest --test-dir OpenMS-build -R "TOPP_ProteomicsLFQ" --output-on-failure`
Expected: All existing tests pass. The default `Seeding:algorithm = "multiplex"` preserves the current behavior exactly.

- [ ] **Step 3: Note results**

If any tests fail, investigate — these use `targeted_only true` or the default multiplex path, neither of which should be affected. Do not update reference files without root cause analysis.

---

### Task 5: Add Biosaur2 Seeding Test

**Files:**
- Modify: `src/tests/topp/CMakeLists.txt:~3157` (after TOPP_ProteomicsLFQ_3 test block)

- [ ] **Step 1: Add a new test that exercises the biosaur2 seeding path**

After the `TOPP_ProteomicsLFQ_3` test block (after line 3158), add a new test. This reuses the same input data as test 3 (which uses `targeted_only false`, meaning seeds are generated) but switches the seeding algorithm:

```cmake
# biosaur2 seeding algorithm (smoke test — verifies pipeline completes with biosaur2 seeds)
add_test("TOPP_ProteomicsLFQ_biosaur2_seeds" ${TOPP_BIN_PATH}/ProteomicsLFQ
         -in
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F1.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F2.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F1.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F2.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA3_F1.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA3_F2.mzML
         -ids
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F1.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F2.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F1.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F2.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA3_F1.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA3_F2.idXML
         -design
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA_design.tsv
         -Seeding:algorithm biosaur2
         -Alignment:align_algorithm:max_rt_shift 0
         -Linking:min_nr_diffs_per_bin 10
         -fasta
         ${DATA_DIR_SHARE}/examples/TOPPAS/data/BSA_Identification/18Protein_SoCe_Tr_detergents_trace_target_decoy.fasta
         -targeted_only false
         -mass_recalibration false
         -out_cxml ${TESTS_TEMP_DIR}/BSA_biosaur2_seeds.tmp.consensusXML
         -out ${TESTS_TEMP_DIR}/BSA_biosaur2_seeds.tmp.mzTab
         -threads 1
         -proteinFDR 0.3
         -test
         )
```

This is a smoke test — it verifies the pipeline completes without error when using biosaur2 seeds. No reference file comparison since biosaur2 produces different seeds than multiplex (different feature counts expected).

- [ ] **Step 2: Build and run the new test**

Run: `cmake --build OpenMS-build -j$(nproc)`
Then: `ctest --test-dir OpenMS-build -R "TOPP_ProteomicsLFQ_biosaur2" --output-on-failure`
Expected: Test passes (exit code 0).

- [ ] **Step 3: Commit**

```bash
git add src/tests/topp/CMakeLists.txt
git commit -m "test(ProteomicsLFQ): add smoke test for biosaur2 seeding algorithm"
```
