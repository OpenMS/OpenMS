# Multimer Detection Test & Functional Gaps Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Fix the same-adduct multimer functional gap and add comprehensive test coverage for all untested multimer detection scenarios.

**Architecture:** Add conditional identity compomers to `MassExplainer::compute()` (gated by `include_identity` parameter, only active when `max_multimer > 1`). Then add 5 targeted test sections to the existing test files covering same-adduct multimers, annotation strings, trimers, penalty effectiveness, and negative mode.

**Tech Stack:** C++, OpenMS test framework (START_SECTION/END_SECTION, TEST_EQUAL, TEST_REAL_SIMILAR)

---

### Task 1: Add `include_identity` parameter to `MassExplainer::compute()`

**Files:**
- Modify: `src/openms/include/OpenMS/DATASTRUCTURES/MassExplainer.h:124`
- Modify: `src/openms/source/DATASTRUCTURES/MassExplainer.cpp:137,220`
- Test: `src/tests/class_tests/openms/source/MassExplainer_test.cpp`

**Step 1: Write failing test for identity compomers**

Add before `END_TEST` in `MassExplainer_test.cpp` (line 174):

```cpp
START_SECTION((void compute(bool include_identity)))
{
  // Without identity compomers (default), no mass=0 net_charge=0 compomers exist
  MassExplainer me_no_id;
  me_no_id.compute();  // include_identity defaults to false

  MassExplainer::CompomerIterator s, e;
  // mass=0, net_charge=0 should have no hits without identity
  SignedSize hits_no_id = me_no_id.query(0, 0.0, 0.5, -100000, s, e);
  TEST_EQUAL(hits_no_id, 0);

  // With identity compomers, should find compomers at mass=0, net_charge=0
  MassExplainer me_id;
  me_id.compute(true);  // include_identity = true

  SignedSize hits_id = me_id.query(0, 0.0, 0.5, -100000, s, e);
  TEST_EQUAL(hits_id > 0, true);

  // Each identity compomer should have getSideMass(LEFT) == getSideMass(RIGHT)
  for (auto it = s; it != e; ++it)
  {
    TEST_REAL_SIMILAR(it->getSideMass(Compomer::LEFT), it->getSideMass(Compomer::RIGHT));
    TEST_EQUAL(it->getNetCharge(), 0);
  }

  // Same-adduct queryMultimer should now work
  // H+ identity: expected = 2*H_mass - 1*H_mass = H_mass
  // Use default adducts (H+, Na+, NH4+, K+)
  std::vector<MassExplainer::CompomerIterator> multimer_hits;
  // Query with approximate H+ mass (proton mass ~ 1.007)
  double H_mass = 1.00728;
  double m1 = 100.0 + H_mass;       // [M+H]+ with M=100
  double m2 = 200.0 + H_mass;       // [2M+H]+ with M=100
  SignedSize n_hits = me_id.queryMultimer(0, m1, m2, 1, 2, 0.1, -100000, multimer_hits);
  TEST_EQUAL(n_hits > 0, true);

  // Without identity, same query should fail
  multimer_hits.clear();
  n_hits = me_no_id.queryMultimer(0, m1, m2, 1, 2, 0.1, -100000, multimer_hits);
  TEST_EQUAL(n_hits, 0);
}
END_SECTION
```

**Step 2: Run test to verify it fails**

Run: `/snap/cmake/1515/bin/cmake --build OpenMS-build --target MassExplainer_test -j$(nproc) && OpenMS-build/src/tests/class_tests/bin/MassExplainer_test`
Expected: Build FAIL — `compute()` doesn't accept `bool` parameter yet

**Step 3: Implement identity compomers**

In `MassExplainer.h` line 124, change the declaration:
```cpp
void compute(bool include_identity = false);
```

In `MassExplainer.cpp` line 137, change the signature:
```cpp
void MassExplainer::compute(bool include_identity)
```

In `MassExplainer.cpp` after line 220 (`explanations_.swap(valids_only);`), add:
```cpp
    // Add identity compomers (same adduct on both LEFT and RIGHT sides).
    // These have net_charge=0, mass=0, and are needed for same-adduct multimer
    // detection (e.g., [M+H]+ paired with [2M+H]+). Gated by include_identity
    // to avoid changing behavior for non-multimer workflows.
    if (include_identity)
    {
      for (const auto& adduct : adduct_charged)
      {
        Compomer cmpi;
        cmpi.add(adduct, Compomer::LEFT);
        cmpi.add(adduct, Compomer::RIGHT);
        if (compomerValid_(cmpi))
        {
          explanations_.push_back(cmpi);
        }
      }
    }
```

**Step 4: Run test to verify it passes**

Run: `/snap/cmake/1515/bin/cmake --build OpenMS-build --target MassExplainer_test -j$(nproc) && OpenMS-build/src/tests/class_tests/bin/MassExplainer_test`
Expected: PASS

**Step 5: Commit**

```bash
git add src/openms/include/OpenMS/DATASTRUCTURES/MassExplainer.h \
        src/openms/source/DATASTRUCTURES/MassExplainer.cpp \
        src/tests/class_tests/openms/source/MassExplainer_test.cpp
git commit -m "feat: add conditional identity compomers to MassExplainer for same-adduct multimers"
```

---

### Task 2: Wire `include_identity` into `candidateEdges_()`

**Files:**
- Modify: `src/openms/source/ANALYSIS/DECHARGING/MetaboliteFeatureDeconvolution.cpp:417`

**Step 1: Change `me.compute()` call**

At line 417, change:
```cpp
    me.compute();
```
to:
```cpp
    me.compute(/*include_identity=*/ max_multimer > 1);
```

Note: `max_multimer` is already declared at line 370 (`Int max_multimer = param_.getValue("max_multimer");`).

**Step 2: Build and run all related tests**

Run: `/snap/cmake/1515/bin/cmake --build OpenMS-build --target MetaboliteFeatureDeconvolution_test -j$(nproc) && ctest --test-dir OpenMS-build -R "DECHARGING|Adduct|Compomer|ChargePair|MassExplainer" --output-on-failure`
Expected: All 12 tests PASS (no behavioral change for existing tests since max_multimer defaults to 1)

**Step 3: Commit**

```bash
git add src/openms/source/ANALYSIS/DECHARGING/MetaboliteFeatureDeconvolution.cpp
git commit -m "feat: enable identity compomers when max_multimer > 1"
```

---

### Task 3: Test same-adduct multimer detection `[M+H]+` / `[2M+H]+`

This is the most important test — it verifies the functional gap is fixed.

**Files:**
- Modify: `src/tests/class_tests/openms/source/MetaboliteFeatureDeconvolution_test.cpp`

**Step 1: Add same-adduct multimer test section**

Add before the final `END_TEST` (after the existing `[EXTRA] Multimer detection` section):

```cpp
START_SECTION(([EXTRA] Same-adduct multimer detection))
{
  // Same-adduct multimer: [M+H]+ (monomer) and [2M+H]+ (dimer)
  // This requires identity compomers (H+LEFT / H+RIGHT) in MassExplainer.
  // With identity, the equation gives: 2*H - H = H, which matches.
  double H_adduct = 1.00728;   // approximate proton mass
  double M = 500.0;

  FeatureMap fm;
  Feature f1;
  f1.setMZ(M + H_adduct);           // monomer [M+H]+
  f1.setRT(100.0);
  f1.setIntensity(1000.0f);
  f1.setCharge(1);
  fm.push_back(f1);

  Feature f2;
  f2.setMZ(2.0 * M + H_adduct);     // dimer [2M+H]+
  f2.setRT(100.0);
  f2.setIntensity(500.0f);
  f2.setCharge(1);
  fm.push_back(f2);

  MetaboliteFeatureDeconvolution mfd;
  Param p = mfd.getDefaults();
  p.setValue("charge_min", 1);
  p.setValue("charge_max", 1);
  p.setValue("charge_span_max", 1);
  p.setValue("max_multimer", 2);
  p.setValue("potential_adducts", std::vector<std::string>{"H:+:1.0"});
  p.setValue("mass_max_diff", 0.5);
  p.setValue("retention_max_diff", 2.0);
  p.setValue("retention_max_diff_local", 2.0);
  mfd.setParameters(p);

  FeatureMap fm_out;
  ConsensusMap cm_out, cm_out2;
  mfd.compute(fm, fm_out, cm_out, cm_out2);

  // The two features should be grouped
  bool found_group = false;
  for (Size i = 0; i < cm_out.size(); ++i)
  {
    if (cm_out[i].size() >= 2) found_group = true;
  }
  TEST_EQUAL(found_group, true);

  // The dimer feature should have mol_multiplier=2
  bool found_dimer = false;
  for (Size i = 0; i < fm_out.size(); ++i)
  {
    if (fm_out[i].metaValueExists("mol_multiplier"))
    {
      TEST_EQUAL((Int)fm_out[i].getMetaValue("mol_multiplier"), 2);
      found_dimer = true;
    }
  }
  TEST_EQUAL(found_dimer, true);
}
END_SECTION
```

**Step 2: Build and run**

Run: `/snap/cmake/1515/bin/cmake --build OpenMS-build --target MetaboliteFeatureDeconvolution_test -j$(nproc) && OpenMS-build/src/tests/class_tests/bin/MetaboliteFeatureDeconvolution_test`
Expected: PASS — identity compomers enable same-adduct detection with single H+ adduct and charge_max=1

**Step 3: Commit**

```bash
git add src/tests/class_tests/openms/source/MetaboliteFeatureDeconvolution_test.cpp
git commit -m "test: verify same-adduct multimer detection [M+H]+ / [2M+H]+"
```

---

### Task 4: Test annotation string content and trimer detection

**Files:**
- Modify: `src/tests/class_tests/openms/source/MetaboliteFeatureDeconvolution_test.cpp`

**Step 1: Add annotation and trimer test section**

```cpp
START_SECTION(([EXTRA] Multimer annotation strings and trimer detection))
{
  // Test that annotation strings are correct and trimers (n=3) work.
  // Three features: monomer, dimer, trimer — all with H+ adduct.
  double H_adduct = 1.00728;
  double M = 300.0;

  FeatureMap fm;
  Feature f1;
  f1.setMZ(M + H_adduct);
  f1.setRT(100.0);
  f1.setIntensity(1000.0f);
  f1.setCharge(1);
  fm.push_back(f1);

  Feature f2;
  f2.setMZ(2.0 * M + H_adduct);
  f2.setRT(100.0);
  f2.setIntensity(500.0f);
  f2.setCharge(1);
  fm.push_back(f2);

  Feature f3;
  f3.setMZ(3.0 * M + H_adduct);
  f3.setRT(100.0);
  f3.setIntensity(300.0f);
  f3.setCharge(1);
  fm.push_back(f3);

  MetaboliteFeatureDeconvolution mfd;
  Param p = mfd.getDefaults();
  p.setValue("charge_min", 1);
  p.setValue("charge_max", 1);
  p.setValue("charge_span_max", 1);
  p.setValue("max_multimer", 3);
  p.setValue("potential_adducts", std::vector<std::string>{"H:+:1.0"});
  p.setValue("mass_max_diff", 0.5);
  p.setValue("retention_max_diff", 2.0);
  p.setValue("retention_max_diff_local", 2.0);
  mfd.setParameters(p);

  FeatureMap fm_out;
  ConsensusMap cm_out, cm_out2;
  mfd.compute(fm, fm_out, cm_out, cm_out2);

  // Check annotation content
  bool found_dimer_annotation = false;
  bool found_trimer_annotation = false;
  for (Size i = 0; i < fm_out.size(); ++i)
  {
    if (fm_out[i].metaValueExists("mol_multiplier"))
    {
      Int mult = (Int)fm_out[i].getMetaValue("mol_multiplier");
      if (mult == 2)
      {
        found_dimer_annotation = true;
        // Verify the adducts string contains "2M"
        if (fm_out[i].metaValueExists("adducts"))
        {
          StringList adducts = fm_out[i].getMetaValue("adducts");
          TEST_EQUAL(adducts[0].hasSubstring("2M"), true);
        }
      }
      if (mult == 3)
      {
        found_trimer_annotation = true;
        // Verify the adducts string contains "3M"
        if (fm_out[i].metaValueExists("adducts"))
        {
          StringList adducts = fm_out[i].getMetaValue("adducts");
          TEST_EQUAL(adducts[0].hasSubstring("3M"), true);
        }
      }
    }
  }
  TEST_EQUAL(found_dimer_annotation, true);
  TEST_EQUAL(found_trimer_annotation, true);
}
END_SECTION
```

**Step 2: Build and run**

Run: `/snap/cmake/1515/bin/cmake --build OpenMS-build --target MetaboliteFeatureDeconvolution_test -j$(nproc) && OpenMS-build/src/tests/class_tests/bin/MetaboliteFeatureDeconvolution_test`
Expected: PASS

**Step 3: Commit**

```bash
git add src/tests/class_tests/openms/source/MetaboliteFeatureDeconvolution_test.cpp
git commit -m "test: verify annotation strings and trimer detection (n=3)"
```

---

### Task 5: Test penalty effectiveness

**Files:**
- Modify: `src/tests/class_tests/openms/source/MetaboliteFeatureDeconvolution_test.cpp`

**Step 1: Add penalty test section**

The test creates a scenario where both a monomer explanation and a dimer explanation are possible for the same feature pair, and verifies that a very strong penalty suppresses the dimer interpretation.

```cpp
START_SECTION(([EXTRA] Multimer penalty suppresses dimer when penalty is extreme))
{
  // Create features where multimer detection is possible.
  // With extreme penalty (e.g., -100), the ILP should not select multimer edges.
  double H_adduct = 1.00728;
  double M = 500.0;

  FeatureMap fm;
  Feature f1;
  f1.setMZ(M + H_adduct);
  f1.setRT(100.0);
  f1.setIntensity(1000.0f);
  f1.setCharge(1);
  fm.push_back(f1);

  Feature f2;
  f2.setMZ(2.0 * M + H_adduct);
  f2.setRT(100.0);
  f2.setIntensity(500.0f);
  f2.setCharge(1);
  fm.push_back(f2);

  // Extreme penalty: exp(-100) ~ 3.7e-44, effectively zero edge score
  MetaboliteFeatureDeconvolution mfd;
  Param p = mfd.getDefaults();
  p.setValue("charge_min", 1);
  p.setValue("charge_max", 1);
  p.setValue("charge_span_max", 1);
  p.setValue("max_multimer", 2);
  p.setValue("multimer_log_penalty", -100.0);
  p.setValue("potential_adducts", std::vector<std::string>{"H:+:1.0"});
  p.setValue("mass_max_diff", 0.5);
  p.setValue("retention_max_diff", 2.0);
  p.setValue("retention_max_diff_local", 2.0);
  mfd.setParameters(p);

  FeatureMap fm_out;
  ConsensusMap cm_out, cm_out2;
  mfd.compute(fm, fm_out, cm_out, cm_out2);

  // With extreme penalty, no grouping should occur
  bool found_group = false;
  for (Size i = 0; i < cm_out.size(); ++i)
  {
    if (cm_out[i].size() >= 2) found_group = true;
  }
  TEST_EQUAL(found_group, false);

  // No mol_multiplier annotation should exist
  bool found_multiplier = false;
  for (Size i = 0; i < fm_out.size(); ++i)
  {
    if (fm_out[i].metaValueExists("mol_multiplier")) found_multiplier = true;
  }
  TEST_EQUAL(found_multiplier, false);
}
END_SECTION
```

**Step 2: Build and run**

Run: `/snap/cmake/1515/bin/cmake --build OpenMS-build --target MetaboliteFeatureDeconvolution_test -j$(nproc) && OpenMS-build/src/tests/class_tests/bin/MetaboliteFeatureDeconvolution_test`
Expected: PASS

**Step 3: Commit**

```bash
git add src/tests/class_tests/openms/source/MetaboliteFeatureDeconvolution_test.cpp
git commit -m "test: verify extreme penalty suppresses multimer edge selection"
```

---

### Task 6: Test negative mode multimer detection

**Files:**
- Modify: `src/tests/class_tests/openms/source/MetaboliteFeatureDeconvolution_test.cpp`

**Step 1: Add negative mode test section**

```cpp
START_SECTION(([EXTRA] Negative mode multimer detection))
{
  // Negative mode: [M-H]- (monomer) and [2M-H]- (dimer)
  // H-1 adduct has charge -1 and mass = -PROTON_MASS_U
  double H_adduct = 1.00728;  // proton mass (will be subtracted in negative mode)
  double M = 500.0;

  FeatureMap fm;
  Feature f1;
  f1.setMZ(M - H_adduct);           // monomer [M-H]- (deprotonated)
  f1.setRT(100.0);
  f1.setIntensity(1000.0f);
  f1.setCharge(1);                   // FFM reports absolute charge
  fm.push_back(f1);

  Feature f2;
  f2.setMZ(2.0 * M - H_adduct);     // dimer [2M-H]-
  f2.setRT(100.0);
  f2.setIntensity(500.0f);
  f2.setCharge(1);
  fm.push_back(f2);

  MetaboliteFeatureDeconvolution mfd;
  Param p = mfd.getDefaults();
  p.setValue("charge_min", 1);
  p.setValue("charge_max", 1);
  p.setValue("charge_span_max", 1);
  p.setValue("max_multimer", 2);
  p.setValue("negative_mode", "true");
  p.setValue("potential_adducts", std::vector<std::string>{"H-1:-:1.0"});
  p.setValue("mass_max_diff", 0.5);
  p.setValue("retention_max_diff", 2.0);
  p.setValue("retention_max_diff_local", 2.0);
  mfd.setParameters(p);

  FeatureMap fm_out;
  ConsensusMap cm_out, cm_out2;
  mfd.compute(fm, fm_out, cm_out, cm_out2);

  // The two features should be grouped as monomer/dimer
  bool found_group = false;
  for (Size i = 0; i < cm_out.size(); ++i)
  {
    if (cm_out[i].size() >= 2) found_group = true;
  }
  TEST_EQUAL(found_group, true);

  // Dimer annotation should exist
  bool found_dimer = false;
  for (Size i = 0; i < fm_out.size(); ++i)
  {
    if (fm_out[i].metaValueExists("mol_multiplier"))
    {
      TEST_EQUAL((Int)fm_out[i].getMetaValue("mol_multiplier"), 2);
      found_dimer = true;
    }
  }
  TEST_EQUAL(found_dimer, true);
}
END_SECTION
```

**Step 2: Build and run**

Run: `/snap/cmake/1515/bin/cmake --build OpenMS-build --target MetaboliteFeatureDeconvolution_test -j$(nproc) && OpenMS-build/src/tests/class_tests/bin/MetaboliteFeatureDeconvolution_test`
Expected: PASS (or FAIL if negative mode identity compomers need special handling — debug if so)

**Step 3: Commit**

```bash
git add src/tests/class_tests/openms/source/MetaboliteFeatureDeconvolution_test.cpp
git commit -m "test: verify negative mode multimer detection [M-H]- / [2M-H]-"
```

---

### Task 7: Test mixed-charge multimer detection

**Files:**
- Modify: `src/tests/class_tests/openms/source/MetaboliteFeatureDeconvolution_test.cpp`

**Step 1: Add mixed-charge test section**

```cpp
START_SECTION(([EXTRA] Mixed-charge multimer detection))
{
  // Mixed charge: [M+H]+ (charge 1) and [2M+2H]2+ (charge 2)
  // Both use H+ adduct but at different charges.
  // mz1 = M + H = 501.007, q1 = 1, m1 = 501.007
  // mz2 = (2M + 2H) / 2 = M + H = 501.007, q2 = 2, m2 = 1002.014
  // For n1=1, n2=2: observed = 2*501.007 - 1*1002.014 = 0.0
  // Identity compomer H+L/H+R: expected = 2*H - 1*H = H ~ 1.007
  // These don't match (0 != 1.007). This is correct: the processHit
  // fills remaining charge with default adduct, creating H1L/H2R.
  // H1L/H2R: expected = 2*H - 1*(2*H) = 0. Match!
  //
  // But H1L/H2R is a regular compomer (not identity), already in the table
  // when charge_max >= 2. So this should work without identity compomers.

  double H_adduct = 1.00728;
  double M = 500.0;

  FeatureMap fm;
  Feature f1;
  f1.setMZ(M + H_adduct);                     // [M+H]+, charge 1
  f1.setRT(100.0);
  f1.setIntensity(1000.0f);
  f1.setCharge(1);
  fm.push_back(f1);

  Feature f2;
  f2.setMZ((2.0 * M + 2.0 * H_adduct) / 2.0);  // [2M+2H]2+, charge 2
  f2.setRT(100.0);
  f2.setIntensity(500.0f);
  f2.setCharge(2);
  fm.push_back(f2);

  MetaboliteFeatureDeconvolution mfd;
  Param p = mfd.getDefaults();
  p.setValue("charge_min", 1);
  p.setValue("charge_max", 2);
  p.setValue("charge_span_max", 2);
  p.setValue("max_multimer", 2);
  p.setValue("potential_adducts", std::vector<std::string>{"H:+:1.0"});
  p.setValue("mass_max_diff", 0.5);
  p.setValue("retention_max_diff", 2.0);
  p.setValue("retention_max_diff_local", 2.0);
  mfd.setParameters(p);

  FeatureMap fm_out;
  ConsensusMap cm_out, cm_out2;
  mfd.compute(fm, fm_out, cm_out, cm_out2);

  // Features should be grouped
  bool found_group = false;
  for (Size i = 0; i < cm_out.size(); ++i)
  {
    if (cm_out[i].size() >= 2) found_group = true;
  }
  TEST_EQUAL(found_group, true);
}
END_SECTION
```

**Step 2: Build and run**

Run: `/snap/cmake/1515/bin/cmake --build OpenMS-build --target MetaboliteFeatureDeconvolution_test -j$(nproc) && OpenMS-build/src/tests/class_tests/bin/MetaboliteFeatureDeconvolution_test`
Expected: PASS (or requires debugging — mixed charge with cross-multiplier is the most complex scenario)

**Step 3: Commit**

```bash
git add src/tests/class_tests/openms/source/MetaboliteFeatureDeconvolution_test.cpp
git commit -m "test: verify mixed-charge multimer detection [M+H]+ / [2M+2H]2+"
```

---

### Task 8: Run full regression suite

**Step 1: Run all decharging-related tests**

Run: `ctest --test-dir OpenMS-build -R "DECHARGING|Adduct|Compomer|ChargePair|MassExplainer" --output-on-failure`
Expected: All 12 tests PASS

**Step 2: Verify TOPP MetaboliteAdductDecharger output unchanged**

The TOPP tests compare output files. If they pass, existing behavior is unchanged.

Run: `ctest --test-dir OpenMS-build -R "TOPP_MetaboliteAdductDecharger" --output-on-failure`
Expected: All TOPP tests PASS (max_multimer defaults to 1, identity compomers not generated)

**Step 3: Final commit if any fixups needed**

```bash
git add -u && git commit -m "fix: address test failures from regression suite"
```
