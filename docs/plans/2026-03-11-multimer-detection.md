# Multimer Detection Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add mass-aware multimer detection (dimers, trimers) to MetaboliteAdductDecharger, controlled by `max_multimer` parameter (default 1 = disabled).

**Architecture:** Multiplier combinations (n1, n2) are tried in `candidateEdges_()` for each feature pair. Cross-multiplier matches use exact algebraic M cancellation: `n2×m1 - n1×m2 = n2×left_adduct_mass - n1×right_adduct_mass`. A new `queryMultimer()` on MassExplainer does linear scan within net_charge groups. Multiplier info is stored on ChargePair edges, not on Adduct.

**Tech Stack:** C++, OpenMS test framework, CMake

**Design doc:** `docs/plans/2026-03-11-multimer-detection-design.md`

**Build commands:**
- Build: `cmake --build OpenMS-build -j$(nproc)`
- Build single test: `cmake --build OpenMS-build --target <TestName> -j$(nproc)`
- Run single test: `ctest --test-dir OpenMS-build -R <TestName> -V`

---

### Task 1: Add `getSideMass()` to Compomer

**Files:**
- Modify: `src/openms/include/OpenMS/DATASTRUCTURES/Compomer.h:238` (before `getLabels`)
- Modify: `src/openms/source/DATASTRUCTURES/Compomer.cpp:291` (before `getLabels` impl)
- Modify: `src/tests/class_tests/openms/source/Compomer_test.cpp`

**Step 1: Write the failing test**

Add before the `END_TEST` line in `Compomer_test.cpp`:

```cpp
START_SECTION((double getSideMass(UInt side) const))
{
  // empty compomer has zero mass on both sides
  Compomer c;
  TEST_REAL_SIMILAR(c.getSideMass(Compomer::LEFT), 0.0);
  TEST_REAL_SIMILAR(c.getSideMass(Compomer::RIGHT), 0.0);

  // single adduct on right side
  Adduct a1(1, 2, 10.5, "Na", -0.3, 0);
  c.add(a1, Compomer::RIGHT);
  TEST_REAL_SIMILAR(c.getSideMass(Compomer::RIGHT), 2 * 10.5); // amount * singleMass
  TEST_REAL_SIMILAR(c.getSideMass(Compomer::LEFT), 0.0);

  // add adduct on left side
  Adduct a2(1, 3, 1.008, "H", -0.1, 0);
  c.add(a2, Compomer::LEFT);
  TEST_REAL_SIMILAR(c.getSideMass(Compomer::LEFT), 3 * 1.008);
  TEST_REAL_SIMILAR(c.getSideMass(Compomer::RIGHT), 2 * 10.5);

  // invalid side throws
  TEST_EXCEPTION(Exception::InvalidValue, c.getSideMass(Compomer::BOTH));
}
END_SECTION
```

**Step 2: Run test to verify it fails**

Run: `cmake --build OpenMS-build --target Compomer_test -j$(nproc) && ctest --test-dir OpenMS-build -R Compomer_test -V`
Expected: Build failure — `getSideMass` not declared

**Step 3: Write the declaration**

In `Compomer.h`, add after the `getLabels` declaration (around line 245):

```cpp
    /**
      @brief Get total adduct mass on a specific side

      Computes sum of amount * singleMass for all adducts on the given side.

      @param side Which side (LEFT or RIGHT)
      @return Total adduct mass contribution on that side
    */
    double getSideMass(const UInt side) const;
```

**Step 4: Write the implementation**

In `Compomer.cpp`, add after `getLabels` implementation:

```cpp
  double Compomer::getSideMass(const UInt side) const
  {
    if (side >= BOTH)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Compomer::getSideMass() only supports LEFT (0) or RIGHT (1), not: ", String(side));
    }
    double mass = 0.0;
    for (const auto& [formula, adduct] : cmp_[side])
    {
      mass += adduct.getAmount() * adduct.getSingleMass();
    }
    return mass;
  }
```

**Step 5: Run test to verify it passes**

Run: `cmake --build OpenMS-build --target Compomer_test -j$(nproc) && ctest --test-dir OpenMS-build -R Compomer_test -V`
Expected: PASS

**Step 6: Commit**

```bash
git add src/openms/include/OpenMS/DATASTRUCTURES/Compomer.h \
        src/openms/source/DATASTRUCTURES/Compomer.cpp \
        src/tests/class_tests/openms/source/Compomer_test.cpp
git commit -m "feat: add Compomer::getSideMass() for per-side adduct mass computation"
```

---

### Task 2: Add multiplier fields to ChargePair

**Files:**
- Modify: `src/openms/include/OpenMS/DATASTRUCTURES/ChargePair.h`
- Modify: `src/openms/source/DATASTRUCTURES/ChargePair.cpp`
- Modify: `src/tests/class_tests/openms/source/ChargePair_test.cpp`

**Step 1: Write the failing tests**

Add before the `END_TEST` line in `ChargePair_test.cpp`:

```cpp
START_SECTION((Int getMolMultiplier(UInt pairID) const))
{
  // default constructor gives multiplier 1
  ChargePair cp;
  TEST_EQUAL(cp.getMolMultiplier(0), 1);
  TEST_EQUAL(cp.getMolMultiplier(1), 1);
}
END_SECTION

START_SECTION((void setMolMultiplier(UInt pairID, Int m)))
{
  ChargePair cp;
  cp.setMolMultiplier(0, 2);
  cp.setMolMultiplier(1, 3);
  TEST_EQUAL(cp.getMolMultiplier(0), 2);
  TEST_EQUAL(cp.getMolMultiplier(1), 3);
}
END_SECTION

START_SECTION(([EXTRA] ChargePair multiplier in equality and hash))
{
  // multiplier affects equality
  ChargePair cp1(34, 45, 4, 5, cmp, 12.34, false);
  ChargePair cp2(34, 45, 4, 5, cmp, 12.34, false);
  cp1.setMolMultiplier(0, 2);
  TEST_FALSE(cp1 == cp2);

  // same multiplier -> equal
  cp2.setMolMultiplier(0, 2);
  TEST_TRUE(cp1 == cp2);

  // hash consistency
  std::hash<ChargePair> hasher;
  TEST_EQUAL(hasher(cp1), hasher(cp2));

  // different multiplier -> different hash (likely)
  cp2.setMolMultiplier(1, 3);
  TEST_FALSE(cp1 == cp2);
  TEST_NOT_EQUAL(hasher(cp1), hasher(cp2));
}
END_SECTION
```

**Step 2: Run test to verify it fails**

Run: `cmake --build OpenMS-build --target ChargePair_test -j$(nproc) && ctest --test-dir OpenMS-build -R ChargePair_test -V`
Expected: Build failure — `getMolMultiplier` not declared

**Step 3: Add declarations to ChargePair.h**

In `ChargePair.h`, add after `setActive` (line 97), before the `//@}` closing:

```cpp
    /// Returns the molecular multiplier (for element 0 or 1). Default is 1 (monomer).
    Int getMolMultiplier(UInt pairID) const;

    /// Set the molecular multiplier (for element 0 or 1). 1=monomer, 2=dimer, etc.
    void setMolMultiplier(UInt pairID, Int m);
```

Add two member fields after `is_active_` (line 124):

```cpp
    /// Molecular multiplier for first feature (1=monomer, 2=dimer, etc.)
    Int feature0_mol_multiplier_;
    /// Molecular multiplier for second feature
    Int feature1_mol_multiplier_;
```

**Step 4: Write the implementation**

In `ChargePair.cpp`, add `feature0_mol_multiplier_(1), feature1_mol_multiplier_(1)` to both constructors' initializer lists.

Default constructor (line 17-27): add before the closing brace of the initializer list:
```cpp
    feature0_mol_multiplier_(1),
    feature1_mol_multiplier_(1)
```

Parameterized constructor (line 30-46): add before the closing brace of the initializer list:
```cpp
    feature0_mol_multiplier_(1),
    feature1_mol_multiplier_(1)
```

Add the field to the assignment operator (after `is_active_` assignment, around line 65):
```cpp
    feature0_mol_multiplier_ = rhs.feature0_mol_multiplier_;
    feature1_mol_multiplier_ = rhs.feature1_mol_multiplier_;
```

Add the accessor implementations (after `setActive`, around line 171):
```cpp
  Int ChargePair::getMolMultiplier(UInt pairID) const
  {
    if (pairID == 0)
    {
      return feature0_mol_multiplier_;
    }
    else
    {
      return feature1_mol_multiplier_;
    }
  }

  void ChargePair::setMolMultiplier(UInt pairID, Int m)
  {
    if (pairID == 0)
    {
      feature0_mol_multiplier_ = m;
    }
    else
    {
      feature1_mol_multiplier_ = m;
    }
  }
```

Add to `operator==` (line 176-185) — add two more comparisons:
```cpp
           (feature0_mol_multiplier_ == i.feature0_mol_multiplier_) &&
           (feature1_mol_multiplier_ == i.feature1_mol_multiplier_);
```

Add to `operator<<` (line 193-201):
```cpp
       << "MolMultiplier: " << cp.getMolMultiplier(0) << " : " << cp.getMolMultiplier(1) << "\n";
```

**Step 5: Update hash function**

In `ChargePair.h`, add to the hash function (after the `is_active` hash, line 149):
```cpp
      OpenMS::hash_combine(seed, OpenMS::hash_int(cp.getMolMultiplier(0)));
      OpenMS::hash_combine(seed, OpenMS::hash_int(cp.getMolMultiplier(1)));
```

**Step 6: Run test to verify it passes**

Run: `cmake --build OpenMS-build --target ChargePair_test -j$(nproc) && ctest --test-dir OpenMS-build -R ChargePair_test -V`
Expected: PASS

**Step 7: Commit**

```bash
git add src/openms/include/OpenMS/DATASTRUCTURES/ChargePair.h \
        src/openms/source/DATASTRUCTURES/ChargePair.cpp \
        src/tests/class_tests/openms/source/ChargePair_test.cpp
git commit -m "feat: add mol_multiplier fields to ChargePair for multimer edges"
```

---

### Task 3: Add `toAdductString` with multiplier parameter

**Files:**
- Modify: `src/openms/include/OpenMS/DATASTRUCTURES/Adduct.h:71`
- Modify: `src/openms/source/DATASTRUCTURES/Adduct.cpp:151`
- Modify: `src/tests/class_tests/openms/source/Adduct_test.cpp`

**Step 1: Write the failing test**

Add before `END_TEST` in `Adduct_test.cpp`:

```cpp
START_SECTION((static String toAdductString(const String& ion_string, const Int& charge, Int mol_multiplier)))
{
  // monomer (multiplier=1) — no prefix
  String r1 = Adduct::toAdductString("H1", 1, 1);
  TEST_EQUAL(r1, "[M+H]+");

  // dimer
  String r2 = Adduct::toAdductString("H1", 1, 2);
  TEST_EQUAL(r2, "[2M+H]+");

  // trimer with Na
  String r3 = Adduct::toAdductString("Na1", 1, 3);
  TEST_EQUAL(r3, "[3M+Na]+");

  // dimer negative mode
  String r4 = Adduct::toAdductString("H-1", -1, 2);
  TEST_EQUAL(r4, "[2M-H]-");
}
END_SECTION
```

**Step 2: Run test to verify it fails**

Run: `cmake --build OpenMS-build --target Adduct_test -j$(nproc) && ctest --test-dir OpenMS-build -R Adduct_test -V`
Expected: Build failure — static overload not declared

**Step 3: Add declaration**

In `Adduct.h`, after line 71 (`toAdductString`), add:

```cpp
    /// Convert to adduct string with explicit multiplier (e.g., mol_multiplier=2 -> "[2M+H]+")
    static String toAdductString(const String& ion_string, const Int& charge, Int mol_multiplier);
```

**Step 4: Write the implementation**

In `Adduct.cpp`, after the existing `toAdductString` (line 176), add:

```cpp
  String Adduct::toAdductString(const String& ion_string, const Int& charge, Int mol_multiplier)
  {
    EmpiricalFormula ef(ion_string);
    String charge_sign = charge >= 0 ? "+" : "-";
    String s("[");

    if (mol_multiplier > 1)
    {
      s += String(mol_multiplier);
    }
    s += "M";

    // elements sorted canonically (by string)
    std::map<String, String> sorted_elem_map;
    for (const auto& element_count : ef)
    {
      String e_symbol(element_count.first->getSymbol());
      String tmp = element_count.second > 0 ? "+" : "-";
      tmp += std::abs(element_count.second) > 1 ? String(std::abs(element_count.second)) : "";
      tmp += e_symbol;
      sorted_elem_map[e_symbol] = std::move(tmp);
    }
    for (const auto& sorted_e_cnt : sorted_elem_map)
    {
      s += sorted_e_cnt.second;
    }
    s += String("]");
    s += std::abs(charge) > 1 ? String(std::abs(charge)) : "";
    s += charge_sign;

    return s;
  }
```

**Step 5: Run test to verify it passes**

Run: `cmake --build OpenMS-build --target Adduct_test -j$(nproc) && ctest --test-dir OpenMS-build -R Adduct_test -V`
Expected: PASS

**Step 6: Commit**

```bash
git add src/openms/include/OpenMS/DATASTRUCTURES/Adduct.h \
        src/openms/source/DATASTRUCTURES/Adduct.cpp \
        src/tests/class_tests/openms/source/Adduct_test.cpp
git commit -m "feat: add static Adduct::toAdductString with multiplier for multimer annotation"
```

---

### Task 4: Add `queryMultimer()` to MassExplainer

**Files:**
- Modify: `src/openms/include/OpenMS/DATASTRUCTURES/MassExplainer.h:177`
- Modify: `src/openms/source/DATASTRUCTURES/MassExplainer.cpp:326`
- Modify: `src/tests/class_tests/openms/source/MassExplainer_test.cpp`

**Step 1: Write the failing test**

Add before `END_TEST` in `MassExplainer_test.cpp`:

```cpp
START_SECTION((SignedSize queryMultimer(const Int net_charge, const double m1, const double m2, const Int n1, const Int n2, const double tolerance, const float thresh_log_p, std::vector<std::vector<Compomer>::const_iterator>& hits) const))
{
  // Build a MassExplainer with H+ adduct only, simple charge range
  Adduct h_plus(1, 1, 1.00728, "H", log(0.9), 0);
  MassExplainer::AdductsType adducts;
  adducts.push_back(h_plus);
  MassExplainer me(adducts, 1, 3, 3, log(1e-10), 0);
  me.compute();

  // Simulate: feature1 = [M+H]+ at m/z=100.00728 (M=99, q1=1, n1=1)
  //           feature2 = [2M+H]+ at m/z=199.00728 (M=99, q2=1, n2=2)
  // m1 = mz1 * |q1| = 100.00728
  // m2 = mz2 * |q2| = 199.00728
  // Cross-mult check: n2*m1 - n1*m2 = 2*100.00728 - 1*199.00728 = 1.00728
  // This should match a compomer with: n2*left_mass - n1*right_mass = 2*H_mass - 1*H_mass = 1.00728
  double m1_val = 100.00728;
  double m2_val = 199.00728;

  std::vector<MassExplainer::CompomerIterator> multimer_hits;
  SignedSize n_hits = me.queryMultimer(0, m1_val, m2_val, 1, 2, 0.01, log(1e-10), multimer_hits);
  TEST_EQUAL(n_hits > 0, true);

  // No match with wrong multiplier combo
  multimer_hits.clear();
  n_hits = me.queryMultimer(0, m1_val, m2_val, 1, 3, 0.01, log(1e-10), multimer_hits);
  TEST_EQUAL(n_hits, 0);
}
END_SECTION
```

**Step 2: Run test to verify it fails**

Run: `cmake --build OpenMS-build --target MassExplainer_test -j$(nproc) && ctest --test-dir OpenMS-build -R MassExplainer_test -V`
Expected: Build failure — `queryMultimer` not declared

**Step 3: Add declaration**

In `MassExplainer.h`, after the `query()` declaration (line 177), add:

```cpp

    /**
      @brief Search for cross-multiplier explanations between two features

      For two features with different molecular multipliers (n1, n2), finds compomers
      where M cancels exactly: n2*(m1 - left_mass) = n1*(m2 - right_mass).

      Uses linear scan within the net_charge group (binary search to find the group).

      @param[in] net_charge Net charge difference (q2 - q1)
      @param[in] m1 Charge-scaled mass of feature 1 (mz1 * |q1|)
      @param[in] m2 Charge-scaled mass of feature 2 (mz2 * |q2|)
      @param[in] n1 Molecular multiplier for feature 1
      @param[in] n2 Molecular multiplier for feature 2 (must differ from n1)
      @param[in] tolerance Absolute mass tolerance
      @param[in] thresh_log_p Minimum log probability
      @param[out] hits Vector of iterators to matching compomers
      @return Number of matching compomers found
    */
    SignedSize queryMultimer(const Int net_charge,
                             const double m1,
                             const double m2,
                             const Int n1,
                             const Int n2,
                             const double tolerance,
                             const float thresh_log_p,
                             std::vector<CompomerIterator>& hits) const;
```

**Step 4: Write the implementation**

In `MassExplainer.cpp`, after the `query()` implementation (line 326), add:

```cpp
  SignedSize MassExplainer::queryMultimer(const Int net_charge,
                                          const double m1,
                                          const double m2,
                                          const Int n1,
                                          const Int n2,
                                          const double tolerance,
                                          const float thresh_log_p,
                                          std::vector<CompomerIterator>& hits) const
  {
    hits.clear();

    // Find the net_charge group using binary search on the sorted explanations_
    // Use a wide mass range to capture the entire net_charge group
    Compomer cmp_low(net_charge, -std::numeric_limits<double>::max(), 1);
    auto group_begin = std::lower_bound(explanations_.begin(), explanations_.end(), cmp_low);

    Compomer cmp_high(net_charge + 1, -std::numeric_limits<double>::max(), 1);
    auto group_end = std::lower_bound(explanations_.begin(), explanations_.end(), cmp_high);

    // The observed value: n2*m1 - n1*m2
    double observed = static_cast<double>(n2) * m1 - static_cast<double>(n1) * m2;

    // Linear scan within net_charge group
    for (auto it = group_begin; it != group_end; ++it)
    {
      if (it->getLogP() < thresh_log_p)
      {
        continue;
      }

      // Compute expected: n2*left_mass - n1*right_mass
      double left_mass = it->getSideMass(Compomer::LEFT);
      double right_mass = it->getSideMass(Compomer::RIGHT);
      double expected = static_cast<double>(n2) * left_mass - static_cast<double>(n1) * right_mass;

      if (fabs(observed - expected) <= tolerance)
      {
        hits.push_back(it);
      }
    }

    return static_cast<SignedSize>(hits.size());
  }
```

**Step 5: Run test to verify it passes**

Run: `cmake --build OpenMS-build --target MassExplainer_test -j$(nproc) && ctest --test-dir OpenMS-build -R MassExplainer_test -V`
Expected: PASS

**Step 6: Commit**

```bash
git add src/openms/include/OpenMS/DATASTRUCTURES/MassExplainer.h \
        src/openms/source/DATASTRUCTURES/MassExplainer.cpp \
        src/tests/class_tests/openms/source/MassExplainer_test.cpp
git commit -m "feat: add MassExplainer::queryMultimer() for cross-multiplier matching"
```

---

### Task 5: Add `max_multimer` and `multimer_log_penalty` parameters

**Files:**
- Modify: `src/openms/source/ANALYSIS/DECHARGING/MetaboliteFeatureDeconvolution.cpp:108`

**Step 1: Add parameter declarations**

In `MetaboliteFeatureDeconvolution.cpp`, after the `max_neutrals` parameter (line 108), add:

```cpp
    defaults_.setValue("max_multimer", 1, "Maximum molecular multiplier for multimer detection (e.g., 2 enables dimer [2M+...] detection, 3 enables trimers). Default 1 disables multimer detection.");
    defaults_.setMinInt("max_multimer", 1);
    defaults_.setValue("multimer_log_penalty", -2.0, "Log-probability penalty per multiplier step for cross-multiplier edges. Applied as (max(n1,n2)-1) * penalty. More negative values make the solver prefer monomer explanations more strongly.");
```

**Step 2: Run existing tests to verify nothing breaks**

Run: `cmake --build OpenMS-build --target MetaboliteFeatureDeconvolution_test -j$(nproc) && ctest --test-dir OpenMS-build -R MetaboliteFeatureDeconvolution_test -V`
Expected: PASS (new parameters have defaults, existing behavior unchanged)

**Step 3: Commit**

```bash
git add src/openms/source/ANALYSIS/DECHARGING/MetaboliteFeatureDeconvolution.cpp
git commit -m "feat: add max_multimer and multimer_log_penalty parameters"
```

---

### Task 6: Wire multimer detection into candidateEdges_()

**Files:**
- Modify: `src/openms/source/ANALYSIS/DECHARGING/MetaboliteFeatureDeconvolution.cpp`

This is the core algorithmic change. Modify `candidateEdges_()` to loop over multiplier combinations and call `queryMultimer()` for cross-multiplier matches.

**Step 1: Read and store new parameters**

At the top of `candidateEdges_()` (around line 361, after `max_neutrals`), add:

```cpp
    Int max_multimer = param_.getValue("max_multimer");
    double multimer_log_penalty = param_.getValue("multimer_log_penalty");
```

**Step 2: Add multiplier loop around the MassExplainer query**

The existing query at line 505 currently does:
```cpp
hits = me.query(q2 - q1, naive_mass_diff, abs_mass_diff, thresh_logp, md_s, md_e);
```

Replace the block from `naive_mass_diff` computation (line 478) through the end of the hit processing loop (line 619) with the multiplier-aware version. The key structure is:

```cpp
            CoordinateType m2 = mz2 * abs(q2);

            double abs_mass_diff;
            // ... existing tolerance computation (lines 480-502) stays the same ...

            // Loop over multiplier combinations
            for (Int n1 = 1; n1 <= max_multimer; ++n1)
            {
              for (Int n2 = n1; n2 <= max_multimer; ++n2)
              {
                if (n1 == n2)
                {
                  // Same-multiplier: existing path, M cancels
                  CoordinateType naive_mass_diff = m2 - m1;
                  hits = me.query(q2 - q1, naive_mass_diff, abs_mass_diff, thresh_logp, md_s, md_e);
                  OPENMS_PRECONDITION(hits >= 0, "MetaboliteFeatureDeconvolution querying #hits got negative result!");
                  overallHits += hits;

                  if (hits > 0)
                  {
                    // ... existing hit processing (lines 513-619) ...
                    // When creating ChargePair (line 606), set multipliers:
                    // cp.setMolMultiplier(0, n1);
                    // cp.setMolMultiplier(1, n2);
                  }
                }
                else
                {
                  // Cross-multiplier: exact algebraic match
                  // Tolerance scales: n2*tol1 + n1*tol2
                  double cross_tol = n2 * (mz_diff_max_abs_q1) + n1 * (mz_diff_max_abs_q2);
                  // where mz_diff_max_abs_q1/q2 are the per-feature tolerance contributions

                  std::vector<MassExplainer::CompomerIterator> multimer_hits;
                  hits = me.queryMultimer(q2 - q1, m1, m2, n1, n2, cross_tol, thresh_logp, multimer_hits);
                  overallHits += hits;

                  // Process multimer hits similarly to same-multiplier hits
                  // Apply multimer_log_penalty: (max(n1,n2) - 1) * multimer_log_penalty
                  // For each hit, process using the same filtering logic
                  // Create ChargePair with cp.setMolMultiplier(0, n1); cp.setMolMultiplier(1, n2);

                  // Also try (n2, n1) — feature 0 as the higher multiplier
                  hits = me.queryMultimer(q2 - q1, m1, m2, n2, n1, cross_tol, thresh_logp, multimer_hits);
                  // Process with cp.setMolMultiplier(0, n2); cp.setMolMultiplier(1, n1);
                }
              }
            }
```

**Important implementation details:**
- The per-feature tolerance contributions need to be computed before the multiplier loop. Currently `abs_mass_diff` combines both. Split into: `double tol_q1 = mz_diff_max * abs(q1)` (Da mode) or `mz1 * mz_diff_max * 1e-6 * abs(q1)` (ppm mode), and similarly `tol_q2`. Then `abs_mass_diff = tol_q1 + tol_q2` for same-multiplier, and `cross_tol = n2 * tol_q1 + n1 * tol_q2` for cross-multiplier.
- The hit processing logic (lines 513-619) should be extracted into a helper or the cross-multiplier path should replicate the essential filtering (RT check, charge check, default adduct fill, ChargePair creation).
- For cross-multiplier hits, apply the penalty: when a hit passes all filters, adjust the compomer's log_p before creating the ChargePair: `cmp.setLogP(cmp.getLogP() + (std::max(n1, n2) - 1) * multimer_log_penalty)`.

**Step 3: Run existing tests**

Run: `cmake --build OpenMS-build --target MetaboliteFeatureDeconvolution_test -j$(nproc) && ctest --test-dir OpenMS-build -R MetaboliteFeatureDeconvolution_test -V`
Expected: PASS (max_multimer defaults to 1, so only n1=n2=1 iteration runs — identical to old behavior)

**Step 4: Commit**

```bash
git add src/openms/source/ANALYSIS/DECHARGING/MetaboliteFeatureDeconvolution.cpp
git commit -m "feat: wire multimer detection into candidateEdges_() with queryMultimer"
```

---

### Task 7: Update annotation to use multiplier from ChargePair

**Files:**
- Modify: `src/openms/source/ANALYSIS/DECHARGING/MetaboliteFeatureDeconvolution.cpp:321-352` (annotate_feature_)
- Modify: `src/openms/source/ANALYSIS/DECHARGING/MetaboliteFeatureDeconvolution.cpp:828-832` (call site)

**Step 1: Modify annotate_feature_ signature**

Change the signature (line 321) to accept a multiplier:

```cpp
void MetaboliteFeatureDeconvolution::annotate_feature_(FeatureMap& fm_out, Adduct& default_adduct,
    Compomer& c, const Size f_idx, const UInt comp_side, const Int new_q, const Int old_q,
    const Int mol_multiplier)
```

Update the header declaration accordingly in `MetaboliteFeatureDeconvolution.h`.

**Step 2: Use multiplier in annotation**

In `annotate_feature_`, change line 336 from:
```cpp
StringList dc_new_adducts = ListUtils::create<String>(adduct.toAdductString(ef_.toString(), new_q));
```
to:
```cpp
StringList dc_new_adducts = ListUtils::create<String>(Adduct::toAdductString(ef_.toString(), new_q, mol_multiplier));
```

After line 337, add:
```cpp
    if (mol_multiplier > 1)
    {
      fm_out[f_idx].setMetaValue("mol_multiplier", mol_multiplier);
    }
```

**Step 3: Update call sites**

At lines 830 and 832, change to pass multiplier from ChargePair:

```cpp
        Int mult0 = feature_relation[i].getMolMultiplier(0);
        Int mult1 = feature_relation[i].getMolMultiplier(1);
        // - left
        annotate_feature_(fm_out, default_adduct, c, f0_idx, Compomer::LEFT, new_q0, old_q0, mult0);
        // - right
        annotate_feature_(fm_out, default_adduct, c, f1_idx, Compomer::RIGHT, new_q1, old_q1, mult1);
```

**Step 4: Run existing tests**

Run: `cmake --build OpenMS-build --target MetaboliteFeatureDeconvolution_test -j$(nproc) && ctest --test-dir OpenMS-build -R MetaboliteFeatureDeconvolution_test -V`
Expected: PASS (all multipliers default to 1, annotation unchanged)

**Step 5: Commit**

```bash
git add src/openms/source/ANALYSIS/DECHARGING/MetaboliteFeatureDeconvolution.cpp \
        src/openms/include/OpenMS/ANALYSIS/DECHARGING/MetaboliteFeatureDeconvolution.h
git commit -m "feat: annotate features with multiplier from ChargePair for multimer labels"
```

---

### Task 8: Integration test — multimer detection end-to-end

**Files:**
- Modify: `src/tests/class_tests/openms/source/MetaboliteFeatureDeconvolution_test.cpp`

**Step 1: Write the integration test**

Add a new test section that creates a FeatureMap with a monomer and dimer feature, runs decharging with `max_multimer=2`, and verifies correct grouping and annotation.

```cpp
START_SECTION(([EXTRA] Multimer detection with max_multimer=2))
{
  // Create features:
  // Feature 1: [M+H]+ at m/z = 100.00728 (M ≈ 99.0, charge 1)
  // Feature 2: [2M+H]+ at m/z = 199.00728 (M ≈ 99.0, charge 1)
  // Both at same RT so they can be paired
  FeatureMap fm;
  Feature f1;
  f1.setMZ(100.00728);
  f1.setRT(100.0);
  f1.setIntensity(1000.0);
  f1.setCharge(1);
  // Add minimal convex hull for RT overlap check
  std::vector<DPosition<2>> hull_points1 = {{99.0, 99.0}, {99.0, 101.0}, {101.0, 101.0}, {101.0, 99.0}};
  ConvexHull2D hull1;
  hull1.addPoints(hull_points1);
  f1.getConvexHulls().push_back(hull1);
  fm.push_back(f1);

  Feature f2;
  f2.setMZ(199.00728);
  f2.setRT(100.0);
  f2.setIntensity(500.0);
  f2.setCharge(1);
  std::vector<DPosition<2>> hull_points2 = {{99.0, 198.0}, {99.0, 200.0}, {101.0, 200.0}, {101.0, 198.0}};
  ConvexHull2D hull2;
  hull2.addPoints(hull_points2);
  f2.getConvexHulls().push_back(hull2);
  fm.push_back(f2);

  // Run decharger with multimer detection enabled
  MetaboliteFeatureDeconvolution mfd;
  Param p = mfd.getDefaults();
  p.setValue("charge_min", 1);
  p.setValue("charge_max", 1);
  p.setValue("charge_span_max", 1);
  p.setValue("max_multimer", 2);
  p.setValue("potential_adducts", std::vector<std::string>{"H:+:0.9"});
  p.setValue("mass_max_diff", 0.01);
  p.setValue("retention_max_diff", 2.0);
  mfd.setParameters(p);

  FeatureMap fm_out;
  ConsensusMap cm_out, cm_out2;
  mfd.compute(fm, fm_out, cm_out, cm_out2);

  // Check that at least one consensus feature groups the two features
  bool found_multimer_group = false;
  for (Size i = 0; i < cm_out.size(); ++i)
  {
    if (cm_out[i].size() >= 2)
    {
      found_multimer_group = true;
    }
  }
  TEST_EQUAL(found_multimer_group, true);

  // Check annotation: one feature should have mol_multiplier=2
  bool found_dimer_annotation = false;
  for (Size i = 0; i < fm_out.size(); ++i)
  {
    if (fm_out[i].metaValueExists("mol_multiplier"))
    {
      TEST_EQUAL((int)fm_out[i].getMetaValue("mol_multiplier"), 2);
      found_dimer_annotation = true;
    }
  }
  TEST_EQUAL(found_dimer_annotation, true);

  // Regression: max_multimer=1 should NOT find multimer relationships
  MetaboliteFeatureDeconvolution mfd2;
  Param p2 = mfd2.getDefaults();
  p2.setValue("charge_min", 1);
  p2.setValue("charge_max", 1);
  p2.setValue("charge_span_max", 1);
  p2.setValue("max_multimer", 1);
  p2.setValue("potential_adducts", std::vector<std::string>{"H:+:0.9"});
  p2.setValue("mass_max_diff", 0.01);
  p2.setValue("retention_max_diff", 2.0);
  mfd2.setParameters(p2);

  FeatureMap fm_out2_feat;
  ConsensusMap cm_out2_a, cm_out2_b;
  mfd2.compute(fm, fm_out2_feat, cm_out2_a, cm_out2_b);

  // With max_multimer=1, same-charge features with large mass diff should not be grouped
  bool found_group_no_multimer = false;
  for (Size i = 0; i < cm_out2_a.size(); ++i)
  {
    if (cm_out2_a[i].size() >= 2)
    {
      found_group_no_multimer = true;
    }
  }
  TEST_EQUAL(found_group_no_multimer, false);
}
END_SECTION
```

**Step 2: Build and run**

Run: `cmake --build OpenMS-build --target MetaboliteFeatureDeconvolution_test -j$(nproc) && ctest --test-dir OpenMS-build -R MetaboliteFeatureDeconvolution_test -V`
Expected: PASS

**Step 3: Commit**

```bash
git add src/tests/class_tests/openms/source/MetaboliteFeatureDeconvolution_test.cpp
git commit -m "test: add integration test for multimer detection with max_multimer=2"
```

---

### Task 9: Run full test suite and verify no regressions

**Files:** None (verification only)

**Step 1: Build everything**

Run: `cmake --build OpenMS-build -j$(nproc)`
Expected: Clean build, no warnings in changed files

**Step 2: Run all decharging-related tests**

Run: `ctest --test-dir OpenMS-build -R "Adduct_test|Compomer_test|ChargePair_test|MassExplainer_test|MetaboliteFeatureDeconvolution_test|FeatureDeconvolution_test" -V`
Expected: All PASS

**Step 3: Commit if any final fixes needed, otherwise done**

```bash
# Final commit if needed
git add -u
git commit -m "fix: address test failures from multimer detection integration"
```
