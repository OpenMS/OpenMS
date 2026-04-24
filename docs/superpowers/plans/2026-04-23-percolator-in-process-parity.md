# Percolator in-process parity tests — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

> **Follow-up migration (done 2026-04-24):** the vendored percolator tree was bumped from `rel-3-08-01` to upstream eb157f7 in [`2026-04-24-percolator-eb157f7-migration.md`](2026-04-24-percolator-eb157f7-migration.md). All tolerances below remained valid — observed values bit-identical to pre-migration on the 5 gated sections. `pep_method="nonparametric"` path was confirmed algorithmically unchanged by eb157f7.

**Goal:** Add three pillars of parity tests — data-diversity extensions at the library layer, a new TOPP-adapter-layer parity test, and four reproducibility invariants — that give high confidence the in-process Percolator matches the bundled `percolator` executable on workloads OpenMS users actually run.

**Architecture:** Pure test additions — no production code is modified. Three files are touched: extend `Percolator_subprocess_parity_test.cpp` (P1.a/b/c), extend `Percolator_test.cpp` (P3.a/b/c/d), and create `PercolatorAdapter_parity_test.cpp` (P2). Subprocess-dependent tests skip silently when `PERCOLATOR_BINARY` (or `PERCOLATOR_ADAPTER_BINARY` for P2) is unset.

**Tech Stack:** OpenMS class-test framework (`START_TEST` / `START_SECTION` / `TEST_EQUAL` / `TEST_TRUE`), C++17, CMake for test registration and env-var plumbing, `std::system` for subprocess invocation, `IdXMLFile` / `PercolatorInfile` / `TextFile` for file I/O.

## Spec reference

See `docs/superpowers/specs/2026-04-23-percolator-in-process-parity-design.md` for the full design. Tolerances in this plan mirror the tolerance table there; deviations (if any) are flagged inline.

## TDD flow for test-only additions

Standard red-green-refactor doesn't cleanly apply when the deliverable *is* the test. For every task below the flow is:

1. Write the test section.
2. Build the test binary.
3. Run it. If parity holds, the test passes on the first run — commit.
4. If the test unexpectedly fails, investigate: either parity is genuinely broken (real finding — stop and flag to the user) or the test has a bug (fix and re-run).

Optional meta-check at the end of the series: revert one line of `Percolator.cpp` (e.g., clear the seed wiring) and confirm the new tests catch it. Not required by this plan.

## File structure

| File | Action | Purpose |
|---|---|---|
| `src/tests/class_tests/openms/source/Percolator_test.cpp` | Modify | Add P3.a/b/c/d reproducibility sections; tighten existing `[EXTRA] reproducibility with fixed seed` (P3.b). |
| `src/tests/class_tests/openms/source/Percolator_subprocess_parity_test.cpp` | Modify | Add P1.a (realistic idXML), P1.b (larger synthetic), P1.c (multi-file PIN) sections. Extend `generateSyntheticData` and `writePinFile` helpers. |
| `src/tests/class_tests/openms/source/PercolatorAdapter_parity_test.cpp` | Create | New P2 test: invoke `PercolatorAdapter` twice with `-use_subprocess true/false` on `PercolatorAdapter_1.idXML`, compare output idXMLs. |
| `src/tests/class_tests/openms/executables.cmake` | Modify | Register `PercolatorAdapter_parity_test` so the class-test framework builds and adds it to CTest. |
| `src/tests/class_tests/openms/CMakeLists.txt` | Modify | Pass `PERCOLATOR_BINARY` and `PERCOLATOR_ADAPTER_BINARY=$<TARGET_FILE:PercolatorAdapter>` env vars to `PercolatorAdapter_parity_test`. |

---

## Task 1: P3.a — Same-instance, consecutive calls, bit-identical

**Files:**
- Modify: `src/tests/class_tests/openms/source/Percolator_test.cpp` — append new `START_SECTION` block at the end, before `END_TEST`.

**What this guards:** a regression where calling `rescore()` twice on the same `Percolator` instance produces different results due to persistent RNG state, leaked global statics in the vendored code, or accumulated scratch buffers.

- [ ] **Step 1: Append new test section**

Locate the line `END_TEST` at the bottom of `src/tests/class_tests/openms/source/Percolator_test.cpp` and insert the following section immediately before it:

```cpp
START_SECTION([EXTRA] same-instance consecutive calls are bit-identical)
{
  // Guards against persistent-RNG / leaked-global-state regressions: a
  // second rescore() on the same instance must produce exactly the same
  // output as the first.
  RescoreInput input = makeModeratelySeparableInput_(500, 123);

  Percolator p;
  Param par = p.getDefaults();
  par.setValue("seed", 123);
  p.setParameters(par);

  RescoreOutput out1 = p.rescore(input);
  RescoreOutput out2 = p.rescore(input);

  TEST_EQUAL(out1.scores.size(), out2.scores.size())
  TEST_EQUAL(out1.q_values.size(), out2.q_values.size())
  TEST_EQUAL(out1.peps.size(), out2.peps.size())

  bool bit_identical = out1.scores.size() == out2.scores.size();
  for (size_t i = 0; bit_identical && i < out1.scores.size(); ++i)
  {
    if (out1.scores[i]   != out2.scores[i])   bit_identical = false;
    if (out1.q_values[i] != out2.q_values[i]) bit_identical = false;
    if (out1.peps[i]     != out2.peps[i])     bit_identical = false;
  }
  TEST_EQUAL(bit_identical, true)
}
END_SECTION
```

- [ ] **Step 2: Build**

```bash
cmake --build OpenMS-build --target Percolator_test -j$(nproc)
```

Expected: clean build, no compiler errors.

- [ ] **Step 3: Run the test**

```bash
ctest --test-dir OpenMS-build -R '^Percolator_test$' --output-on-failure
```

Expected: `Percolator_test` passes, output lists one additional `[EXTRA]` section.

- [ ] **Step 4: Commit**

```bash
git add src/tests/class_tests/openms/source/Percolator_test.cpp
git commit -m "$(cat <<'EOF'
test(Percolator): same-instance consecutive calls bit-identical

New [EXTRA] section on Percolator_test: two rescore() calls on the same
instance must produce bit-identical scores/q-values/PEPs. Guards against
persistent-RNG or leaked-global-state regressions.
EOF
)"
```

---

## Task 2: P3.b — Tighten existing cross-instance same-seed section to bit-identical

**Files:**
- Modify: `src/tests/class_tests/openms/source/Percolator_test.cpp` — replace the existing `[EXTRA] reproducibility with fixed seed` section body (currently at lines 120-155 before this task lands; locate by its start-section string).

**What this guards:** the existing test used `1e-9` slack that could hide small drift. Tightening to strict `==` and extending to q-values/PEPs surfaces any residual non-determinism.

- [ ] **Step 1: Replace the existing section body**

Find this block in `Percolator_test.cpp`:

```cpp
START_SECTION([EXTRA] reproducibility with fixed seed)
{
  // ... current body uses 1e-9 tolerance on scores only ...
}
END_SECTION
```

Replace the whole block (from `START_SECTION` to `END_SECTION` inclusive) with:

```cpp
START_SECTION([EXTRA] reproducibility with fixed seed is bit-identical)
{
  // Cross-instance, same-seed invariance. Tightened from 1e-9 tolerance on
  // scores only to strict equality on scores, q-values, and PEPs.
  RescoreInput input;
  input.feature_names = StringList{"f"};
  const size_t n = 400;
  input.features.reserve(n);
  input.is_decoy.reserve(n);
  for (size_t i = 0; i < n; ++i)
  {
    const bool is_decoy = (i % 2 == 1);
    input.features.push_back({(is_decoy ? -1.0 : +1.0) + 0.01 * static_cast<double>(i)});
    input.is_decoy.push_back(is_decoy);
  }

  Percolator p1;
  Percolator p2;
  for (auto* p : {&p1, &p2})
  {
    Param par = p->getDefaults();
    par.setValue("seed", 17);
    p->setParameters(par);
  }

  RescoreOutput out1 = p1.rescore(input);
  RescoreOutput out2 = p2.rescore(input);

  TEST_EQUAL(out1.scores.size(), out2.scores.size())
  TEST_EQUAL(out1.q_values.size(), out2.q_values.size())
  TEST_EQUAL(out1.peps.size(), out2.peps.size())

  bool bit_identical = out1.scores.size() == out2.scores.size();
  for (size_t i = 0; bit_identical && i < out1.scores.size(); ++i)
  {
    if (out1.scores[i]   != out2.scores[i])   bit_identical = false;
    if (out1.q_values[i] != out2.q_values[i]) bit_identical = false;
    if (out1.peps[i]     != out2.peps[i])     bit_identical = false;
  }
  TEST_EQUAL(bit_identical, true)
}
END_SECTION
```

- [ ] **Step 2: Build**

```bash
cmake --build OpenMS-build --target Percolator_test -j$(nproc)
```

Expected: clean build.

- [ ] **Step 3: Run**

```bash
ctest --test-dir OpenMS-build -R '^Percolator_test$' --output-on-failure
```

Expected: all sections pass, including the tightened one. If this fails, the vendored code has cross-instance drift that was hiding under `1e-9` — stop and flag.

- [ ] **Step 4: Commit**

```bash
git add src/tests/class_tests/openms/source/Percolator_test.cpp
git commit -m "$(cat <<'EOF'
test(Percolator): tighten cross-instance same-seed test to bit-identical

Upgrades [EXTRA] reproducibility with fixed seed from 1e-9 tolerance on
scores only to strict equality on scores/q-values/PEPs. If residual drift
exists across instances, it now surfaces instead of being hidden.
EOF
)"
```

---

## Task 3: P3.c — Input-order invariance (Jaccard characterization)

**Files:**
- Modify: `src/tests/class_tests/openms/source/Percolator_test.cpp` — append new section before `END_TEST`.
- Headers already included in this file: `<algorithm>`, `<numeric>`. We also need `<set>` and `<random>`. Add them to the top of the file if missing.

**What this guards:** the docstring says results depend on input ordering; this section documents the magnitude of that dependence so an unexpected increase is caught.

- [ ] **Step 1: Ensure required headers are present**

Near the top of `Percolator_test.cpp` (after the existing `#include`s), ensure these are present — add any that are missing:

```cpp
#include <random>
#include <set>
#include <iterator>
```

- [ ] **Step 2: Append new section**

Insert this section immediately before `END_TEST`:

```cpp
START_SECTION([EXTRA] input-order invariance Jaccard characterization)
{
  // Documents the input-order dependence the docstring calls out. Asserts
  // Pearson r >= 0.99 on scores (after un-permuting) and Jaccard >= 0.95
  // on accepted-target sets at q <= 0.01. A regression that increases the
  // drift (e.g., Jaccard dropping from 0.95 to 0.60) is caught.
  RescoreInput base = makeModeratelySeparableInput_(800, 456);

  std::vector<size_t> perm(base.features.size());
  std::iota(perm.begin(), perm.end(), 0);
  std::mt19937 rng(456);
  std::vector<size_t> perm_shuf = perm;
  std::shuffle(perm_shuf.begin(), perm_shuf.end(), rng);

  auto permuteInput = [](const RescoreInput& in, const std::vector<size_t>& p) {
    RescoreInput out;
    out.feature_names = in.feature_names;
    out.features.reserve(p.size());
    out.is_decoy.reserve(p.size());
    for (size_t i : p) {
      out.features.push_back(in.features[i]);
      out.is_decoy.push_back(in.is_decoy[i]);
    }
    return out;
  };

  Percolator p1, p2;
  for (auto* p : {&p1, &p2})
  {
    Param par = p->getDefaults();
    par.setValue("seed", 456);
    p->setParameters(par);
  }

  RescoreOutput out1 = p1.rescore(permuteInput(base, perm));
  RescoreOutput out2 = p2.rescore(permuteInput(base, perm_shuf));

  // Un-permute out2 so row i refers to base row i.
  const size_t n = out1.scores.size();
  std::vector<double> scores2_aligned(n), qvals2_aligned(n);
  for (size_t i = 0; i < perm_shuf.size(); ++i)
  {
    scores2_aligned[perm_shuf[i]] = out2.scores[i];
    qvals2_aligned[perm_shuf[i]] = out2.q_values[i];
  }

  // Pearson r on aligned scores.
  double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
  for (size_t i = 0; i < n; ++i)
  {
    const double x = out1.scores[i];
    const double y = scores2_aligned[i];
    sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y;
  }
  const double nf = static_cast<double>(n);
  const double num = nf * sxy - sx * sy;
  const double den = std::sqrt((nf * sxx - sx * sx) * (nf * syy - sy * sy));
  const double r = (den > 1e-15) ? num / den : 0.0;
  TEST_TRUE(r >= 0.99)

  // Jaccard on accepted-target set at q <= 0.01.
  std::set<size_t> a, b;
  for (size_t i = 0; i < n; ++i)
  {
    if (base.is_decoy[i]) continue;
    if (out1.q_values[i]   <= 0.01) a.insert(i);
    if (qvals2_aligned[i]  <= 0.01) b.insert(i);
  }
  std::set<size_t> inter, uni;
  std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                        std::inserter(inter, inter.begin()));
  std::set_union(a.begin(), a.end(), b.begin(), b.end(),
                 std::inserter(uni, uni.begin()));
  const double jaccard = uni.empty() ? 1.0
    : static_cast<double>(inter.size()) / static_cast<double>(uni.size());
  TEST_TRUE(jaccard >= 0.95)
}
END_SECTION
```

- [ ] **Step 3: Build & run**

```bash
cmake --build OpenMS-build --target Percolator_test -j$(nproc) && \
ctest --test-dir OpenMS-build -R '^Percolator_test$' --output-on-failure
```

Expected: passes. If Jaccard falls below 0.95 on moderately-separable data, input-order sensitivity is larger than documented — flag to user before relaxing the threshold.

- [ ] **Step 4: Commit**

```bash
git add src/tests/class_tests/openms/source/Percolator_test.cpp
git commit -m "$(cat <<'EOF'
test(Percolator): input-order invariance Jaccard characterization

New [EXTRA] section: shuffle input rows under fixed seed, run twice,
assert Pearson r >= 0.99 on aligned scores and Jaccard >= 0.95 on
accepted-target set at q <= 0.01. Documents the input-order drift the
docstring calls out so regressions that increase it are caught.
EOF
)"
```

---

## Task 4: P3.d — Seed-change negative control

**Files:**
- Modify: `src/tests/class_tests/openms/source/Percolator_test.cpp` — append new section before `END_TEST`.

**What this guards:** (a) catches a bug where the `seed` param isn't wired through (scores would then be bit-identical across seeds), (b) documents that different seeds converge to similar but not identical decision boundaries.

- [ ] **Step 1: Append new section**

Insert immediately before `END_TEST`:

```cpp
START_SECTION([EXTRA] seed change produces different but statistically equivalent scores)
{
  // Negative control: two different seeds must produce different scores
  // (catches seed-not-wired-through bugs) but Pearson r >= 0.9 and
  // |cA - cB| / max(cA, cB) <= 0.20 at q <= 0.01 (statistical equivalence).
  RescoreInput input = makeModeratelySeparableInput_(800, 789);

  Percolator p17, p42;
  {
    Param par = p17.getDefaults();
    par.setValue("seed", 17);
    p17.setParameters(par);
  }
  {
    Param par = p42.getDefaults();
    par.setValue("seed", 42);
    p42.setParameters(par);
  }

  RescoreOutput out17 = p17.rescore(input);
  RescoreOutput out42 = p42.rescore(input);

  // Not bit-identical.
  bool any_diff = false;
  for (size_t i = 0; i < out17.scores.size(); ++i)
  {
    if (out17.scores[i] != out42.scores[i]) { any_diff = true; break; }
  }
  TEST_EQUAL(any_diff, true)

  // Pearson r >= 0.9.
  const size_t n = out17.scores.size();
  double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
  for (size_t i = 0; i < n; ++i)
  {
    const double x = out17.scores[i];
    const double y = out42.scores[i];
    sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y;
  }
  const double nf = static_cast<double>(n);
  const double num = nf * sxy - sx * sy;
  const double den = std::sqrt((nf * sxx - sx * sx) * (nf * syy - sy * sy));
  const double r = (den > 1e-15) ? num / den : 0.0;
  TEST_TRUE(r >= 0.9)

  // Count ratio at q <= 0.01 (targets only).
  int c17 = 0, c42 = 0;
  for (size_t i = 0; i < n; ++i)
  {
    if (input.is_decoy[i]) continue;
    if (out17.q_values[i] <= 0.01) ++c17;
    if (out42.q_values[i] <= 0.01) ++c42;
  }
  const int maxc = std::max(c17, c42);
  const double ratio = (maxc > 0)
    ? static_cast<double>(std::abs(c17 - c42)) / static_cast<double>(maxc)
    : 0.0;
  TEST_TRUE(ratio <= 0.20)
}
END_SECTION
```

- [ ] **Step 2: Build & run**

```bash
cmake --build OpenMS-build --target Percolator_test -j$(nproc) && \
ctest --test-dir OpenMS-build -R '^Percolator_test$' --output-on-failure
```

Expected: passes. If `any_diff` is false, the `seed` param isn't wired through — serious bug, flag.

- [ ] **Step 3: Commit**

```bash
git add src/tests/class_tests/openms/source/Percolator_test.cpp
git commit -m "$(cat <<'EOF'
test(Percolator): seed-change negative control

New [EXTRA] section: runs with seed 17 vs 42, asserts scores are NOT
bit-identical (seed must be wired through) but Pearson r >= 0.9 and
|cA-cB|/max <= 0.20 at q <= 0.01 (seeds find statistically equivalent
answers).
EOF
)"
```

---

## Task 5: P1.a — Realistic-idXML library parity

**Files:**
- Modify: `src/tests/class_tests/openms/source/Percolator_subprocess_parity_test.cpp` — append new `START_SECTION` block before `END_TEST`.

**What this guards:** the existing library-layer parity tests all use synthetic data. P1.a brings a realistic feature distribution (from `PercolatorAdapter_1.idXML`) into the loop.

**Dataset caveat (spec):** `PercolatorAdapter_1.idXML` has only 40 PeptideHits. Metrics may be statistically thin; use lower `TEST_TRUE(matches >= N)` floors than synthetic sections.

- [ ] **Step 1: Append new section**

Insert this section before `END_TEST` in `Percolator_subprocess_parity_test.cpp`:

```cpp
///////////////////////////////////////////////////////////////////////////////
// Test 6: Realistic-idXML library parity (P1.a).
//
// Feeds PercolatorAdapter_1.idXML through both paths at the library layer.
// The .pin is the single source of truth for row order: we write it via
// PercolatorInfile::store, parse it back to build the RescoreInput, run
// subprocess on the same .pin, and align by SpecId. Ensures the in-process
// and subprocess paths see identical input, independent of any asymmetry
// between stampPinFeaturesOnHits and buildRescoreInput.
///////////////////////////////////////////////////////////////////////////////

START_SECTION([EXTRA] realistic idXML parity at library layer)
{
  const String bin = percolatorBinary();
  if (bin.empty())
  {
    TEST_EQUAL(true, true);  // skip silently
  }
  else
  {
    vector<ProteinIdentification> prs;
    PeptideIdentificationList pids;
    IdXMLFile().load(inputIdxml(), prs, pids);
    TEST_FALSE(pids.empty())

    int min_charge, max_charge;
    deriveChargeRange(pids, min_charge, max_charge);
    StringList feature_set = buildFeatureSet(prs, min_charge, max_charge);
    const std::string enz = "no_enzyme";

    // Write .pin (also stamps meta values under the hood).
    const String pin_path = File::getTemporaryFile();
    PercolatorInfile::store(pin_path, pids, feature_set, enz, min_charge, max_charge);

    // Parse .pin back to build RescoreInput that matches subprocess input
    // exactly. Avoids buildRescoreInput / stampPinFeaturesOnHits alignment
    // pitfalls because .pin is the ground truth for row order.
    TextFile pin;
    pin.load(pin_path);
    const size_t n_lines = std::distance(pin.begin(), pin.end());
    TEST_TRUE(n_lines >= 2)

    StringList header = ListUtils::create<String>(*pin.begin(), '\t');
    std::map<String, size_t> col_index;
    for (size_t c = 0; c < header.size(); ++c) col_index[header[c]] = c;

    StringList numeric = numericFeatures(feature_set);
    std::vector<String> pin_specids;
    RescoreInput ri;
    ri.feature_names = numeric;
    for (auto it = pin.begin() + 1; it != pin.end(); ++it)
    {
      StringList cols = ListUtils::create<String>(*it, '\t');
      // Tolerate Proteins column with embedded tabs by checking up to the
      // last numeric-feature column.
      bool ok = true;
      for (const auto& f : numeric)
      {
        if (col_index[f] >= cols.size()) { ok = false; break; }
      }
      if (!ok) continue;
      pin_specids.push_back(cols[col_index["SpecId"]]);
      ri.is_decoy.push_back(cols[col_index["Label"]] == "-1");
      ri.scan_numbers.push_back(cols[col_index["ScanNr"]].toInt());
      ri.spec_file_numbers.push_back(0);
      ri.exp_masses.push_back(cols[col_index["ExpMass"]].toDouble());
      ri.calc_masses.push_back(cols[col_index["CalcMass"]].toDouble());
      std::vector<double> feats;
      feats.reserve(numeric.size());
      for (const auto& f : numeric) feats.push_back(cols[col_index[f]].toDouble());
      ri.features.push_back(std::move(feats));
    }
    TEST_TRUE(ri.features.size() >= 20)

    // Run subprocess.
    SubprocessOut sub = runSubprocess(bin, pin_path, "-S 1");
    TEST_EQUAL(sub.exit_code, 0)

    // Run in-process, matching subprocess's isotonic PEP default.
    Percolator p;
    Param par = p.getDefaults();
    par.setValue("seed", 1);
    par.setValue("pep_method", "isotonic");
    p.setParameters(par);
    RescoreOutput out = p.rescore(ri);
    TEST_EQUAL(out.scores.size(), ri.features.size())

    // Pearson r, max |Δpep|, accepted-target sets at q in {0.01, 0.05, 0.10}.
    double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
    double max_dpep = 0;
    size_t matches = 0;
    for (size_t i = 0; i < out.scores.size(); ++i)
    {
      auto it = sub.triplets.find(pin_specids[i]);
      if (it == sub.triplets.end()) continue;
      matches++;
      const double x = out.scores[i];
      const double y = it->second.score;
      sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y;
      max_dpep = std::max(max_dpep, std::abs(out.peps[i] - it->second.pep));
    }
    const double nf = static_cast<double>(matches);
    const double num = nf * sxy - sx * sy;
    const double den = std::sqrt((nf * sxx - sx * sx) * (nf * syy - sy * sy));
    const double r = (den > 1e-15) ? num / den : 0.0;

    TEST_TRUE(matches >= 20)
    TEST_TRUE(r >= 0.99)  // TODO(percolator-3.08): tighten to 0.999 after upgrade
    TEST_TRUE(max_dpep <= 0.05)

    for (double thr : {0.01, 0.05, 0.10})
    {
      std::set<String> a_in, a_sub;
      for (size_t i = 0; i < out.scores.size(); ++i)
      {
        if (ri.is_decoy[i]) continue;
        auto it = sub.triplets.find(pin_specids[i]);
        if (it == sub.triplets.end()) continue;
        if (out.q_values[i] <= thr) a_in.insert(pin_specids[i]);
        if (it->second.qval <= thr) a_sub.insert(pin_specids[i]);
      }
      if (a_in != a_sub)
      {
        TEST_EQUAL(String("q=") + String(thr) + " accepted-set mismatch: in="
                   + String(a_in.size()) + " sub=" + String(a_sub.size()),
                   String("sets match"))
      }
      TEST_EQUAL(a_in == a_sub, true)
    }
  }
}
END_SECTION
```

- [ ] **Step 2: Build**

```bash
cmake --build OpenMS-build --target Percolator_subprocess_parity_test -j$(nproc)
```

Expected: clean build.

- [ ] **Step 3: Run**

```bash
ctest --test-dir OpenMS-build -R '^Percolator_subprocess_parity_test$' --output-on-failure
```

Expected: passes. If `r < 0.99` on the real data, that's a bigger 3.06↔3.08 algorithmic drift than expected — flag, then consider loosening to 0.95 with a `TODO(percolator-3.08)` note.

- [ ] **Step 4: Commit**

```bash
git add src/tests/class_tests/openms/source/Percolator_subprocess_parity_test.cpp
git commit -m "$(cat <<'EOF'
test(Percolator): library-layer parity on realistic idXML (P1.a)

New [EXTRA] section: writes PercolatorAdapter_1.idXML to .pin, parses
.pin back as the source of truth for row order, runs in-process and
subprocess, aligns by SpecId, asserts Pearson r >= 0.99 on scores,
accepted-target set equality at q <= 0.01/0.05/0.10, max |Δpep| <= 0.05.
EOF
)"
```

---

## Task 6: P1.b — Larger synthetic + reservoir-sampling parity

**Files:**
- Modify: `src/tests/class_tests/openms/source/Percolator_subprocess_parity_test.cpp` — extend `generateSyntheticData` signature with a size multiplier, append new `START_SECTION` block before `END_TEST`.

**What this guards:** the existing param-matrix §4 exercises `subset_max_train=200` on the 2000-row synthetic dataset — thin by comparison to production scale where the reservoir sampler actually matters. P1.b uses 20 000 rows with `subset_max_train=5000`.

- [ ] **Step 1: Extend `generateSyntheticData` with a size multiplier**

Locate this function in `Percolator_subprocess_parity_test.cpp` (currently ~line 209). Change its signature and body:

```cpp
void generateSyntheticData(RescoreInput& ri, std::mt19937& rng,
                           double size_mult = 1.0)
{
  ri.features.clear();
  ri.is_decoy.clear();
  ri.feature_names = StringList{"f0", "f1", "noise"};

  std::normal_distribution<double> noise(0.0, 0.5);
  std::normal_distribution<double> unit(0.0, 1.0);

  const size_t n_easy = static_cast<size_t>(500 * size_mult);
  const size_t n_hard = static_cast<size_t>(500 * size_mult);
  const size_t n_dec  = static_cast<size_t>(1000 * size_mult);
  for (size_t i = 0; i < n_easy; ++i)
  {
    ri.features.push_back({+2.0 + noise(rng), +1.0 + noise(rng), unit(rng)});
    ri.is_decoy.push_back(false);
  }
  for (size_t i = 0; i < n_hard; ++i)
  {
    ri.features.push_back({+0.5 + noise(rng), +0.5 * noise(rng), unit(rng)});
    ri.is_decoy.push_back(false);
  }
  for (size_t i = 0; i < n_dec; ++i)
  {
    ri.features.push_back({+0.0 + noise(rng), +0.0 * noise(rng), unit(rng)});
    ri.is_decoy.push_back(true);
  }

  // (Same PIN-OrderScanHash bookkeeping as before; unchanged.)
  const size_t n = ri.features.size();
  ri.scan_numbers.assign(n, 0);
  ri.spec_file_numbers.assign(n, 0);
  ri.exp_masses.assign(n, 1000.0);
  ri.calc_masses.assign(n, 1000.0);
  for (size_t i = 0; i < n; ++i)
  {
    ri.scan_numbers[i] = static_cast<int>(i + 1);
  }
}
```

The change from the existing version is: `double size_mult = 1.0` added to the signature, and `n_easy / n_hard / n_dec` now multiplied. Default `size_mult = 1.0` preserves the 2000-row behavior for all existing callers.

- [ ] **Step 2: Append new section**

Insert before `END_TEST`:

```cpp
///////////////////////////////////////////////////////////////////////////////
// Test 7: Reservoir-sampling parity at larger scale (P1.b).
//
// 20 000 rows with subset_max_train=5000 forces the reservoir sampler on
// both paths. Existing §4 case uses 200 on 2000 rows, which is informative
// but below the natural scale where sampling matters.
///////////////////////////////////////////////////////////////////////////////

START_SECTION([EXTRA] reservoir-sampling parity at 20k rows)
{
  const String bin = percolatorBinary();
  if (bin.empty())
  {
    TEST_EQUAL(true, true);  // skip
  }
  else
  {
    std::mt19937 rng(2026);
    RescoreInput ri;
    generateSyntheticData(ri, rng, /*size_mult=*/10.0);
    TEST_EQUAL(ri.features.size(), 20000)

    const String pin_path = File::getTemporaryFile();
    writePinFile(pin_path, ri);

    SubprocessOut sub = runSubprocess(bin, pin_path, "-S 1 -N 5000");
    TEST_EQUAL(sub.exit_code, 0)

    Percolator p;
    Param par = p.getDefaults();
    par.setValue("seed", 1);
    par.setValue("pep_method", "isotonic");
    par.setValue("subset_max_train", 5000);
    p.setParameters(par);
    RescoreOutput out = p.rescore(ri);

    double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
    size_t matches = 0;
    int tgt_q01_in = 0, tgt_q01_sub = 0;
    int tgt_q05_in = 0, tgt_q05_sub = 0;
    for (size_t i = 0; i < out.scores.size(); ++i)
    {
      char k[32]; std::snprintf(k, sizeof(k), "row_%08zu", i);
      auto it = sub.triplets.find(k);
      if (it == sub.triplets.end()) continue;
      matches++;
      const double x = out.scores[i];
      const double y = it->second.score;
      sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y;
      if (!ri.is_decoy[i])
      {
        if (out.q_values[i]   <= 0.01) tgt_q01_in++;
        if (it->second.qval <= 0.01) tgt_q01_sub++;
        if (out.q_values[i]   <= 0.05) tgt_q05_in++;
        if (it->second.qval <= 0.05) tgt_q05_sub++;
      }
    }
    const double nf = static_cast<double>(matches);
    const double num = nf * sxy - sx * sy;
    const double den = std::sqrt((nf * sxx - sx * sx) * (nf * syy - sy * sy));
    const double r = (den > 1e-15) ? num / den : 0.0;

    TEST_TRUE(matches >= 18000)
    TEST_TRUE(r >= 0.999)
    TEST_EQUAL(tgt_q01_in, tgt_q01_sub)
    TEST_EQUAL(tgt_q05_in, tgt_q05_sub)
  }
}
END_SECTION
```

- [ ] **Step 3: Build & run**

```bash
cmake --build OpenMS-build --target Percolator_subprocess_parity_test -j$(nproc) && \
ctest --test-dir OpenMS-build -R '^Percolator_subprocess_parity_test$' --output-on-failure
```

Expected: passes. Note: this section runs subprocess on 20 000 rows — may take noticeably longer than existing sections. If the default CTest timeout is too tight, flag.

- [ ] **Step 4: Commit**

```bash
git add src/tests/class_tests/openms/source/Percolator_subprocess_parity_test.cpp
git commit -m "$(cat <<'EOF'
test(Percolator): reservoir-sampling parity at 20k rows (P1.b)

Extends generateSyntheticData with size_mult; new [EXTRA] section runs
at 20k rows with subset_max_train=5000 to exercise the reservoir sampler
at a realistic scale. Tolerances unchanged from existing §2 (Pearson
r >= 0.999, exact target-count at q <= 0.01/0.05).
EOF
)"
```

---

## Task 7: P1.c — Multi-file PIN parity

**Files:**
- Modify: `src/tests/class_tests/openms/source/Percolator_subprocess_parity_test.cpp` — extend `writePinFile` with an optional `emit_filename` flag; append new section before `END_TEST`.

**What this guards:** the `OrderScanHash` CV-fold hash incorporates `specFileNr`. The existing tests keep `specFileNr=0` everywhere and write no `filename` column. P1.c fabricates a two-file case and confirms the CV fold hash agrees across paths.

**Subprocess behavior note:** `SetHandler.cpp:110-111` (upstream comment at `generateSyntheticData`) picks up `specFileNr` from an optional `filename` or `spectrafile` PIN column. Position of the column in the header is not fixed by upstream — parsing is column-name-driven.

- [ ] **Step 1: Extend `writePinFile`**

Locate this function in `Percolator_subprocess_parity_test.cpp` (currently ~line 263). Replace with:

```cpp
/// Write a minimal .pin file that percolator's command-line tool will accept.
/// Columns: SpecId, Label, ScanNr, ExpMass, CalcMass, <features>[, filename],
///          Peptide, Proteins.
/// Row ids `row_%08zu` so subprocess PSMId strings match in-process row indices.
/// When @p emit_filename is true, a `filename` column is appended after the
/// features with value `file_<spec_file_numbers[i]>`; subprocess then populates
/// PSMDescription::specFileNr via SetHandler.cpp:110-111.
void writePinFile(const String& path, const RescoreInput& ri,
                  bool emit_filename = false)
{
  std::ofstream f(path.c_str());
  if (!f) throw std::runtime_error("cannot open pin file for writing");
  f.setf(std::ios::scientific);
  f.precision(17);

  f << "SpecId\tLabel\tScanNr\tExpMass\tCalcMass";
  for (const auto& fn : ri.feature_names) f << "\t" << fn;
  if (emit_filename) f << "\tfilename";
  f << "\tPeptide\tProteins\n";

  const size_t n = ri.features.size();
  for (size_t i = 0; i < n; ++i)
  {
    char row_id[32];
    std::snprintf(row_id, sizeof(row_id), "row_%08zu", i);
    const int label = ri.is_decoy[i] ? -1 : +1;
    const int scan = (i < ri.scan_numbers.size())
      ? ri.scan_numbers[i]
      : static_cast<int>(i + 1);
    f << row_id << "\t" << label << "\t" << scan
      << "\t" << 1000.0 << "\t" << 1000.0;
    for (double v : ri.features[i]) f << "\t" << v;
    if (emit_filename)
    {
      const int fnr = (i < ri.spec_file_numbers.size())
        ? ri.spec_file_numbers[i] : 0;
      f << "\tfile_" << fnr;
    }
    f << "\t-.PEPTIDE.-\t"
      << (ri.is_decoy[i] ? "decoy_protein" : "target_protein") << "\n";
  }
}
```

The only difference from the existing function: optional `emit_filename` arg, and a `\tfile_<fnr>` column emitted after features when true. Default `false` preserves existing §2-§5 behavior.

- [ ] **Step 2: Append new section**

Insert before `END_TEST`:

```cpp
///////////////////////////////////////////////////////////////////////////////
// Test 8: Multi-file PIN parity (P1.c).
//
// Two sub-batches of synthetic rows with distinct specFileNr. The CV fold
// hash (OrderScanHash) incorporates specFileNr; this section confirms it
// agrees across paths when more than one file is present — a code path the
// other library-layer tests don't exercise (they all use specFileNr=0).
///////////////////////////////////////////////////////////////////////////////

START_SECTION([EXTRA] multi-file PIN parity)
{
  const String bin = percolatorBinary();
  if (bin.empty())
  {
    TEST_EQUAL(true, true);  // skip
  }
  else
  {
    std::mt19937 rng(2026);
    RescoreInput ri;
    generateSyntheticData(ri, rng);  // 2000 rows, specFileNr=0 by default

    // Split into two "files": first half -> specFileNr=0, second -> specFileNr=1.
    const size_t half = ri.features.size() / 2;
    for (size_t i = 0; i < ri.features.size(); ++i)
    {
      ri.spec_file_numbers[i] = (i < half) ? 0 : 1;
    }

    const String pin_path = File::getTemporaryFile();
    writePinFile(pin_path, ri, /*emit_filename=*/true);

    SubprocessOut sub = runSubprocess(bin, pin_path, "-S 1");
    TEST_EQUAL(sub.exit_code, 0)

    Percolator p;
    Param par = p.getDefaults();
    par.setValue("seed", 1);
    par.setValue("pep_method", "isotonic");
    p.setParameters(par);
    RescoreOutput out = p.rescore(ri);

    double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
    size_t matches = 0;
    int tgt_q01_in = 0, tgt_q01_sub = 0;
    int tgt_q05_in = 0, tgt_q05_sub = 0;
    for (size_t i = 0; i < out.scores.size(); ++i)
    {
      char k[32]; std::snprintf(k, sizeof(k), "row_%08zu", i);
      auto it = sub.triplets.find(k);
      if (it == sub.triplets.end()) continue;
      matches++;
      sx += out.scores[i];
      sy += it->second.score;
      sxx += out.scores[i] * out.scores[i];
      syy += it->second.score * it->second.score;
      sxy += out.scores[i] * it->second.score;
      if (!ri.is_decoy[i])
      {
        if (out.q_values[i] <= 0.01) ++tgt_q01_in;
        if (it->second.qval <= 0.01) ++tgt_q01_sub;
        if (out.q_values[i] <= 0.05) ++tgt_q05_in;
        if (it->second.qval <= 0.05) ++tgt_q05_sub;
      }
    }
    const double nf = static_cast<double>(matches);
    const double num = nf * sxy - sx * sy;
    const double den = std::sqrt((nf * sxx - sx * sx) * (nf * syy - sy * sy));
    const double r = (den > 1e-15) ? num / den : 0.0;

    TEST_TRUE(matches > 1000)
    TEST_TRUE(r >= 0.999)
    TEST_EQUAL(tgt_q01_in, tgt_q01_sub)
    TEST_EQUAL(tgt_q05_in, tgt_q05_sub)
  }
}
END_SECTION
```

- [ ] **Step 3: Build & run**

```bash
cmake --build OpenMS-build --target Percolator_subprocess_parity_test -j$(nproc) && \
ctest --test-dir OpenMS-build -R '^Percolator_subprocess_parity_test$' --output-on-failure
```

Expected: passes. If subprocess doesn't pick up the `filename` column (SetHandler column parsing differs from assumption), the CV-fold split will diverge and `r` will drop. If that happens: verify the column name (`filename` vs `spectrafile`) by inspecting SetHandler.cpp; or drop the multi-file section as not-portably-testable and flag.

- [ ] **Step 4: Commit**

```bash
git add src/tests/class_tests/openms/source/Percolator_subprocess_parity_test.cpp
git commit -m "$(cat <<'EOF'
test(Percolator): multi-file PIN parity (P1.c)

Extends writePinFile with an optional emit_filename flag; new [EXTRA]
section runs a two-file synthetic PIN through both paths, confirming the
OrderScanHash CV fold assignment agrees when specFileNr != 0 on some rows.
EOF
)"
```

---

## Task 8: P2 — Create adapter parity test file and wire it into CMake

**Files:**
- Create: `src/tests/class_tests/openms/source/PercolatorAdapter_parity_test.cpp`
- Modify: `src/tests/class_tests/openms/executables.cmake` — register the new test.
- Modify: `src/tests/class_tests/openms/CMakeLists.txt` — pass both `PERCOLATOR_BINARY` and `PERCOLATOR_ADAPTER_BINARY` env vars to the new test.

**What this guards:** no side-by-side regression gate exists today for `PercolatorAdapter -use_subprocess true/false`. Task 8 stands up the scaffolding; Task 9 fills in the parity logic.

- [ ] **Step 1: Create the test file with a skip-only body**

Write `src/tests/class_tests/openms/source/PercolatorAdapter_parity_test.cpp`:

```cpp
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
//
// Side-by-side parity test for the PercolatorAdapter TOPP tool: runs it twice
// on the same idXML with -use_subprocess true/false, aligns the output hits
// by (spectrum_reference, sequence, charge), and compares the stamped
// Percolator outputs (score, q-value, PEP).
//
// Gated on two env vars set by CMake:
//   PERCOLATOR_BINARY         -> forwarded as -percolator_executable
//   PERCOLATOR_ADAPTER_BINARY -> absolute path to bin/PercolatorAdapter
// Skips silently when either is unset.

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/File.h>

#include <cstdlib>
#include <string>

using namespace OpenMS;
using namespace std;

START_TEST(PercolatorAdapter_parity, "$Id$")

START_SECTION([EXTRA] adapter parity: -use_subprocess true vs false on same idXML)
{
  const char* perc = std::getenv("PERCOLATOR_BINARY");
  const char* adap = std::getenv("PERCOLATOR_ADAPTER_BINARY");
  if (!perc || !*perc || !adap || !*adap)
  {
    TEST_EQUAL(true, true);  // skip silently
  }
  else
  {
    // Parity logic lands in Task 9.
    TEST_EQUAL(true, true);
  }
}
END_SECTION

END_TEST
```

- [ ] **Step 2: Register the test in `executables.cmake`**

Open `src/tests/class_tests/openms/executables.cmake`, find the block listing the `Percolator*` tests (around line 248-251):

```cmake
  Percolator_test
  Percolator_subprocess_parity_test
  PercolatorInfile_test
  PercolatorOutfile_test
```

Add `PercolatorAdapter_parity_test` in alphabetical order:

```cmake
  Percolator_test
  Percolator_subprocess_parity_test
  PercolatorAdapter_parity_test
  PercolatorInfile_test
  PercolatorOutfile_test
```

- [ ] **Step 3: Pass env vars to the new test in `CMakeLists.txt`**

Open `src/tests/class_tests/openms/CMakeLists.txt`, find the `if (PERCOLATOR_BINARY_FOR_TEST)` block at ~line 75-81:

```cmake
if (PERCOLATOR_BINARY_FOR_TEST)
  set_tests_properties(Percolator_subprocess_parity_test PROPERTIES
    ENVIRONMENT "PERCOLATOR_BINARY=${PERCOLATOR_BINARY_FOR_TEST}")
  message(STATUS "Percolator_subprocess_parity_test: enabled (found ${PERCOLATOR_BINARY_FOR_TEST})")
else()
  message(STATUS "Percolator_subprocess_parity_test: subprocess sections skipped (percolator binary not on PATH nor in THIRDPARTY)")
endif()
```

Extend it to also set env on the new adapter-parity test:

```cmake
if (PERCOLATOR_BINARY_FOR_TEST)
  set_tests_properties(Percolator_subprocess_parity_test PROPERTIES
    ENVIRONMENT "PERCOLATOR_BINARY=${PERCOLATOR_BINARY_FOR_TEST}")
  set_tests_properties(PercolatorAdapter_parity_test PROPERTIES
    ENVIRONMENT "PERCOLATOR_BINARY=${PERCOLATOR_BINARY_FOR_TEST};PERCOLATOR_ADAPTER_BINARY=$<TARGET_FILE:PercolatorAdapter>")
  message(STATUS "Percolator_subprocess_parity_test: enabled (found ${PERCOLATOR_BINARY_FOR_TEST})")
  message(STATUS "PercolatorAdapter_parity_test: enabled (adapter=${PERCOLATOR_ADAPTER_BINARY_PATH_STATUS_MSG})")
else()
  message(STATUS "Percolator_subprocess_parity_test: subprocess sections skipped (percolator binary not on PATH nor in THIRDPARTY)")
  message(STATUS "PercolatorAdapter_parity_test: subprocess sections skipped (percolator binary not on PATH nor in THIRDPARTY)")
endif()
```

(Note: `$<TARGET_FILE:PercolatorAdapter>` resolves at generation time to the absolute path of the built PercolatorAdapter binary — so no extra `find_program` call is needed.)

- [ ] **Step 4: Build (re-run CMake first so the new test is picked up)**

```bash
cmake -B OpenMS-build && \
cmake --build OpenMS-build --target PercolatorAdapter_parity_test -j$(nproc)
```

Expected: clean CMake configure, clean build.

- [ ] **Step 5: Run the (trivial) test**

```bash
ctest --test-dir OpenMS-build -R '^PercolatorAdapter_parity_test$' --output-on-failure
```

Expected: passes. Confirms the env-var wiring is correct even before parity logic exists.

- [ ] **Step 6: Commit**

```bash
git add src/tests/class_tests/openms/source/PercolatorAdapter_parity_test.cpp \
        src/tests/class_tests/openms/executables.cmake \
        src/tests/class_tests/openms/CMakeLists.txt
git commit -m "$(cat <<'EOF'
test(PercolatorAdapter): add parity test scaffolding (P2)

Creates PercolatorAdapter_parity_test.cpp as a skip-only stub, registers
it in executables.cmake, and wires PERCOLATOR_BINARY and
PERCOLATOR_ADAPTER_BINARY env vars through CMake. Parity logic follows.
EOF
)"
```

---

## Task 9: P2 — Implement adapter parity logic

**Files:**
- Modify: `src/tests/class_tests/openms/source/PercolatorAdapter_parity_test.cpp` — replace the placeholder body with the full parity implementation.

**What this does:** invokes `PercolatorAdapter` twice on the same idXML with `-use_subprocess true` and `-use_subprocess false`, loads both output idXMLs, aligns hits by (spectrum_reference, sequence, charge), asserts both paths stamp the same Percolator meta keys, compares SVM score / q-value / PEP.

- [ ] **Step 1: Add includes and helpers, then rewrite the section**

Replace the entire contents of `PercolatorAdapter_parity_test.cpp` with:

```cpp
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
//
// Side-by-side parity test for the PercolatorAdapter TOPP tool: runs it twice
// on the same idXML with -use_subprocess true/false, aligns output hits by
// (spectrum_reference, sequence, charge), and compares the stamped Percolator
// outputs (score, q-value, PEP).
//
// Gated on two env vars set by CMake:
//   PERCOLATOR_BINARY         -> forwarded as -percolator_executable
//   PERCOLATOR_ADAPTER_BINARY -> absolute path to bin/PercolatorAdapter
// Skips silently when either is unset.

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace OpenMS;
using namespace std;

namespace
{

/// Composite key for 1:1 cross-run matching of adapter output hits.
String rowKey(const PeptideIdentification& pid, const PeptideHit& hit)
{
  String sr = String(pid.getSpectrumReference());
  if (sr.empty() && pid.metaValueExists("spectrum_id"))
  {
    sr = pid.getMetaValue("spectrum_id").toString();
  }
  return sr + "|" + hit.getSequence().toString() + "|" + String(hit.getCharge());
}

struct OutputTriplet
{
  double score = 0.0;
  double qval  = 0.0;
  double pep   = 0.0;
  bool   is_decoy = false;
  String key_for_score = "";  ///< meta key actually found on the hit
};

/// Extract the Percolator triplet from a hit, probing both meta-key schemes.
/// Returns {triplet, key_for_score}. If neither scheme is present, score
/// defaults to 0.0 and key_for_score is empty — caller treats as missing.
OutputTriplet extractTriplet(const PeptideIdentification& pid, const PeptideHit& hit)
{
  OutputTriplet t;
  t.is_decoy = hit.isDecoy();
  if (hit.metaValueExists("MS:1001492"))
  {
    t.score = static_cast<double>(hit.getMetaValue("MS:1001492"));
    t.qval  = hit.metaValueExists("MS:1001491")
      ? static_cast<double>(hit.getMetaValue("MS:1001491")) : 0.0;
    t.pep   = hit.metaValueExists("MS:1001493")
      ? static_cast<double>(hit.getMetaValue("MS:1001493")) : 0.0;
    t.key_for_score = "MS:1001492";
  }
  else if (hit.metaValueExists("percolator_score"))
  {
    t.score = static_cast<double>(hit.getMetaValue("percolator_score"));
    t.qval  = hit.metaValueExists("percolator_q_value")
      ? static_cast<double>(hit.getMetaValue("percolator_q_value")) : 0.0;
    t.pep   = hit.metaValueExists("percolator_pep")
      ? static_cast<double>(hit.getMetaValue("percolator_pep")) : 0.0;
    t.key_for_score = "percolator_score";
  }
  (void)pid;
  return t;
}

/// Invoke PercolatorAdapter on @p in_idxml, writing @p out_idxml. Returns 0
/// on success (exit code passed through).
int runAdapter(const String& adapter_bin,
               const String& percolator_bin,
               bool use_subprocess,
               const String& in_idxml,
               const String& out_idxml)
{
  const String stdout_log = File::getTemporaryFile();
  const String stderr_log = File::getTemporaryFile();

  std::ostringstream cmd;
  cmd << "\"" << adapter_bin << "\""
      << " -test"
      << " -in \""           << in_idxml     << "\""
      << " -out \""          << out_idxml    << "\""
      << " -out_type idXML"
      << " -percolator_executable \"" << percolator_bin << "\""
      << " -use_subprocess " << (use_subprocess ? "true" : "false")
      << " > \"" << stdout_log << "\""
      << " 2> \"" << stderr_log << "\"";
  return std::system(cmd.str().c_str());
}

/// Build a map[rowKey -> triplet] from an adapter output idXML.
std::map<String, OutputTriplet> loadTripletsByRowKey(const String& idxml)
{
  vector<ProteinIdentification> prs;
  PeptideIdentificationList pids;
  IdXMLFile().load(idxml, prs, pids);
  std::map<String, OutputTriplet> out;
  for (const auto& pid : pids)
  {
    for (const auto& hit : pid.getHits())
    {
      out[rowKey(pid, hit)] = extractTriplet(pid, hit);
    }
  }
  return out;
}

} // namespace

///////////////////////////////////////////////////////////////////////////////

START_TEST(PercolatorAdapter_parity, "$Id$")

START_SECTION([EXTRA] adapter parity: -use_subprocess true vs false on same idXML)
{
  const char* perc = std::getenv("PERCOLATOR_BINARY");
  const char* adap = std::getenv("PERCOLATOR_ADAPTER_BINARY");
  if (!perc || !*perc || !adap || !*adap)
  {
    TEST_EQUAL(true, true);  // skip silently
  }
  else
  {
    const String percolator_bin = String(perc);
    const String adapter_bin    = String(adap);
    const String in_idxml =
      OPENMS_GET_TEST_DATA_PATH("../../../topp/THIRDPARTY/PercolatorAdapter_1.idXML");

    const String out_sub = File::getTemporaryFile();
    const String out_inp = File::getTemporaryFile();

    TEST_EQUAL(runAdapter(adapter_bin, percolator_bin, /*use_subprocess=*/true,
                          in_idxml, out_sub), 0)
    TEST_EQUAL(runAdapter(adapter_bin, percolator_bin, /*use_subprocess=*/false,
                          in_idxml, out_inp), 0)

    auto sub = loadTripletsByRowKey(out_sub);
    auto inp = loadTripletsByRowKey(out_inp);

    TEST_TRUE(!sub.empty())
    TEST_TRUE(!inp.empty())

    // Prelude: both paths must stamp the same meta-key scheme. If one stamps
    // MS:1001492 and the other stamps percolator_score, comparison would
    // silently read zeros and pass spuriously.
    std::set<String> sub_keys, inp_keys;
    for (const auto& kv : sub) if (!kv.second.key_for_score.empty())
      sub_keys.insert(kv.second.key_for_score);
    for (const auto& kv : inp) if (!kv.second.key_for_score.empty())
      inp_keys.insert(kv.second.key_for_score);
    if (sub_keys != inp_keys)
    {
      std::ostringstream msg;
      msg << "meta-key mismatch: sub={ ";
      for (const auto& k : sub_keys) msg << k << " ";
      msg << "} vs inp={ ";
      for (const auto& k : inp_keys) msg << k << " ";
      msg << "}";
      TEST_EQUAL(String(msg.str()), String("meta keys match"))
    }
    TEST_EQUAL(sub_keys == inp_keys, true)

    // Same surviving hit set (symmetric difference = empty).
    std::set<String> only_sub, only_inp;
    for (const auto& kv : sub) if (!inp.count(kv.first)) only_sub.insert(kv.first);
    for (const auto& kv : inp) if (!sub.count(kv.first)) only_inp.insert(kv.first);
    TEST_EQUAL(only_sub.size(), 0)
    TEST_EQUAL(only_inp.size(), 0)

    // Pearson r, max |Δpep|, accepted-target sets at q in {0.01, 0.05, 0.10}.
    double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
    double max_dpep = 0;
    size_t matches = 0;
    struct Mismatch { String key; double s_sub, s_inp, q_sub, q_inp, p_sub, p_inp; };
    std::vector<Mismatch> mismatches;

    for (const auto& kv : inp)
    {
      auto it = sub.find(kv.first);
      if (it == sub.end()) continue;
      matches++;
      const double x = kv.second.score;
      const double y = it->second.score;
      sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y;
      const double dp = std::abs(kv.second.pep - it->second.pep);
      max_dpep = std::max(max_dpep, dp);
      if (std::abs(x - y) > 0.5 * std::max({std::abs(x), std::abs(y), 1.0})
          && mismatches.size() < 5)
      {
        mismatches.push_back({kv.first, y, x,
                              it->second.qval, kv.second.qval,
                              it->second.pep,  kv.second.pep});
      }
    }
    const double nf = static_cast<double>(matches);
    const double num = nf * sxy - sx * sy;
    const double den = std::sqrt((nf * sxx - sx * sx) * (nf * syy - sy * sy));
    const double r = (den > 1e-15) ? num / den : 0.0;

    TEST_TRUE(matches >= 20)
    if (r < 0.99)
    {
      std::ostringstream msg;
      msg << "r=" << r << " (expected >= 0.99); first mismatches: ";
      for (const auto& m : mismatches)
      {
        msg << "[" << m.key << " s_sub=" << m.s_sub << " s_inp=" << m.s_inp
            << " q_sub=" << m.q_sub << " q_inp=" << m.q_inp
            << " p_sub=" << m.p_sub << " p_inp=" << m.p_inp << "] ";
      }
      TEST_EQUAL(String(msg.str()), String("r ok"))
    }
    TEST_TRUE(r >= 0.99)  // TODO(percolator-3.08): tighten to 0.999 after upgrade
    TEST_TRUE(max_dpep <= 0.05)

    // Accepted-target set equality at q in {0.01, 0.05, 0.10}.
    for (double thr : {0.01, 0.05, 0.10})
    {
      std::set<String> acc_sub, acc_inp;
      for (const auto& kv : inp)
      {
        if (kv.second.is_decoy) continue;
        auto it = sub.find(kv.first);
        if (it == sub.end()) continue;
        if (kv.second.qval <= thr)  acc_inp.insert(kv.first);
        if (it->second.qval <= thr) acc_sub.insert(kv.first);
      }
      if (acc_sub != acc_inp)
      {
        TEST_EQUAL(String("q=") + String(thr) + " accepted-set mismatch: sub="
                   + String(acc_sub.size()) + " inp=" + String(acc_inp.size()),
                   String("sets match"))
      }
      TEST_EQUAL(acc_sub == acc_inp, true)
    }
  }
}
END_SECTION

END_TEST
```

- [ ] **Step 2: Build**

```bash
cmake --build OpenMS-build --target PercolatorAdapter_parity_test -j$(nproc)
```

Expected: clean build.

- [ ] **Step 3: Run**

```bash
ctest --test-dir OpenMS-build -R '^PercolatorAdapter_parity_test$' --output-on-failure
```

Expected: passes on a machine with `PERCOLATOR_BINARY` set. If it skips (no binary), that's also acceptable.

**Likely failure modes and triage:**
- Meta-key mismatch (`sub_keys != inp_keys`): the in-process path's adapter wiring stamps `percolator_*` while subprocess stamps `MS:*` (or vice versa). Real finding — not something to "loosen"; flag and stop. This is the spec check working as designed.
- Different surviving hit set: the adapter applies path-specific filtering (e.g., one path drops TDC-weeded rows, the other keeps them at `q=1`). Investigate which filter differs.
- `r < 0.99` with no meta-key issue: adapter-level drift larger than expected. First examine the `mismatches` output for patterns (all targets? all decoys? score sign flipped?). Then decide whether to tighten the adapter or loosen the tolerance with a `TODO(percolator-3.08)`.

- [ ] **Step 4: Commit**

```bash
git add src/tests/class_tests/openms/source/PercolatorAdapter_parity_test.cpp
git commit -m "$(cat <<'EOF'
test(PercolatorAdapter): implement -use_subprocess parity assertions (P2)

Runs PercolatorAdapter twice on PercolatorAdapter_1.idXML with
-use_subprocess true/false, aligns hits by (spectrum_reference, sequence,
charge), and asserts:

- Both paths stamp the same Percolator meta-key scheme (MS:1001491/2/3
  vs percolator_score/_q_value/_pep). Prelude check — prevents a silent
  pass where both sides read zeros for missing keys.
- Symmetric difference of surviving hit sets is empty.
- Pearson r on SVM scores >= 0.99 (TODO(percolator-3.08): 0.999).
- Accepted-target set equality at q <= 0.01/0.05/0.10.
- max |Δpep| <= 0.05.

Dumps the first 5 mismatches on failure for triage.
EOF
)"
```

---

## Task 10: Full-suite validation + README note

**Files:**
- None created. Just run everything; if the spec says something is expected to be exercised, confirm it is.

- [ ] **Step 1: Clean build from scratch, then run all three tests**

```bash
cmake --build OpenMS-build -j$(nproc) && \
ctest --test-dir OpenMS-build -R 'Percolator_test|Percolator_subprocess_parity_test|PercolatorAdapter_parity_test' --output-on-failure
```

Expected: all three pass. If `PERCOLATOR_BINARY` is discoverable, the gated sections run too.

- [ ] **Step 2: Confirm skip-path still works (if possible)**

On a system or with an env that has no `percolator` binary, confirm the subprocess-gated sections skip silently and still pass. If the current host has the binary discoverable, this is optional — the skip path is covered by the existing `Percolator_subprocess_parity_test` wiring unchanged from before.

```bash
PERCOLATOR_BINARY= PERCOLATOR_ADAPTER_BINARY= \
  ctest --test-dir OpenMS-build -R 'Percolator_subprocess_parity_test|PercolatorAdapter_parity_test' --output-on-failure
```

Expected: both tests pass (the sections inside skip silently).

- [ ] **Step 3: No commit needed**

This is validation only. If everything passes, the feature is complete.

---

## Self-review

**Spec coverage check:**

| Spec section | Task(s) | Covered |
|---|---|---|
| P1.a Realistic idXML parity | Task 5 | ✓ |
| P1.b Larger synthetic (reservoir) | Task 6 | ✓ |
| P1.c Multi-file PIN | Task 7 | ✓ |
| P2 Adapter parity | Tasks 8–9 | ✓ |
| P2 meta-key prelude | Task 9 | ✓ |
| P2 first-5-mismatches dump | Task 9 | ✓ |
| P3.a Same-instance consecutive | Task 1 | ✓ |
| P3.b Cross-instance bit-identical | Task 2 | ✓ |
| P3.c Input-order Jaccard | Task 3 | ✓ |
| P3.d Seed-change negative | Task 4 | ✓ |
| Tolerance table | All tasks | ✓ (per-task assertions match) |
| Helper: `writePinFile` emit_filename | Task 7 | ✓ |
| Helper: `generateSyntheticData` size_mult | Task 6 | ✓ |
| CMake wiring for new test + env vars | Task 8 | ✓ |
| `TODO(percolator-3.08)` markers | Tasks 5 & 9 | ✓ |
| Full-suite green run | Task 10 | ✓ |

No gaps.

**Placeholder scan:** no "TBD" / "TODO later" / "similar to Task N"-style references. Each task shows the actual code. The `TODO(percolator-3.08)` markers in the code are planned future-maintenance tags per the spec, not placeholders.

**Type consistency:**
- `RescoreInput`, `RescoreOutput`, `Percolator`, `PercolatorModel` — used consistently across tasks, matching the production API.
- `OutputTriplet` — new struct in Task 9 only; not used elsewhere.
- Helper signatures (`writePinFile`, `generateSyntheticData`, `runAdapter`, `rowKey`) — consistent between their definition task and their use.
- `generateSyntheticData(..., size_mult=1.0)` default keeps all existing callers binary-compatible — double-checked in Task 6.

No inconsistencies.
