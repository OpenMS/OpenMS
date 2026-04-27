# SNES X/B/Z Mother Truncation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stop dropping whole SNES mother peptides when their proposed window overlaps an ambiguous residue (X/B/Z); instead, truncate to unambiguous spans so shorter valid realizations at the same anchor survive (issue #9192 item 2).

**Architecture:** In `FragmentIndex::generateSNESMothers_`, keep the existing protein-level fast path for the no-ambiguous-residue case unchanged. For proteins containing X/B/Z, compute contiguous unambiguous spans up front and run the existing Single-N / Single-C sweeps within each span. The `emitMother` lambda loses its inner X/B/Z check (now redundant — slow path is gated structurally by span boundaries).

**Tech Stack:** C++ 20, OpenMS test framework (`START_TEST` / `START_SECTION` / `TEST_EQUAL`), CMake/ctest.

**Spec:** `docs/superpowers/specs/2026-04-27-snes-xbz-truncate-design.md`

---

## File Structure

| File | Role | Change |
|------|------|--------|
| `src/openms/source/ANALYSIS/ID/FragmentIndex.cpp` | SNES mother generation logic | Modify `generateSNESMothers_` body (lines 737-830). Simplify `emitMother` (remove X/B/Z check, lines 774-784). Add per-span sweep on slow path. |
| `src/tests/class_tests/openms/source/FragmentIndex_test.cpp` | Test suite | Add two `START_SECTION` blocks beside the existing X/B/Z rejection test (line 2419). |

No header changes. No build-system changes (the test binary already includes the new sections automatically).

---

## Task 1: Add failing tests for X/B/Z mother truncation

**Files:**
- Modify: `src/tests/class_tests/openms/source/FragmentIndex_test.cpp` — append two new sections immediately after the existing block ending at line 2462 (`END_SECTION` of the X/B/Z rejection test).

- [ ] **Step 1: Add the two new test sections**

After the existing `END_SECTION` for `(SNES mother generation rejects ambiguous residue spans (X/B/Z))` (line 2462), insert:

```cpp
START_SECTION((SNES mother generation truncates Single-N mother to unambiguous prefix on X/B/Z))
{
  // Issue #9192 item 2: a Single-N mother anchored at position 0 with proposed
  // length 12 spans the X at position 8. The whole mother used to be dropped;
  // now the unambiguous prefix [0, 8) length 8 must still be emitted.
  const std::vector<FASTAFile::FASTAEntry> entries{
      {"p", "p", "ACDEFGHIXKLMNPQSTVWY"}}; // X at 0-based position 8

  FragmentIndex_test fi;
  auto p = fi.getParameters();
  p.setValue("peptide:enzyme_specificity", "none");
  p.setValue("peptide:min_size", 8);
  p.setValue("peptide:max_size", 12);
  p.setValue("peptide:min_mass", 0);
  p.setValue("peptide:max_mass", 50000);
  p.setValue("modifications:variable", std::vector<std::string>{});
  p.setValue("modifications:fixed", std::vector<std::string>{});
  p.setValue("snes_enabled", "true");
  fi.setParameters(p);
  fi.build(entries);

  // Truncated Single-N mother [0, 8) of length 8 must exist.
  bool found_truncated_prefix = false;
  for (const auto& mother : fi.getPeptides())
  {
    const bool is_single_c = FragmentIndex::isSingleCMother(mother.mod_bitmask_);
    if (!is_single_c
        && mother.sequence_.first == 0
        && mother.sequence_.second == 8)
    {
      found_truncated_prefix = true;
      break;
    }
  }
  TEST_EQUAL(found_truncated_prefix, true)

  // Invariant: no kept mother spans the X at position 8.
  for (const auto& mother : fi.getPeptides())
  {
    const size_t start = mother.sequence_.first;
    const size_t end = start + mother.sequence_.second;
    TEST_EQUAL(start > 8u || end <= 8u, true)
  }
}
END_SECTION

START_SECTION((SNES mother generation truncates Single-C mother to unambiguous suffix on X/B/Z))
{
  // Issue #9192 item 2: a Single-C mother anchored at the last residue (j=19)
  // with proposed length 12 spans [8, 20) and covers the X. The whole mother
  // used to be dropped; now the unambiguous suffix [9, 20) length 11 must
  // still be emitted (length capped per-span at e - s = 11).
  const std::vector<FASTAFile::FASTAEntry> entries{
      {"p", "p", "ACDEFGHIXKLMNPQSTVWY"}}; // X at 0-based position 8

  FragmentIndex_test fi;
  auto p = fi.getParameters();
  p.setValue("peptide:enzyme_specificity", "none");
  p.setValue("peptide:min_size", 8);
  p.setValue("peptide:max_size", 12);
  p.setValue("peptide:min_mass", 0);
  p.setValue("peptide:max_mass", 50000);
  p.setValue("modifications:variable", std::vector<std::string>{});
  p.setValue("modifications:fixed", std::vector<std::string>{});
  p.setValue("snes_enabled", "true");
  fi.setParameters(p);
  fi.build(entries);

  // Truncated Single-C mother [9, 20) of length 11 must exist.
  bool found_truncated_suffix = false;
  for (const auto& mother : fi.getPeptides())
  {
    const bool is_single_c = FragmentIndex::isSingleCMother(mother.mod_bitmask_);
    if (is_single_c
        && mother.sequence_.first == 9
        && mother.sequence_.second == 11)
    {
      found_truncated_suffix = true;
      break;
    }
  }
  TEST_EQUAL(found_truncated_suffix, true)
}
END_SECTION
```

- [ ] **Step 2: Re-run cmake to pick up the modified test source**

Run: `cmake --build OpenMS-build --target FragmentIndex_test -j$(nproc)`
Expected: build succeeds; the binary is rebuilt with the two new sections compiled in.

- [ ] **Step 3: Run the new tests to confirm they fail**

Run from project root:
```bash
cd OpenMS-build && ctest -R "^FragmentIndex_test$" --output-on-failure
```
or equivalently run the binary directly:
```bash
OpenMS-build/src/tests/class_tests/openms/bin/FragmentIndex_test 2>&1 | grep -A2 "truncates Single"
```
Expected: both new sections fail. The Single-N section reports `found_truncated_prefix == false`; the Single-C section reports `found_truncated_suffix == false`. Existing sections (including the original X/B/Z rejection test at line 2419) still pass.

- [ ] **Step 4: Commit the failing tests**

```bash
git add src/tests/class_tests/openms/source/FragmentIndex_test.cpp
git commit -m "$(cat <<'EOF'
test(FragmentIndex): SNES mother truncation on X/B/Z (#9192)

Add two failing tests asserting that Single-N / Single-C mother peptides
whose proposed window spans an X/B/Z residue still emit a truncated mother
covering the unambiguous prefix or suffix at the same anchor. Currently the
whole mother is dropped — see issue #9192 item 2.
EOF
)"
```

---

## Task 2: Implement span-split sweep on the slow path

**Files:**
- Modify: `src/openms/source/ANALYSIS/ID/FragmentIndex.cpp` — replace the body of `generateSNESMothers_` between lines 750 and 829 (per-protein loop body up to and including the Single-C sweep). Keep the OMP wrapper, thread bookkeeping, mass-constant precompute, and post-loop merge / sort / log unchanged.

- [ ] **Step 1: Replace the per-protein loop body**

Locate this block (currently FragmentIndex.cpp:750-829):

```cpp
      // Fast path: for the overwhelmingly common case of a protein with no
      // ambiguous residues, a single linear scan tells us we can skip the per-
      // mother X/B/Z check entirely. For the rare protein with X/B/Z, we fall
      // back to checking only the next-bad-position via std::string::find_first_of
      // starting from the mother's start index — still O(length) in the worst
      // case, but only if the mother actually overlaps a bad region.
      //
      // Follow-up: a protein-wide skip here is overly conservative. A max-length
      // mother spanning an `X` is discarded, but shorter realizations with the
      // same anchor *not* covering the `X` could be valid. CodeRabbit #1.
      // Correct fix: split the protein into contiguous unambiguous spans and
      // generate mothers per-span. Deferred to v1.1.
      const bool protein_has_ambiguous = seq.find_first_of("XBZ") != std::string::npos;

      // Honor peptide:max_size=0 as "no maximum" (the documented semantics of
      // the non-SNES path). Using raw peptide_max_length_ in std::min would give
      // length 0 and an empty SNES index.
      const size_t effective_max_length = (peptide_max_length_ == 0) ? L : peptide_max_length_;

      auto emitMother = [&](size_t start, size_t length, bool is_single_c)
      {
        if (length < peptide_min_length_) return;
        const char* seq_ptr = seq.c_str() + start;

        // Reject mothers containing ambiguous codes — any realized sub-peptide
        // spanning an X/B/Z would fail AASequence::fromString downstream.
        if (protein_has_ambiguous)
        {
          const size_t bad = seq.find_first_of("XBZ", start);
          if (bad != std::string::npos && bad < start + length)
          {
            skipped_peptides.fetch_add(1);
            return;
          }
        }

        double mass = base_sum_constants;
        for (size_t k = 0; k < length; ++k)
        {
          const unsigned char aa = static_cast<unsigned char>(seq_ptr[k]);
          mass += residue_mass_table_[aa] + fixed_mod_deltas_[aa];
        }
        const float mz = static_cast<float>(mass);
        // Only the lower bound is safe at mother-generation time: shorter
        // realizations of a mother whose total mass exceeds peptide_max_mass_
        // can still fall within the user's configured mass range. Enforce the
        // upper bound at realization time via the precursor-tolerance window
        // (which is always <= peptide_max_mass_ for observed spectra). CodeRabbit #5.
        if (mz < peptide_min_mass_) return;

        const uint32_t kind_bits = is_single_c ? SNES_KIND_BIT_MASK : 0u;
        thread_peptides[tid].emplace_back(
            static_cast<UInt32>(protein_idx),
            kind_bits,
            std::make_pair(static_cast<uint16_t>(start), static_cast<uint16_t>(length)),
            mz);
      };

      // Single-N mothers: anchored at position i, span the longest possible peptide
      // starting there (capped at effective_max_length). i sweeps [0, L - min_length].
      for (size_t i = 0; i + peptide_min_length_ <= L; ++i)
      {
        const size_t length = std::min<size_t>(effective_max_length, L - i);
        emitMother(i, length, /*is_single_c=*/false);
      }

      // Single-C mothers: anchored at position j (last residue), span the longest
      // possible peptide ending there. j sweeps [min_length - 1, L - 1].
      // When j + 1 <= effective_max_length the mother happens to coincide with a
      // Single-N mother at position 0 — that's harmless redundancy: both emit a
      // different ion series into the index, so there's no duplicate fragment.
      // Guard against min_length=0: j would wrap to SIZE_MAX. Clamp to 1 locally;
      // this is the only code path sensitive to the min_length=0 edge case.
      const size_t snes_min_length = std::max<size_t>(1, peptide_min_length_);
      for (size_t j = snes_min_length - 1; j < L; ++j)
      {
        const size_t length = std::min<size_t>(effective_max_length, j + 1);
        const size_t start = j + 1 - length;
        emitMother(start, length, /*is_single_c=*/true);
      }
```

Replace it with:

```cpp
      // Honor peptide:max_size=0 as "no maximum" (the documented semantics of
      // the non-SNES path). Using raw peptide_max_length_ in std::min would give
      // length 0 and an empty SNES index.
      const size_t effective_max_length = (peptide_max_length_ == 0) ? L : peptide_max_length_;

      // Mass-compute + filter + emit. No X/B/Z check here: the fast path runs
      // only when the protein is unambiguous, and the slow path drives
      // correctness via span boundaries.
      auto emitMother = [&](size_t start, size_t length, bool is_single_c)
      {
        if (length < peptide_min_length_) return;
        const char* seq_ptr = seq.c_str() + start;

        double mass = base_sum_constants;
        for (size_t k = 0; k < length; ++k)
        {
          const unsigned char aa = static_cast<unsigned char>(seq_ptr[k]);
          mass += residue_mass_table_[aa] + fixed_mod_deltas_[aa];
        }
        const float mz = static_cast<float>(mass);
        // Only the lower bound is safe at mother-generation time: shorter
        // realizations of a mother whose total mass exceeds peptide_max_mass_
        // can still fall within the user's configured mass range. Enforce the
        // upper bound at realization time via the precursor-tolerance window
        // (which is always <= peptide_max_mass_ for observed spectra). CodeRabbit #5.
        if (mz < peptide_min_mass_) return;

        const uint32_t kind_bits = is_single_c ? SNES_KIND_BIT_MASK : 0u;
        thread_peptides[tid].emplace_back(
            static_cast<UInt32>(protein_idx),
            kind_bits,
            std::make_pair(static_cast<uint16_t>(start), static_cast<uint16_t>(length)),
            mz);
      };

      // Single-N mothers anchored at every position in [s, e - min_length], length
      // capped at effective_max_length and at the span end. Single-C mothers
      // anchored at every position j in [s + snes_min_length - 1, e - 1] with the
      // same length cap. snes_min_length guards the peptide_min_length_=0 corner
      // case (j would wrap to SIZE_MAX otherwise).
      const size_t snes_min_length = std::max<size_t>(1, peptide_min_length_);
      auto sweepSpan = [&](size_t s, size_t e)
      {
        if (e <= s || e - s < peptide_min_length_)
        {
          if (e > s) skipped_peptides.fetch_add(1);
          return;
        }
        for (size_t i = s; i + peptide_min_length_ <= e; ++i)
        {
          const size_t length = std::min<size_t>(effective_max_length, e - i);
          emitMother(i, length, /*is_single_c=*/false);
        }
        for (size_t j = s + snes_min_length - 1; j < e; ++j)
        {
          const size_t length = std::min<size_t>(effective_max_length, j + 1 - s);
          const size_t start = j + 1 - length;
          emitMother(start, length, /*is_single_c=*/true);
        }
      };

      // Fast path: no X/B/Z anywhere — sweep the whole protein as a single span.
      // Slow path: walk find_first_of("XBZ") to split into contiguous unambiguous
      // spans. Issue #9192 item 2: previously the whole mother was dropped on any
      // X/B/Z overlap; truncating to the unambiguous prefix/suffix at the same
      // anchor preserves valid shorter realizations.
      const size_t first_bad = seq.find_first_of("XBZ");
      if (first_bad == std::string::npos)
      {
        sweepSpan(0, L);
      }
      else
      {
        size_t p = 0;
        size_t bad = first_bad;
        while (true)
        {
          sweepSpan(p, bad);
          p = bad + 1;
          if (p >= L) break;
          bad = seq.find_first_of("XBZ", p);
          if (bad == std::string::npos) { sweepSpan(p, L); break; }
        }
      }
```

- [ ] **Step 2: Build**

Run: `cmake --build OpenMS-build --target FragmentIndex_test -j$(nproc)`
Expected: build succeeds with no warnings on the modified file. If you see an "unused variable `protein_has_ambiguous`" warning, you missed deleting the old declaration.

- [ ] **Step 3: Run the previously-failing tests — they must now pass**

```bash
cd OpenMS-build && ctest -R "^FragmentIndex_test$" --output-on-failure
```
Expected: ALL `FragmentIndex_test` sections pass, including the two new ones from Task 1 (`(SNES mother generation truncates Single-N mother to unambiguous prefix on X/B/Z)` and `(SNES mother generation truncates Single-C mother to unambiguous suffix on X/B/Z)`).

- [ ] **Step 4: Verify the existing X/B/Z rejection test still passes**

The original test at line 2419 (`(SNES mother generation rejects ambiguous residue spans (X/B/Z))`) uses `peptide_min=peptide_max=8`, so no mother is truncatable — both invariants (no kept mother spans the X; the start-0 length-8 prefix is kept) still hold under the new code. Confirm by reading the ctest output: that section's `TEST_EQUAL` lines must report passes, not failures.

- [ ] **Step 5: Commit**

```bash
git add src/openms/source/ANALYSIS/ID/FragmentIndex.cpp
git commit -m "$(cat <<'EOF'
fix(FragmentIndex): truncate SNES mothers to unambiguous spans on X/B/Z

Previously, when a SNES mother peptide's proposed window overlapped an
X/B/Z residue, the whole mother was dropped — losing valid shorter
realizations anchored at the same terminus that did not span the
ambiguous residue.

Now generateSNESMothers_ splits each protein containing X/B/Z into
contiguous unambiguous spans and runs the Single-N / Single-C sweeps
within each span, with the length cap enforced per-span. Proteins
without X/B/Z take the same single-span fast path as before. The
emitMother lambda no longer carries an X/B/Z check — span boundaries
are now structural.

Closes #9192 (item 2; item 1 tracked separately).
EOF
)"
```

---

## Task 3: Run the full FragmentIndex test suite to verify no regressions

**Files:** none modified.

- [ ] **Step 1: Run the complete `FragmentIndex_test` binary**

```bash
cd OpenMS-build && ctest -R "^FragmentIndex_test$" --output-on-failure -V
```
Expected: every `START_SECTION` reports a pass. No `FAILED` line. The `Generated N SNES mothers (M skipped...)` log message may show a different M than before for the X-containing test (counter semantic change) — that is informational, not a test assertion.

- [ ] **Step 2: Run the broader SNES-related test set (PeptideDataBaseSearchFI / ProSE end-to-end)**

```bash
cd OpenMS-build && ctest -R "FragmentIndex|PeptideDataBaseSearch|ProSE" --output-on-failure
```
Expected: all pass. The change is local to `generateSNESMothers_` and only affects mother enumeration on X/B/Z-containing proteins, so unrelated tests are unaffected.

- [ ] **Step 3: If anything fails, STOP and investigate**

Per CLAUDE.md: never bulk-replace expected values or regenerate test outputs. Investigate the root cause. If a SNES test using a FASTA with X/B/Z asserts on a specific mother count, the count may need updating because previously-rejected mothers are now (correctly) emitted. Confirm the new count by hand before changing any expected value.

- [ ] **Step 4: No commit needed if tests already pass.** This task is purely verification.

---

## Self-review notes

- **Spec coverage:** Span computation → Task 2 Step 1. Sweep semantics → Task 2 Step 1 (`sweepSpan` lambda). Counter semantics → Task 2 Step 1 (`if (e > s) skipped_peptides.fetch_add(1)` on too-short spans). Edge cases (whole protein ambiguous, span shorter than min, max=0, min=0, multiple X, adjacent X) → all covered by the early-out + `sweepSpan` body. Tests for both Single-N and Single-C truncation → Task 1. Existing test 2419 unchanged → Task 3 Step 1.
- **Placeholder scan:** none.
- **Type consistency:** `sweepSpan(size_t, size_t)`, `emitMother(size_t, size_t, bool)`, `Peptide::sequence_` uses `pair<uint16_t, uint16_t>` (matches existing). `FragmentIndex::isSingleCMother` is the public static accessor declared at `FragmentIndex.h:254` and already used throughout the test file (e.g. lines 1282, 1344, 1392, 1508, 1749, 1880).
