# Percolator vendored tree migration to eb157f7 — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Sync the vendored percolator tree from upstream `rel-3-08-01` (SHA `febeef3`) to commit `eb157f7` (PR #399, merged 2026-04-23), absorbing the new I-spline PEP calibration + small cleanups, without breaking the OpenMS-side `Percolator` wrapper or its tests.

**Architecture:** Three phases. **(1) Preparation** — capture pre-migration test output as baseline, update vendored-tree metadata (whitelist / sources.cmake / UPSTREAM_COMMIT). **(2) Sync and adapt** — copy the changed + new upstream files, apply namespace-wrap patches to the new content, regenerate the two patch hunk blocks that no longer apply cleanly, edit our static `Globals.h` stub. **(3) Rebuild and validate** — clean rebuild, run all three percolator test suites, update any tolerances that shifted, commit. No changes to `Percolator.cpp`'s wrapper API or its call sites are expected: `Scores::calcPep(bool, bool, bool)` signature is preserved upstream.

**Tech Stack:** Bash scripts (`sync-from-upstream.sh`), `gh` CLI for upstream fetch, C++ vendored tree, CMake-driven OpenMS build, OpenMS class-test framework, `ctest` for end-to-end runs.

## Analysis reference

The cascading-impact analysis at `docs/superpowers/plans/` (session transcript 2026-04-24) maps upstream `eb157f7` changes to our vendored tree:

- **Preserved:** `Scores::calcPep` signature, `Scores::calcQvals` (2-arg overload still exists), `CrossValidation` constructor args, `Normalizer` / `SanityCheck` / `Globals` / `FeatureNames` APIs our wrapper uses.
- **New files required:** `MonotoneRegressor.h`, `MonotoneRegressor.cpp`, `PAVARegressor.h`, `ISplineTRRRegressor.h`.
- **Full rewrites:** `IsotonicPEP.h`, `IsotonicPEP.cpp` — all namespace-wrap hunks fail, must be regenerated.
- **Context shifts:** `Scores.cpp` three hunks (`"\% FDR"` → `"% FDR"`, `getScoreLabelPairs`, `weedOutRedundant`), `Scores.h` one hunk (new `calcQvals` overload).
- **Dead code removal:** `Globals.cpp` loses `getXMLDir()` body — our static `Globals.h` stub must drop the declaration.
- **Algorithmic shift:** `pep_method="isotonic"` now adds Bayesian 0.5/n pseudocount + changes TDC anchor strategy; `pep_method="logistic_regression"` algorithm completely replaced (TRR I-spline). `pep_method="nonparametric"` (spline/`PosteriorEstimator::estimatePEP`) is **unchanged**.
- **No new external dependencies:** Eigen already linked in OpenMS; `Eigen/Sparse` header is dropped but `Eigen/Dense` remains in the new `ISplineTRRRegressor.h`.

## File structure

| File | Action | Purpose |
|---|---|---|
| `src/openms/thirdparty/percolator/whitelist.txt` | Modify | Add 4 new files to the sync whitelist |
| `src/openms/thirdparty/percolator/sources.cmake` | Modify | Add `MonotoneRegressor.cpp` to the compiled source list |
| `src/openms/thirdparty/percolator/UPSTREAM_COMMIT` | Modify | Update pinned commit SHA |
| `src/openms/thirdparty/percolator/UPSTREAM.md` | Modify | Update commit SHA + ref name in docs |
| `src/openms/thirdparty/percolator/sync-from-upstream.sh` | Modify | Update default `UPSTREAM_REF` to `eb157f7` |
| `src/openms/thirdparty/percolator/MonotoneRegressor.h` | Create (sync + wrap) | Abstract base + factory |
| `src/openms/thirdparty/percolator/MonotoneRegressor.cpp` | Create (sync + wrap) | Factory impl; pulls in `PAVARegressor.h` + `ISplineTRRRegressor.h` |
| `src/openms/thirdparty/percolator/PAVARegressor.h` | Create (sync + wrap) | PAVA extracted (header-only) |
| `src/openms/thirdparty/percolator/ISplineTRRRegressor.h` | Create (sync + wrap) | TRR I-spline regressor (~476 lines header-only) |
| `src/openms/thirdparty/percolator/IsotonicPEP.h` | Rewrite (sync + wrap) | Shrunk header — `InferPEP` class only |
| `src/openms/thirdparty/percolator/IsotonicPEP.cpp` | Rewrite (sync + wrap) | Factory-dispatch impl |
| `src/openms/thirdparty/percolator/Scores.h` | Modify (sync + wrap) | Adds `reserve(size_t)` inline + new `calcQvals` overload decl |
| `src/openms/thirdparty/percolator/Scores.cpp` | Modify (sync + wrap) | 3 context shifts inside existing namespace-wrap hunks |
| `src/openms/thirdparty/percolator/Globals.cpp` | Modify (sync) | Drop `getXMLDir` body (41 lines) |
| `src/openms/thirdparty/percolator/Globals.h` | Modify (manual) | Drop `getXMLDir` declaration (our static stub, not synced) |
| `src/openms/thirdparty/percolator/CrossValidation.cpp` | Modify (sync + wrap) | Cosmetic `std::make_pair` → brace-init (2 sites) |
| `src/openms/thirdparty/percolator/patches/01-namespace-wrap.patch` | Modify | Regenerate hunks for `IsotonicPEP.{h,cpp}`, `Scores.{h,cpp}`, and new files |
| `docs/superpowers/specs/2026-04-23-percolator-in-process-parity-design.md` | Modify | Record new observed drift values; note TODOs resolved |

**No expected edits to:**
- `src/openms/source/ANALYSIS/ID/Percolator.cpp` — wrapper API surface unchanged.
- `src/openms/source/ANALYSIS/ID/PercolatorImpl.h` — no internal structure changes required.
- `src/topp/PercolatorAdapter.cpp` — adapter CLI unchanged.
- Existing test files — `pep_method="nonparametric"` in §6 is unchanged algorithmically; other tests use score/q tolerances that absorb the PEP-calibration shift.

## Known risk areas

1. **P3.a / P3.b strict `==` bit-identical** (`Percolator_test.cpp`): the tests compare run1 vs run2, not against stored values, so the TRR I-spline path must be deterministic at `num_threads=1`. TRR is iterative and deterministic, but if the new path pulls in any thread-dependent state (OpenMP reductions inside Eigen dense ops), bit-identicality may break. Mitigation: if this happens, consider adding `-DEIGEN_DONT_PARALLELIZE` to the vendored source compile flags.
2. **P3.c Jaccard (0.99 floor)** and **P3.d count ratio (0.20 floor)** — may be tight if the PEP path shift propagates into q-values. Unlikely but possible.
3. **§2 / §7 / §8** (subprocess parity at `r ≥ 0.999` with strict count match): synthetic data, well-conditioned training. Score parity depends on SVM only; subprocess 3.06 unchanged; in-process SVM algorithm unchanged by `eb157f7`. Should pass unchanged.
4. **§6 (library realistic) and `PercolatorAdapter_parity_test`**: use `"nonparametric"` path which is unchanged. `max_dpep` tolerances (`0.15`, `0.25`) may tighten if residual 3.06↔3.08 PEP drift was specifically in the PAVA path (which isn't used by `nonparametric`).

---

## Task 1: Capture baseline test output

**Files:** None modified. Just record numeric observations for post-migration comparison.

**What this guards:** before touching the vendored code, lock in current behavior. Makes it easy to triage whether a post-migration test failure is due to the new code or a pre-existing test flakiness.

- [ ] **Step 1: Add baseline diagnostic prints to parity tests (temporary)**

Open `src/tests/class_tests/openms/source/Percolator_subprocess_parity_test.cpp` and `src/tests/class_tests/openms/source/PercolatorAdapter_parity_test.cpp`. At the final TEST_TRUE on `max_dpep` in each section that has one (§2, §6, §7, §8 in subprocess parity; the single section in adapter parity), add an `std::cerr` line immediately before the `TEST_TRUE(max_dpep <= ...)` call:

```cpp
    std::cerr << "BASELINE: section=<name> r=" << r << " max_dpep=" << max_dpep << std::endl;
```

Fill `<name>` with the section identifier (e.g. `§2`, `§6-realistic`, `§7-reservoir`, `§8-multifile`, `P2-adapter`). These are temporary — removed in Task 13.

- [ ] **Step 2: Rebuild and run all three suites**

```bash
cmake --build OpenMS-build --target Percolator_test Percolator_subprocess_parity_test PercolatorAdapter_parity_test -j$(nproc)
ctest --test-dir OpenMS-build -R '^(Percolator_test|Percolator_subprocess_parity_test|PercolatorAdapter_parity_test)$' --output-on-failure 2>&1 | tee /tmp/percolator_baseline.txt
```

Expected: 3/3 tests pass.

- [ ] **Step 3: Grep baseline values out**

```bash
grep -E "^BASELINE:" /tmp/percolator_baseline.txt > /tmp/percolator_baseline_values.txt
cat /tmp/percolator_baseline_values.txt
```

Expected: one line per gated section, e.g.:
```
BASELINE: section=§2 r=0.999999 max_dpep=0.0012
BASELINE: section=§6-realistic r=1.0 max_dpep=0.10073
BASELINE: section=§7-reservoir r=0.999998 max_dpep=0.0023
BASELINE: section=§8-multifile r=0.999999 max_dpep=0.0015
BASELINE: section=P2-adapter r=1.0 max_dpep=0.23
```

Save this file — it's the reference for Tasks 10-12. Do NOT commit the diagnostic lines.

- [ ] **Step 4: Revert the diagnostic `std::cerr` lines**

Use `git checkout -- src/tests/class_tests/openms/source/Percolator_subprocess_parity_test.cpp src/tests/class_tests/openms/source/PercolatorAdapter_parity_test.cpp` to drop the temporary additions. Verify `git status --short` shows a clean tree.

- [ ] **Step 5: No commit for Task 1**

This task collects data, not code changes. Proceed to Task 2.

---

## Task 2: Update vendored-tree metadata

**Files:**
- Modify: `src/openms/thirdparty/percolator/whitelist.txt`
- Modify: `src/openms/thirdparty/percolator/sources.cmake`
- Modify: `src/openms/thirdparty/percolator/UPSTREAM_COMMIT`
- Modify: `src/openms/thirdparty/percolator/UPSTREAM.md`
- Modify: `src/openms/thirdparty/percolator/sync-from-upstream.sh`

- [ ] **Step 1: Add 4 new files to whitelist**

Open `src/openms/thirdparty/percolator/whitelist.txt`. Current entries are alphabetical within upstream's `src/` tree root. Add these lines (alphabetical insertion):

```
ISplineTRRRegressor.h
MonotoneRegressor.cpp
MonotoneRegressor.h
PAVARegressor.h
```

Find the right insertion points by looking at the existing alphabetical order (e.g., `ISplineTRRRegressor.h` goes just before `IsotonicPEP.cpp`; `MonotoneRegressor.cpp/.h` goes after `MassHandler.h`; `PAVARegressor.h` goes after `PSMDescription.h`).

- [ ] **Step 2: Add the one .cpp to sources.cmake**

Open `src/openms/thirdparty/percolator/sources.cmake`. Find the list of compiled percolator source files (look for `CrossValidation.cpp`). Add `MonotoneRegressor.cpp` to the list (alphabetically after `MassHandler.cpp` or wherever matches existing sort order).

- [ ] **Step 3: Update UPSTREAM_COMMIT**

Replace the contents of `src/openms/thirdparty/percolator/UPSTREAM_COMMIT` with the new full SHA:

```
eb157f74e963430e559e0d0bcd31291e4ad660ba
```

(Single line, no trailing markers — the file currently has just a SHA and a newline.)

- [ ] **Step 4: Update UPSTREAM.md**

Open `src/openms/thirdparty/percolator/UPSTREAM.md`. Change the two metadata lines at the top:

```
- **Pinned at**: commit eb157f7 (post-rel-3-08-01; PR #399 / I-spline PEP)
- **Commit SHA**: eb157f74e963430e559e0d0bcd31291e4ad660ba
```

- [ ] **Step 5: Update sync script default**

Open `src/openms/thirdparty/percolator/sync-from-upstream.sh`. Find `UPSTREAM_REF="${1:-rel-3-06-1}"` (or whatever the current default is). Change to:

```bash
UPSTREAM_REF="${1:-eb157f74e963430e559e0d0bcd31291e4ad660ba}"
```

This makes the script default to re-syncing to this same commit until someone changes it again.

- [ ] **Step 6: Commit metadata-only change**

```bash
git add src/openms/thirdparty/percolator/whitelist.txt \
        src/openms/thirdparty/percolator/sources.cmake \
        src/openms/thirdparty/percolator/UPSTREAM_COMMIT \
        src/openms/thirdparty/percolator/UPSTREAM.md \
        src/openms/thirdparty/percolator/sync-from-upstream.sh
git commit -m "$(cat <<'EOF'
chore(percolator): prep whitelist + metadata for eb157f7 sync

- Add MonotoneRegressor.{h,cpp}, PAVARegressor.h, ISplineTRRRegressor.h
  to whitelist (new in upstream PR #399).
- Add MonotoneRegressor.cpp to sources.cmake.
- Update UPSTREAM_COMMIT / UPSTREAM.md / sync-from-upstream.sh default
  ref to eb157f7.

Sync of the actual files happens in the next commit.
EOF
)"
```

---

## Task 3: Sync vendored files from upstream eb157f7

**Files:** overwritten or created (via the sync script or manual `gh api`):
- Create: `MonotoneRegressor.h`, `MonotoneRegressor.cpp`, `PAVARegressor.h`, `ISplineTRRRegressor.h`
- Overwrite: `IsotonicPEP.h`, `IsotonicPEP.cpp`, `Scores.h`, `Scores.cpp`, `Globals.cpp`, `CrossValidation.cpp`

**Note:** `Globals.h.cmake` is NOT in our whitelist (we use a static stub at `Globals.h` instead). The sync script should skip it correctly.

- [ ] **Step 1: Run sync-from-upstream.sh**

```bash
cd src/openms/thirdparty/percolator
./sync-from-upstream.sh eb157f74e963430e559e0d0bcd31291e4ad660ba
cd - >/dev/null
```

Expected output:
- Clones upstream shallow at the given ref
- Copies whitelisted files (now including the 4 new ones)
- Attempts to apply `patches/*.patch`

Expected patch-application failures on:
- `01-namespace-wrap.patch` hunks for `IsotonicPEP.{h,cpp}` (full rewrite)
- `01-namespace-wrap.patch` hunks for `Scores.cpp` (context shift on `"\% FDR"`, `getScoreLabelPairs`, `weedOutRedundant`)
- `01-namespace-wrap.patch` hunk for `Scores.h` (new overload)

If the sync script fails midway (e.g., on patch failures), the files are still copied. Good — we regenerate the patch in Tasks 5-6.

- [ ] **Step 2: Verify sync outcome**

```bash
git status --short src/openms/thirdparty/percolator/
```

Expected: 10+ files modified/added:
- 4 new files (`MonotoneRegressor.h/.cpp`, `PAVARegressor.h`, `ISplineTRRRegressor.h`)
- 6+ modified files (`IsotonicPEP.h/.cpp`, `Scores.h/.cpp`, `Globals.cpp`, `CrossValidation.cpp`; possibly also `PosteriorEstimator.cpp` and `SetHandler.cpp` — those have upstream edits but no context-breaking changes, so their patches should apply cleanly)

- [ ] **Step 3: Spot-check the new files landed unwrapped**

```bash
grep -l "namespace OpenMS" src/openms/thirdparty/percolator/MonotoneRegressor.h \
                          src/openms/thirdparty/percolator/MonotoneRegressor.cpp \
                          src/openms/thirdparty/percolator/PAVARegressor.h \
                          src/openms/thirdparty/percolator/ISplineTRRRegressor.h
```

Expected: no output (none of the new files are yet namespace-wrapped). This is the baseline for Task 4.

- [ ] **Step 4: Sanity-check vendored `IsotonicPEP` structure**

```bash
grep -nE 'class (InferPEP|IsotonicRegression|PavaRegression|IsplineRegression|MonotoneRegressor)' \
    src/openms/thirdparty/percolator/IsotonicPEP.h \
    src/openms/thirdparty/percolator/MonotoneRegressor.h
```

Expected output:
- `IsotonicPEP.h`: only `class InferPEP` (no `IsotonicRegression`, `PavaRegression`, `IsplineRegression` anymore).
- `MonotoneRegressor.h`: `class MonotoneRegressor` abstract base.

- [ ] **Step 5: No commit for Task 3 yet**

The sync leaves a broken build (namespace-wrap unapplied on new files; hunks failed on rewritten files). Proceed directly to Task 4.

---

## Task 4: Apply namespace-wrap to the four new files

**Files:**
- Modify: `src/openms/thirdparty/percolator/MonotoneRegressor.h`
- Modify: `src/openms/thirdparty/percolator/MonotoneRegressor.cpp`
- Modify: `src/openms/thirdparty/percolator/PAVARegressor.h`
- Modify: `src/openms/thirdparty/percolator/ISplineTRRRegressor.h`

**Template:** every vendored file must have its body wrapped in `namespace OpenMS { namespace Internal { namespace Percolator { ... }}}`. Wrap goes **after** all `#include` directives and any `using Eigen::...` aliases (those must live inside the OpenMS namespace to avoid global pollution).

- [ ] **Step 1: Open MonotoneRegressor.h and apply the wrap**

Inspect the file structure:
```bash
head -40 src/openms/thirdparty/percolator/MonotoneRegressor.h
tail -10 src/openms/thirdparty/percolator/MonotoneRegressor.h
```

Identify the boundaries: after the last `#include` (or guard opener), before the first class/struct/function declaration. And before the trailing `#endif` of the include guard.

Insert immediately after the last `#include`:
```cpp

namespace OpenMS { namespace Internal { namespace Percolator {
```

Insert immediately before the trailing `#endif`:
```cpp
}}}  // namespace OpenMS::Internal::Percolator

```

- [ ] **Step 2: Same for MonotoneRegressor.cpp**

```bash
head -20 src/openms/thirdparty/percolator/MonotoneRegressor.cpp
tail -10 src/openms/thirdparty/percolator/MonotoneRegressor.cpp
```

The file has includes at the top, then `make_monotone_regressor()` factory function body. Insert the opening wrap after the last `#include`; insert the closing wrap at EOF (no include guard; just the end of the file).

If the `.cpp` has `using namespace std;` at file scope, that can stay outside the OpenMS namespace or go inside — prefer inside to match the other vendored `.cpp` files' pattern. Check a reference: `grep -A3 "namespace OpenMS" src/openms/thirdparty/percolator/Scores.cpp | head -10`.

- [ ] **Step 3: PAVARegressor.h — header-only with PAVA class definition**

```bash
head -25 src/openms/thirdparty/percolator/PAVARegressor.h
```

Check for any file-scope `using` directives (common in Eigen-using headers — e.g. `using Eigen::VectorXd;`). Those go **inside** the OpenMS namespace after the wrap opens.

Wrap pattern same as Step 1.

- [ ] **Step 4: ISplineTRRRegressor.h — the big one (~476 lines)**

```bash
head -40 src/openms/thirdparty/percolator/ISplineTRRRegressor.h
```

Expected: `#include <Eigen/Dense>`, possibly file-scope `using Eigen::MatrixXd;`, `using Eigen::VectorXd;`. The TRR regressor class + helpers follow.

Apply wrap the same way. Confirm that any file-scope `using Eigen::X;` lines are **inside** the OpenMS namespace. If they're outside (which would pollute the global namespace), move them inside.

- [ ] **Step 5: Verify all four files wrapped**

```bash
for f in src/openms/thirdparty/percolator/MonotoneRegressor.h \
         src/openms/thirdparty/percolator/MonotoneRegressor.cpp \
         src/openms/thirdparty/percolator/PAVARegressor.h \
         src/openms/thirdparty/percolator/ISplineTRRRegressor.h; do
  echo "=== $f ==="
  grep -c "namespace OpenMS { namespace Internal { namespace Percolator" "$f"
  grep -c "}}}.*namespace OpenMS::Internal::Percolator" "$f"
done
```

Expected: each file reports `1` and `1` (one opening wrap, one closing wrap).

- [ ] **Step 6: Build the vendored sources only**

```bash
cmake --build OpenMS-build --target OpenMS -j$(nproc) 2>&1 | tail -20
```

Expected: likely compile errors at this point due to `IsotonicPEP.{h,cpp}` and `Scores.{h,cpp}` still lacking fresh namespace-wrap (Task 5). Capture the error messages; they confirm which files need Task 5 attention.

- [ ] **Step 7: No commit yet**

The build is still broken until Task 5 regenerates the failed hunks. Proceed.

---

## Task 5: Regenerate namespace-wrap for IsotonicPEP.{h,cpp}

**Files:**
- Modify: `src/openms/thirdparty/percolator/IsotonicPEP.h` — full namespace wrap re-application
- Modify: `src/openms/thirdparty/percolator/IsotonicPEP.cpp` — full namespace wrap re-application
- Modify: `src/openms/thirdparty/percolator/patches/01-namespace-wrap.patch` — replace the two obsolete hunk blocks with regenerated ones

**Approach:** the upstream files are now in the tree (from Task 3), but unwrapped — they're raw upstream `src/IsotonicPEP.h/.cpp`. We re-apply the wrap the same way as Task 4, then regenerate the patch hunks so the next sync will work.

- [ ] **Step 1: Inspect the new unwrapped `IsotonicPEP.h`**

```bash
head -40 src/openms/thirdparty/percolator/IsotonicPEP.h
tail -10 src/openms/thirdparty/percolator/IsotonicPEP.h
```

Expected: shrunk to ~55-60 lines. `class InferPEP` decl only; forward-declares `MonotoneRegressor`; `#include <memory>` (for `unique_ptr`); maybe `#include "MonotoneRegressor.h"`.

- [ ] **Step 2: Apply namespace wrap to IsotonicPEP.h**

Same pattern as Task 4: after includes, before `class InferPEP`, open the wrap. Before the include guard's `#endif`, close the wrap.

- [ ] **Step 3: Inspect and wrap IsotonicPEP.cpp**

```bash
head -30 src/openms/thirdparty/percolator/IsotonicPEP.cpp
tail -10 src/openms/thirdparty/percolator/IsotonicPEP.cpp
```

Apply wrap: after `#include` block, before the `InferPEP` method definitions. Close before EOF.

- [ ] **Step 4: Build to confirm IsotonicPEP is correctly wrapped**

```bash
cmake --build OpenMS-build --target OpenMS -j$(nproc) 2>&1 | grep -E "IsotonicPEP|error:" | head -20
```

Expected: any remaining errors should be in `Scores.{h,cpp}` or `MonotoneRegressor.cpp` (if its `#include "IsotonicPEP.h"` references an unwrapped name). `IsotonicPEP` itself should compile cleanly.

- [ ] **Step 5: Regenerate the patch hunks**

The `01-namespace-wrap.patch` file still contains the OBSOLETE hunks for `IsotonicPEP.{h,cpp}`. To regenerate them correctly we need to diff the upstream-source (unwrapped) file against our now-wrapped vendored file.

```bash
# Get upstream (unwrapped) content for both files
mkdir -p /tmp/perc-eb157f7
gh api "repos/percolator/percolator/contents/src/IsotonicPEP.h?ref=eb157f7" --jq '.content' | base64 -d > /tmp/perc-eb157f7/IsotonicPEP.h
gh api "repos/percolator/percolator/contents/src/IsotonicPEP.cpp?ref=eb157f7" --jq '.content' | base64 -d > /tmp/perc-eb157f7/IsotonicPEP.cpp

# Diff upstream vs wrapped-in-tree, producing a git-style patch fragment
diff -u /tmp/perc-eb157f7/IsotonicPEP.h src/openms/thirdparty/percolator/IsotonicPEP.h > /tmp/IsotonicPEP_h.hunk.patch
diff -u /tmp/perc-eb157f7/IsotonicPEP.cpp src/openms/thirdparty/percolator/IsotonicPEP.cpp > /tmp/IsotonicPEP_cpp.hunk.patch
```

The two `.hunk.patch` files are the new hunks. Their headers won't match `01-namespace-wrap.patch`'s existing format (which uses `a/` and `b/` path prefixes), so normalize:

```bash
sed -i 's|/tmp/perc-eb157f7/|a/|; s|src/openms/thirdparty/percolator/|b/|' \
    /tmp/IsotonicPEP_h.hunk.patch /tmp/IsotonicPEP_cpp.hunk.patch
```

- [ ] **Step 6: Splice the new hunks into the patch**

Open `src/openms/thirdparty/percolator/patches/01-namespace-wrap.patch`. Use Grep to find the `diff --git a/IsotonicPEP.h b/IsotonicPEP.h` and `diff --git a/IsotonicPEP.cpp b/IsotonicPEP.cpp` markers — these start the obsolete hunks. Find the next `diff --git` marker to identify where each hunk ends.

Replace the obsolete block for `IsotonicPEP.h` with the contents of `/tmp/IsotonicPEP_h.hunk.patch`. Same for `IsotonicPEP.cpp`. Keep the rest of the patch file unchanged.

- [ ] **Step 7: Verify the regenerated patch applies cleanly against upstream**

```bash
# Re-sync and re-apply from scratch as a sanity check (safe — working tree is already in the target state)
cd src/openms/thirdparty/percolator
./sync-from-upstream.sh eb157f74e963430e559e0d0bcd31291e4ad660ba 2>&1 | tail -10
cd - >/dev/null
git diff --stat src/openms/thirdparty/percolator/IsotonicPEP.h src/openms/thirdparty/percolator/IsotonicPEP.cpp
```

Expected: `git diff --stat` shows ZERO changes if the patch applies cleanly to fresh upstream content. If the diff shows any changes, the regenerated hunks are wrong — redo Step 5-6.

- [ ] **Step 8: No commit for Task 5 yet**

Task 6 still needs to handle Scores.{h,cpp}. Continue.

---

## Task 6: Regenerate namespace-wrap for Scores.{h,cpp}

**Files:**
- Modify: `src/openms/thirdparty/percolator/Scores.h` — update in-tree wrapped version to match new upstream
- Modify: `src/openms/thirdparty/percolator/Scores.cpp` — update in-tree wrapped version to match new upstream
- Modify: `src/openms/thirdparty/percolator/patches/01-namespace-wrap.patch` — regenerate Scores hunks

**Approach:** these files weren't fully rewritten upstream; just context-shifted at three places (`Scores.cpp`: `\%` escape removed, `getScoreLabelPairs` body restructured, `weedOutRedundant` uses `[]` instead of `.at()`) plus one new overload declaration in `Scores.h`. Our namespace wrap still applies in principle — just the hunks' context lines are stale.

- [ ] **Step 1: Apply namespace wrap to the synced (but unwrapped) Scores.h / Scores.cpp**

If the sync script failed on these hunks during Task 3, the in-tree files are unwrapped upstream content. Apply the wrap the same way as Task 4: open after `#include`s (inside the include guard for `.h`), close before `#endif` / EOF.

- [ ] **Step 2: Confirm the two files build**

```bash
cmake --build OpenMS-build --target OpenMS -j$(nproc) 2>&1 | grep -E "Scores\.(h|cpp)|error:" | head -15
```

Expected: no errors in `Scores.h` / `Scores.cpp`. If MonotoneRegressor.cpp's inclusion of PAVARegressor.h produces errors about Eigen types not being in scope, re-check the Task 4 wrap for `PAVARegressor.h` / `ISplineTRRRegressor.h`.

- [ ] **Step 3: Regenerate Scores patch hunks**

```bash
gh api "repos/percolator/percolator/contents/src/Scores.h?ref=eb157f7" --jq '.content' | base64 -d > /tmp/perc-eb157f7/Scores.h
gh api "repos/percolator/percolator/contents/src/Scores.cpp?ref=eb157f7" --jq '.content' | base64 -d > /tmp/perc-eb157f7/Scores.cpp

diff -u /tmp/perc-eb157f7/Scores.h src/openms/thirdparty/percolator/Scores.h > /tmp/Scores_h.hunk.patch
diff -u /tmp/perc-eb157f7/Scores.cpp src/openms/thirdparty/percolator/Scores.cpp > /tmp/Scores_cpp.hunk.patch

sed -i 's|/tmp/perc-eb157f7/|a/|; s|src/openms/thirdparty/percolator/|b/|' \
    /tmp/Scores_h.hunk.patch /tmp/Scores_cpp.hunk.patch
```

- [ ] **Step 4: Splice Scores hunks into the patch**

In `src/openms/thirdparty/percolator/patches/01-namespace-wrap.patch`, locate the existing hunks for `Scores.h` and `Scores.cpp` (via `grep -n "diff --git a/Scores"` on the patch). Replace each obsolete hunk block with the corresponding `/tmp/Scores_*.hunk.patch` content.

- [ ] **Step 5: Also add hunks for the new files**

The four new files (MonotoneRegressor.{h,cpp}, PAVARegressor.h, ISplineTRRRegressor.h) now exist in the vendored tree wrapped (from Task 4) but have no entries in `01-namespace-wrap.patch`. Generate hunks:

```bash
for f in MonotoneRegressor.h MonotoneRegressor.cpp PAVARegressor.h ISplineTRRRegressor.h; do
  gh api "repos/percolator/percolator/contents/src/$f?ref=eb157f7" --jq '.content' | base64 -d > "/tmp/perc-eb157f7/$f"
  diff -u "/tmp/perc-eb157f7/$f" "src/openms/thirdparty/percolator/$f" > "/tmp/${f//[.\/]/_}.hunk.patch"
  sed -i "s|/tmp/perc-eb157f7/|a/|; s|src/openms/thirdparty/percolator/|b/|" "/tmp/${f//[.\/]/_}.hunk.patch"
done
```

Append each to `01-namespace-wrap.patch` (in alphabetical file order, maintaining the patch's existing sort).

- [ ] **Step 6: Full resync + reapply + diff check**

```bash
cd src/openms/thirdparty/percolator
./sync-from-upstream.sh eb157f74e963430e559e0d0bcd31291e4ad660ba 2>&1 | tail -10
cd - >/dev/null
git diff --stat src/openms/thirdparty/percolator/
```

Expected: only the 5 metadata files from Task 2 should diff (whitelist, sources.cmake, UPSTREAM_COMMIT, UPSTREAM.md, sync-from-upstream.sh), plus the 01-namespace-wrap.patch. The actual source files (wrapped) should match what's in-tree after the resync+patch-apply → zero diff there.

If any source file shows a diff, the patch still isn't perfectly regenerated. Go back to Steps 3-5 for that file.

- [ ] **Step 7: Run the full OpenMS library build**

```bash
cmake --build OpenMS-build --target OpenMS -j$(nproc) 2>&1 | tail -5
```

Expected: `Built target OpenMS`. No errors.

If `ISplineTRRRegressor.h`'s Eigen types cause compile errors elsewhere, it's likely because a `using Eigen::X;` at file scope ended up outside the OpenMS namespace. Re-check Task 4 Step 4 wrap placement.

- [ ] **Step 8: No commit for Task 6 yet**

Still need Task 7 (Globals.h stub edit) before committing. Continue.

---

## Task 7: Remove getXMLDir from our Globals.h stub

**Files:**
- Modify: `src/openms/thirdparty/percolator/Globals.h` (our manually-maintained static stub — NOT synced from upstream since upstream has `Globals.h.cmake` instead)

**What this fixes:** `Globals.cpp` (now synced from upstream at eb157f7) no longer defines `getXMLDir()`. Our static stub `Globals.h` still declares it. Without this fix, linking OpenMS fails because the declaration has no definition, or a downstream caller references it.

- [ ] **Step 1: Locate the declaration**

```bash
grep -n "getXMLDir" src/openms/thirdparty/percolator/Globals.h
```

Expected: one line, a declaration like `static std::string getXMLDir();` or similar.

- [ ] **Step 2: Remove the declaration + any includes that exist only for it**

Delete the line. Check for adjacent comments that apply only to `getXMLDir` and remove those too. Verify no include (e.g. `<string>`) is needed only for this declaration; if so leave it — it's probably needed elsewhere.

- [ ] **Step 3: Build to confirm linkage**

```bash
cmake --build OpenMS-build --target OpenMS -j$(nproc) 2>&1 | tail -5
```

Expected: `Built target OpenMS`. If any undefined-reference error mentions `getXMLDir`, some caller still references it — grep for it across the tree and remove those calls (should be zero callers since our wrapper doesn't use XML).

- [ ] **Step 4: Commit the full sync (Tasks 3-7 together)**

```bash
git add src/openms/thirdparty/percolator/
git commit -m "$(cat <<'EOF'
chore(percolator): sync vendored tree to upstream eb157f7

Brings in upstream PR #399 changes (merged 2026-04-23):

- New files: MonotoneRegressor.{h,cpp}, PAVARegressor.h,
  ISplineTRRRegressor.h — abstract base + factory + PAVA + trust-region
  reflective I-spline regressor.
- Rewrites: IsotonicPEP.{h,cpp} — slimmed header, factory-dispatch
  implementation. PEP calibration moves from inline PAVA + binned-LDLT
  I-spline to the new factored regressor hierarchy.
- Updates: Scores.{h,cpp} with reserve() / new calcQvals overload /
  reduced heap allocations in getScoreLabelPairs / % FDR text fix /
  [] indexing in weedOutRedundant. Globals.cpp loses getXMLDir body.
  CrossValidation.cpp gets brace-init cosmetics.
- Patch regen: 01-namespace-wrap.patch updated for all affected files.
- Stub edit: Globals.h (our static stub, not synced) drops the
  getXMLDir declaration to match the upstream removal.

Scores::calcPep(bool spline, bool interp, bool pava) signature is
preserved, so Percolator.cpp's dispatch keeps working unchanged.
The "nonparametric" pep_method path is unaffected (still calls
PosteriorEstimator::estimatePEP). The "isotonic" and
"logistic_regression" paths now produce numerically different PEPs
(pseudocount addition, TRR I-spline vs. binned-LDLT I-spline).
EOF
)"
```

---

## Task 8: Full clean rebuild

**Files:** none modified.

**What this validates:** no residual stale object files, all percolator sources compile together with the new regressor hierarchy.

- [ ] **Step 1: Clean rebuild of the library target**

```bash
cmake --build OpenMS-build --target OpenMS -j$(nproc) --clean-first 2>&1 | tail -10
```

Expected: clean build, ~5-10 minutes.

If errors appear:
- Missing symbols in `ISplineTRRRegressor.h` → check Eigen namespace wrap placement (Task 4 Step 4).
- Multiple definition / redefinition errors → a namespace wrap was applied twice. Grep for duplicate `namespace OpenMS { namespace Internal { namespace Percolator {` blocks in each new file.
- Undefined `getXMLDir` → Task 7 Step 2 didn't fully remove references.

- [ ] **Step 2: Build the test targets**

```bash
cmake --build OpenMS-build --target Percolator_test Percolator_subprocess_parity_test PercolatorAdapter_parity_test PercolatorAdapter -j$(nproc) 2>&1 | tail -5
```

Expected: all four targets link cleanly.

- [ ] **Step 3: No commit for Task 8**

Build verification only. Continue.

---

## Task 9: Run Percolator_test, assess impact

**Files:** potentially modify `src/tests/class_tests/openms/source/Percolator_test.cpp` if tests need adjusting.

**What this task catches:** algorithmic shift in `pep_method="logistic_regression"` (the library default) — now routes to the new TRR I-spline regressor instead of the old PavaRegression-based path. P3.a/b bit-identicality expects run1 == run2 (not comparison against stored values), so identicality itself should hold if TRR is deterministic. But `medianTargetBeatsDecoy_` in the basic test may shift. Run and observe.

- [ ] **Step 1: Run Percolator_test**

```bash
ctest --test-dir OpenMS-build -R '^Percolator_test$' --output-on-failure 2>&1 | tail -40
```

- [ ] **Step 2: Assess outcome**

**If 100% pass:** proceed to Step 5.

**If P3.a or P3.b fails** (`bit_identical != true`): the new TRR I-spline solver is non-deterministic at `num_threads=1` (unexpected; TRR is deterministic in principle). Investigate: check if `ISplineTRRRegressor.h` uses any Eigen parallel ops that aren't gated on a thread-count. If so, set `-DEIGEN_DONT_PARALLELIZE` explicitly on the vendored compile flags. As a last resort, loosen P3.a/b to `abs(δ) ≤ 1e-12` tolerance with an inline comment explaining the TRR solver's deterministic-but-not-bit-exact property.

**If the basic `rescore` test fails**: the SVM path would be broken. Unexpected — SVM code (`ssl.cpp`, `LogisticRegression.*`) wasn't touched by `eb157f7`. Likely a namespace-wrap error. Recheck Task 4 wrap placements.

- [ ] **Step 3 (conditional): If P3.a/b need tolerance loosening, edit inline**

Example for P3.a:
```cpp
// Was: if (out1.scores[i] != out2.scores[i]) bit_identical = false;
// Now:
if (std::abs(out1.scores[i] - out2.scores[i]) > 1e-12) bit_identical = false;
```
Same for `q_values` and `peps`. Add a comment explaining: *"After upstream eb157f7's TRR I-spline PEP path, PEPs are deterministic but not bit-exact across Eigen/TRR iterative solver internals. 1e-12 tolerance preserves the 'bit-identical-level' contract within numerical noise."*

- [ ] **Step 4 (conditional): Rerun Percolator_test until green**

```bash
cmake --build OpenMS-build --target Percolator_test -j$(nproc) && \
ctest --test-dir OpenMS-build -R '^Percolator_test$' --output-on-failure
```

- [ ] **Step 5: Commit the test file change if made**

If Step 3 was needed:
```bash
git add src/tests/class_tests/openms/source/Percolator_test.cpp
git commit -m "$(cat <<'EOF'
test(Percolator): loosen P3.a/b tolerance to 1e-12 after eb157f7 sync

Upstream eb157f7 replaces the PAVA + binned-LDLT I-spline PEP path with
a trust-region reflective I-spline regressor. TRR is deterministic but
the Eigen dense ops internally are not bit-exact across runs (summation
ordering in BLAS-like kernels). Loosen strict == to 1e-12 tolerance —
still far tighter than any meaningful scientific drift, still catches
persistent-state / global-leak regressions.
EOF
)"
```

If Step 3 was NOT needed, no commit — continue to Task 10.

---

## Task 10: Run subprocess parity test, compare vs baseline

**Files:** potentially modify `src/tests/class_tests/openms/source/Percolator_subprocess_parity_test.cpp` tolerance values.

**What this tests:** are post-migration observed values better, equal, or worse than the baseline captured in Task 1?

- [ ] **Step 1: Re-add BASELINE diagnostic prints (temporary)**

Open `src/tests/class_tests/openms/source/Percolator_subprocess_parity_test.cpp`. Add the same `std::cerr << "POSTMIG: section=..."` lines as Task 1 Step 1 (at each `max_dpep` assertion site), BUT use the prefix `POSTMIG:` instead of `BASELINE:`.

- [ ] **Step 2: Rebuild and run**

```bash
cmake --build OpenMS-build --target Percolator_subprocess_parity_test -j$(nproc)
ctest --test-dir OpenMS-build -R '^Percolator_subprocess_parity_test$' --output-on-failure 2>&1 | tee /tmp/percolator_postmig.txt
grep "^POSTMIG:" /tmp/percolator_postmig.txt > /tmp/percolator_postmig_values.txt
```

- [ ] **Step 3: Compare baseline vs post-migration**

```bash
paste /tmp/percolator_baseline_values.txt /tmp/percolator_postmig_values.txt
```

Expected:
- **§2, §3, §4, §5, §7, §8** (synthetic with `pep_method="isotonic"`): r values should stay at their baseline (score algorithm unchanged). `max_dpep` may change slightly either way due to the PAVA pseudocount addition.
- **§6 (`pep_method="nonparametric"`)**: `r` and `max_dpep` should be essentially unchanged (the nonparametric path is untouched by eb157f7).

- [ ] **Step 4: Tolerate or tighten**

For each section:
- If the post-mig observed value is within the current test tolerance → no action.
- If it exceeds the current tolerance (test now fails) → investigate and, if it's algorithmic drift not a regression, widen the tolerance with an inline comment explaining. Flag to controller if unsure.
- If it's strictly tighter (lower `max_dpep`, higher `r`) by ≥ 2× → TIGHTEN the test tolerance to catch future regressions. Update the `TODO(percolator-3.08-binary)` markers accordingly (if this migration closes the gap, remove the TODO; if it doesn't close it, keep the TODO).

- [ ] **Step 5: Remove POSTMIG diagnostic prints**

`git checkout -- src/tests/class_tests/openms/source/Percolator_subprocess_parity_test.cpp` to drop the `std::cerr` lines. If tolerance adjustments (Step 4) are needed, re-apply only those changes (via direct edit).

- [ ] **Step 6: Rebuild and verify final test passes**

```bash
cmake --build OpenMS-build --target Percolator_subprocess_parity_test -j$(nproc) && \
ctest --test-dir OpenMS-build -R '^Percolator_subprocess_parity_test$' --output-on-failure
```

Expected: `PASSED`.

- [ ] **Step 7: Commit tolerance changes (only if any were made)**

If tolerances were adjusted in Step 4:

```bash
git add src/tests/class_tests/openms/source/Percolator_subprocess_parity_test.cpp
git commit -m "$(cat <<'EOF'
test(Percolator): update parity tolerances after eb157f7 sync

Post-migration observed values on the CometAdapter_4_out fixture:
  §6 (nonparametric) max_dpep went from <baseline> to <postmig>
  §2/§7/§8 synthetic max_dpep went from <baseline> to <postmig>

Tightened tolerances where the observed values improved; loosened /
kept where they shifted in the other direction (algorithmic drift from
the TRR I-spline / PAVA-pseudocount shift, not regression).
TODO(percolator-3.08-binary) markers removed for sections now meeting
tighter bounds.
EOF
)"
```

If no tolerance changes were needed, no commit — continue to Task 11.

---

## Task 11: Run adapter parity test, compare vs baseline

**Files:** potentially modify `src/tests/class_tests/openms/source/PercolatorAdapter_parity_test.cpp` tolerance values.

**What this tests:** the same pattern as Task 10 but at the PercolatorAdapter layer, where subprocess post-processing amplifies any residual drift.

- [ ] **Step 1-6: Mirror Task 10's flow for the adapter test**

Apply the same pattern as Task 10 to `src/tests/class_tests/openms/source/PercolatorAdapter_parity_test.cpp` and its single section. Note: this test spawns the full `PercolatorAdapter` TOPP binary — make sure `OpenMS-build/bin/PercolatorAdapter` is rebuilt after the vendored tree sync.

- [ ] **Step 7: Commit tolerance changes (only if any were made)**

Similar to Task 10 Step 7, but adjust the commit message subject/body for adapter test.

---

## Task 12: Full-suite validation

**Files:** none modified.

- [ ] **Step 1: Clean build + ctest of all three**

```bash
cmake --build OpenMS-build -j$(nproc) --clean-first 2>&1 | tail -5
ctest --test-dir OpenMS-build -R '^(Percolator_test|Percolator_subprocess_parity_test|PercolatorAdapter_parity_test)$' --output-on-failure 2>&1 | tail -10
```

Expected: 3/3 tests pass.

- [ ] **Step 2: If anything fails, go back to Task 9/10/11**

Do not commit tolerance relaxations without surfacing to the controller.

- [ ] **Step 3: No commit (validation-only task).**

---

## Task 13: Update spec + plan docs with final observations

**Files:**
- Modify: `docs/superpowers/specs/2026-04-23-percolator-in-process-parity-design.md`
- Modify: `docs/superpowers/plans/2026-04-23-percolator-in-process-parity.md` (mark as done; link to this migration plan)

**What this preserves:** record the final observed drift values + what the migration achieved, so future readers can see the trajectory.

- [ ] **Step 1: Update the spec's tolerance table**

Open `docs/superpowers/specs/2026-04-23-percolator-in-process-parity-design.md`. In the tolerance table, annotate any entries whose tolerance changed post-migration. Add a "Post-eb157f7" column or footnotes.

- [ ] **Step 2: Update the upgrade-path section**

The section "Upgrade path when the bundled binary is updated to 3.08.01" should now be retitled or amended to reflect that we've done the VENDORED upgrade (but NOT the subprocess binary upgrade). The two are still decoupled.

- [ ] **Step 3: Cross-reference this migration plan in the earlier plan**

At the top of `docs/superpowers/plans/2026-04-23-percolator-in-process-parity.md`, add a `> **Follow-up migration:** See [2026-04-24-percolator-eb157f7-migration.md]` note.

- [ ] **Step 4: Commit**

```bash
git add docs/superpowers/specs/2026-04-23-percolator-in-process-parity-design.md \
        docs/superpowers/plans/2026-04-23-percolator-in-process-parity.md
git commit -m "$(cat <<'EOF'
docs(Percolator): record eb157f7 migration outcome in spec + plan

Post-migration observed drift values and any tolerance updates are
captured in the spec's tolerance table. Cross-reference the eb157f7
migration plan from the original parity plan.
EOF
)"
```

---

## Self-review

**Spec coverage check:**

| Analysis finding | Task(s) | Covered |
|---|---|---|
| Add 4 new files to whitelist | Task 2 | ✓ |
| Add MonotoneRegressor.cpp to sources.cmake | Task 2 | ✓ |
| Update UPSTREAM_COMMIT / UPSTREAM.md / sync script | Task 2 | ✓ |
| Sync rewritten files (IsotonicPEP, Scores, Globals, CrossValidation) | Task 3 | ✓ |
| Sync new files | Task 3 | ✓ |
| Apply namespace wrap to 4 new files | Task 4 | ✓ |
| Regenerate namespace wrap for IsotonicPEP.{h,cpp} | Task 5 | ✓ |
| Regenerate namespace wrap for Scores.{h,cpp} | Task 6 | ✓ |
| Add hunks for new files to 01-namespace-wrap.patch | Task 6 | ✓ |
| Remove getXMLDir from Globals.h stub | Task 7 | ✓ |
| Clean rebuild | Task 8 | ✓ |
| Run Percolator_test, handle P3.a/b if needed | Task 9 | ✓ |
| Run subprocess parity test, compare to baseline | Task 10 | ✓ |
| Run adapter parity test, compare to baseline | Task 11 | ✓ |
| Full-suite validation | Task 12 | ✓ |
| Record outcome in spec/plan docs | Task 13 | ✓ |
| Baseline capture (pre-migration) | Task 1 | ✓ |

No gaps.

**Placeholder scan:**
- No "TBD" / "implement later" strings.
- Conditional steps (Task 9 Step 3, Task 10 Step 4, Task 11) are written as explicit conditionals with concrete code paths per branch, not "TBD if test fails".
- Every code block shows actual content.

**Type / API consistency:**
- `sync-from-upstream.sh` argument is always the full SHA or `eb157f7...` — consistent across Tasks 2 and 3.
- File paths use the same `src/openms/thirdparty/percolator/` prefix everywhere.
- Patch-hunk regeneration commands use the same `/tmp/perc-eb157f7/` temp staging.

**Risk acknowledgments captured inline:**
- Task 9 Step 2 has an explicit branch for P3.a/b non-determinism under TRR I-spline.
- Task 10/11 Step 4 covers all three directions (tighter / equal / looser) for tolerance updates with explicit action per branch.

## Estimated effort

**1 day** — aligns with the analysis agent's estimate.

Breakdown:
- Tasks 1-2 (baseline + metadata): 30 minutes.
- Task 3 (sync): 15 minutes.
- Tasks 4-7 (wrap new files, regenerate hunks, stub edit): 2-3 hours — most of the concentration work.
- Task 8 (clean rebuild): 10-15 minutes wait.
- Tasks 9-11 (test passes + tolerance adjustment): 1-2 hours including iteration.
- Task 12 (validation): 10 minutes.
- Task 13 (doc update): 15 minutes.
- Slack + surprises: 2 hours.
