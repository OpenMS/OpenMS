# Percolator in-process vs external-binary parity — test plan design

**Date:** 2026-04-23
**Scope:** OpenMS in-process Percolator (`Percolator::rescore/train/score`, vendored at `src/openms/thirdparty/percolator/`, pinned to upstream tag `rel-3-08-01`, SHA `febeef346327ff3adaf6712c7b8b250499aecc63`).

## Goal

Establish a test plan that gives high confidence that OpenMS's in-process Percolator behaves equivalently to the bundled `percolator` executable on the workloads OpenMS users actually run — in particular, on the paths exercised by `PercolatorAdapter`'s `-use_subprocess true/false` toggle.

## Baseline

- Reference binary: `THIRDPARTY/Linux/x86_64/Percolator/percolator`, version **3.06.0**. (Note: mismatch with vendored source version 3.08.01; see "Version note" below.)
- Reference invocation: `std::system()` subprocess, mirrors how `PercolatorAdapter` calls it today.
- All subprocess-dependent tests skip silently when the binary is not discoverable (`PERCOLATOR_BINARY` env var unset). This keeps CI green on systems without Percolator installed.

## Version note (3.06.0 baseline vs 3.08.01 vendored)

The bundled binary is 3.06.0, while the vendored in-process code is from 3.08.01. Where P1/P2 tests against realistic data see drift that's tighter than observed-against-3.06.0 tolerances once a 3.08.01 binary is swapped in, we'll tighten those tolerances in the next pass. Each relaxed tolerance in the test code carries an inline `TODO(percolator-3.08)` comment naming the expected tightening.

## Scope

### In scope
1. **Data diversity** at the library layer — realistic idXML inputs, larger-scale synthetic data, multi-file PIN structure.
2. **TOPP-adapter-layer parity** — side-by-side `PercolatorAdapter -use_subprocess true/false` on realistic idXML, score-by-score comparison of output idXMLs.
3. **Reproducibility / determinism** — invariance properties of the in-process path (same-instance consecutive, cross-instance same seed, input order, seed change as negative control).

### Out of scope (future work)
- Parameter-matrix expansion beyond the 5 flags already covered in `Percolator_subprocess_parity_test.cpp` §4.
- Scale / performance characterization (100K–500K PSMs, wall-clock, memory).
- Negative / boundary cases (too-few-decoys, constant feature, NaN/Inf) — assert same exception semantics both paths.
- Model persistence round-trip parity (`saveModel`/`loadModel` + `train()`-on-A / `score()`-on-B).
- Thread-count drift characterization (docstring explicitly flags it as expected).
- `.osw`-input adapter parity (`TOPP_PercolatorAdapter_2..4` territory).
- Peptide-level and protein-level output parity.
- Version-alignment track: rebuild a 3.08.01 binary from the pinned SHA, re-run the full plan, expect tighter tolerances.

## Three-pillar architecture

| Pillar | File | Layer | Comparison pair |
|---|---|---|---|
| P1 Data diversity | **extend** `src/tests/class_tests/openms/source/Percolator_subprocess_parity_test.cpp` | Library (`Percolator::rescore`) | in-process vs subprocess on hand-written .pin |
| P2 Adapter parity | **new** `src/tests/class_tests/openms/source/PercolatorAdapter_parity_test.cpp` | TOPP adapter | `PercolatorAdapter -use_subprocess true` vs `false` on same idXML |
| P3 Reproducibility | **extend** `src/tests/class_tests/openms/source/Percolator_test.cpp` | Library (self) | in-process vs itself across invariance axes |

## Pillar P1 — Data diversity (library layer)

Three new sections appended to `Percolator_subprocess_parity_test.cpp`. All gated on `PERCOLATOR_BINARY`.

### P1.a — Realistic-idXML library parity

Drive both paths with the real `PercolatorAdapter_1.idXML`:

1. Load idXML via `IdXMLFile`.
2. Build the numeric feature set (standard + `extra_features`), as the existing §1 (PIN stamp) test already does.
3. Build a `RescoreInput` out of meta values.
4. Call `Percolator::fillPINCompatibleFields(..., /*flatten_hits=*/true, ri)` so `scan` / `specFileNr` / `expMass` / `calcMass` match the PIN path.
5. Write the same data to a `.pin` via `PercolatorInfile::store`, run the subprocess.
6. Align subprocess rows with in-process rows by PIN `SpecId` string (also used as row id in subprocess output).
7. Compare: Pearson on scores, accepted-target set equality at q ∈ {0.01, 0.05, 0.10}, max `|Δpep|`.

Rationale: the most important addition — the real feature distribution users will hit, with real feature names / extra_features from PSMFeatureExtractor.

**Dataset caveat:** `PercolatorAdapter_1.idXML` has only 40 PeptideHits — known-good for the existing TOPP test but thin for statistical metrics (high variance in accepted-set comparisons at low q thresholds, few degrees of freedom for Pearson). If P1.a gets flaky or uninformative at this size, fallback is to fabricate a synthetic-but-PIN-compatible idXML (synthesize features + `fillPINCompatibleFields` + `IdXMLFile::store`) at ~1000 hits. The asset itself is still valuable because it carries the search-engine-specific feature names; scale is the only constraint.

### P1.b — Larger synthetic (reservoir-sampling path)

Extend `generateSyntheticData` with an optional `n_rows` argument; add a section that generates 20 000 rows and sets `subset_max_train = 5000` on both paths (`-N 5000` subprocess-side, `par.setValue("subset_max_train", 5000)` in-process).

Rationale: confirms the reservoir-sampled training subset agrees between paths. The existing §4 case `subset_max_train=200` is informative but below the natural dataset scale where the reservoir algorithm matters.

### P1.c — Multi-file PIN parity

Two sub-batches of 1000 synthetic rows each with distinct `spec_file_numbers` values. Extend `writePinFile` to optionally emit a `filename` column (triggering subprocess `specFileNr` population — see the comment at `generateSyntheticData` explaining why it's currently defaulted to 0). Populate `ri.spec_file_numbers` with matching values on the in-process side.

Rationale: confirms the CV fold hash (`OrderScanHash` over `specFileNr`+`scanNr`) agrees across paths when more than one file is present — a currently uncovered code path.

### Shared helper additions (inside anonymous namespace of `Percolator_subprocess_parity_test.cpp`)

- `writePinFile(path, ri, emit_filename=false)` — optional filename column.
- `buildRescoreInputFromIdXML(...)` — extracted from existing §1 helpers.

## Pillar P2 — Adapter-level parity (new file)

New file `src/tests/class_tests/openms/source/PercolatorAdapter_parity_test.cpp`. Wired into `src/tests/class_tests/openms/CMakeLists.txt` the same way `Percolator_subprocess_parity_test` is — skip silently when either binary is missing.

### Mechanism

1. CMake passes two env vars to the test executable:
   - `PERCOLATOR_BINARY` — already wired; used for `-percolator_executable …` on the subprocess path.
   - `PERCOLATOR_ADAPTER_BINARY` — new; points at the in-tree `bin/PercolatorAdapter` (CMake knows the path via `$<TARGET_FILE:PercolatorAdapter>`). Test skips silently if unset.
2. Test driver invokes `PercolatorAdapter` twice via `std::system(...)` on the same input idXML, identical flags except `-use_subprocess true` vs `false`. Input: `PercolatorAdapter_1.idXML` (same asset the TOPP test uses). Output paths from `File::getTemporaryFile()`.
3. Load both output idXMLs via `IdXMLFile`; extract (score, q-value, PEP) per hit. Row key: `(spectrum_reference, sequence, charge)` — reuses the `rowKey(...)` helper from P1.

### Prelude check — assert both paths stamp the same meta-key set

Iterate the union of meta keys that carry Percolator output on hits from each path. Either both stamp the `MS:1001491/492/493` CV keys, or both stamp `percolator_score/_q_value/_pep`, or the test fails with a readable message.

Rationale: prevents the parity test from silently degrading to "keys differed, but their absent values were zero, so everything was equal". Also turns the parity test into a spec for the adapter.

### Comparison

- Same surviving hit set on both sides (symmetric difference = ∅ after aligning by row key).
- Pearson r on SVM scores ≥ 0.99.
- Accepted-target set equality at q ∈ {0.01, 0.05, 0.10}.
- `max |Δpep| ≤ 0.05`.

### Deliberately not asserted

- No bit-equality on scores — adapter-level PEP method / CV-merge rescaling amplifies numerical drift from the library layer.
- No idXML byte-diff — would tangle the test with serializer ordering changes.

### Failure output

On failure, dump the first 5 mismatching rows with `(row_key, score_sub, score_inp, q_sub, q_inp, pep_sub, pep_inp)` so bisects are tractable.

### CMake sketch

```cmake
add_test(NAME PercolatorAdapter_parity_test ...)
if (PERCOLATOR_BINARY_FOR_TEST)
  set_tests_properties(PercolatorAdapter_parity_test PROPERTIES
    ENVIRONMENT "PERCOLATOR_BINARY=${PERCOLATOR_BINARY_FOR_TEST};PERCOLATOR_ADAPTER_BINARY=$<TARGET_FILE:PercolatorAdapter>")
endif()
```

## Pillar P3 — Reproducibility (extend class unit test)

New sections appended to `src/tests/class_tests/openms/source/Percolator_test.cpp`. No `PERCOLATOR_BINARY` dependency — these are in-process-only invariance tests. All four use a shared moderate-signal helper (reuse existing `makeModeratelySeparableInput_` or a variant at ~1000 rows for speed).

### P3.a — Same-instance, consecutive calls, bit-identical

Call `p.rescore(input)` twice on the same `Percolator` instance. Assert scores, q-values, and PEPs are bit-identical (`==`, not `abs(δ) < eps`). Guards against a persistent-RNG / leaked-global-state regression where a second call diverges from the first.

### P3.b — Cross-instance, same seed, bit-identical

Tightens the existing `[EXTRA] reproducibility with fixed seed` section: current code checks scores within `1e-9`. Upgrade to exact equality on scores **and** extend the assertion to q-values and PEPs. If there's silent drift hidden under the 1e-9 tolerance, this surfaces it.

### P3.c — Input-order invariance characterization

Shuffle the input rows under a fixed seed; run both permutations. The docstring explicitly says results depend on input ordering, so this is a **characterization** test, not an equivalence test:
- Pearson r ≥ 0.99 on scores between the two permutations (after un-permuting to align the same row in both runs).
- Accepted-target sets A, B at q ≤ 0.01 satisfy `|A ∩ B| / |A ∪ B| ≥ 0.95` (Jaccard).

Rationale: documents the observed drift so a regression that increases it (e.g., Jaccard dropping from 0.95 to 0.6) is caught. Jaccard (set-level) is preferred over count delta because counts can match while membership shuffles.

### P3.d — Seed change (negative control)

Run with seed 17 and seed 42. Assert:
- Scores are **not** bit-identical across the two seeds (catches a bug where `seed` param isn't wired through).
- Pearson r ≥ 0.9 across the two seeds (both converge to similar decision boundaries).
- Accepted-target counts at q ≤ 0.01 satisfy `|cA − cB| / max(cA, cB) ≤ 0.20` (counts, not sets — different seeds legitimately shuffle members near the threshold).

## Tolerance table (single source of truth)

| Test family | Layer | Scores | q-values | PEPs | Weights |
|---|---|---|---|---|---|
| P1 existing §2 (synth, defaults) | library | Pearson r ≥ 0.999; max rel `|Δ|` ≤ 1e-4 | exact target-count at q ≤ 0.01/0.05 | — | — |
| P1 existing §3 (synth, ranking) | library | — | accepted-set equality at q ≤ 0.01/0.05/0.10 | — | — |
| P1 existing §4 (param matrix) | library | Pearson r ≥ 0.99 | exact target-count at q ≤ 0.01 | — | — |
| P1 existing §5 (weights) | library | — | — | — | max `|Δw_feature|` ≤ 1e-3, bias excluded |
| **P1.a new (real idXML)** | library | Pearson r ≥ 0.99 | accepted-set equality at q ≤ 0.01/0.05/0.10 | max `|Δ|` ≤ 0.05 | — |
| **P1.b new (larger, reservoir)** | library | Pearson r ≥ 0.999 | exact target-count at q ≤ 0.01/0.05 | — | — |
| **P1.c new (multi-file PIN)** | library | Pearson r ≥ 0.999 | exact target-count at q ≤ 0.01/0.05 | — | — |
| **P2 adapter parity** | adapter | Pearson r ≥ 0.99 | accepted-set equality at q ≤ 0.01/0.05/0.10 | max `|Δ|` ≤ 0.05 | — |
| **P3.a / P3.b reproducibility** | in-proc | bit-identical | bit-identical | bit-identical | — |
| **P3.c input-order** | in-proc | Pearson r ≥ 0.99 | Jaccard ≥ 0.95 at q ≤ 0.01 | — | — |
| **P3.d seed change** | in-proc | NOT equal; Pearson r ≥ 0.9 | count ratio `|cA−cB|/max ≤ 0.20` at q ≤ 0.01 | — | — |

P1.a uses the looser `r ≥ 0.99` band because 3.06↔3.08 algorithmic drift is expected on real feature distributions. If observed drift is tighter in practice, we tighten the threshold in the next iteration.

### Post-eb157f7 migration (2026-04-24)

The vendored tree was bumped from `rel-3-08-01` (febeef3) to upstream eb157f7
(PR #399 — TRR I-spline PEP calibration); see
[`2026-04-24-percolator-eb157f7-migration.md`](../plans/2026-04-24-percolator-eb157f7-migration.md).
Post-migration observed values on the tolerance-relevant sections:

| Section | r (observed) | `max |Δ|` (observed) | Tolerance status |
|---|---|---|---|
| P1 existing §2 | 1.0 | max_rel_ds = 4.54e-6 | unchanged |
| **P1.a new (real idXML)** | 1.0 | max_dpep = 0.10073 | unchanged; TODO(percolator-3.08-binary) still open |
| **P1.b new (reservoir)** | 0.999998 | — | unchanged |
| **P1.c new (multi-file)** | 1.0 | — | unchanged |
| **P2 adapter parity** | 1.0 | max_dpep = 0.10073 | unchanged; TODO(percolator-3.08-binary) still open |

All values bit-identical to the pre-migration baseline — upstream eb157f7 does
not affect the `pep_method="nonparametric"` path the tests exercise. The
`TODO(percolator-3.08-binary)` markers remain valid: they track the
independent goal of upgrading the bundled subprocess `percolator` binary
from 3.06 to 3.08.01, which is a separate task from the vendored-tree sync.

## Execution model

- All P1 additions and P2 are CTest targets; gated on their respective env vars; green on systems without the binaries.
- P3 runs unconditionally (class test, no external dep).
- Enable locally: `ctest --test-dir OpenMS-build -R 'Percolator.*parity|Percolator_test'`.
- CMake wiring: extend the existing `find_program(percolator ...)` branch in `src/tests/class_tests/openms/CMakeLists.txt` to also set `PERCOLATOR_ADAPTER_BINARY=$<TARGET_FILE:PercolatorAdapter>` for the new adapter-parity test.

## CI considerations

OpenMS's existing CI either has or doesn't have `percolator` on PATH; this plan doesn't change that policy. If CI has the binary, all gated tests run. If not, only the always-run subset runs — P1 §1 (PIN stamp) plus P3 (reproducibility) — matching current behaviour.

## Upgrade path when the bundled binary is updated to 3.08.01

The relaxed P1.a tolerance (`Pearson r ≥ 0.99`) and the existing P1 §4 param-matrix `min_r = 0.99` floors become tightenable candidates. Each relaxed tolerance carries a `TODO(percolator-3.08)` comment inline, so the next sync from upstream knows where to re-check. The rest of the plan is version-agnostic.

## Deliverables summary

- Extensions to `src/tests/class_tests/openms/source/Percolator_subprocess_parity_test.cpp` (P1.a, P1.b, P1.c, plus `writePinFile` and `buildRescoreInputFromIdXML` helpers).
- New file `src/tests/class_tests/openms/source/PercolatorAdapter_parity_test.cpp` (P2).
- Extensions to `src/tests/class_tests/openms/source/Percolator_test.cpp` (P3.a–d).
- CMake wiring changes in `src/tests/class_tests/openms/CMakeLists.txt` to register `PercolatorAdapter_parity_test` and pass both `PERCOLATOR_BINARY` and `PERCOLATOR_ADAPTER_BINARY` to it.
- No changes to production code (`Percolator.cpp`, `PercolatorAdapter.cpp`, vendored `src/openms/thirdparty/percolator/`) are part of this plan. If the P2 "same meta-key set" prelude fails, that's a separate finding to be tracked as a follow-up.
