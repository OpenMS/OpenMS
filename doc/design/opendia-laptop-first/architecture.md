# OpenDIA architecture and shared invariants

[Master plan](README.md) · [Historical evidence](evidence.md)

**Proposed, 20 October 2026.** This is the shared contract, not an implemented engine. WP documents contain only package-specific work and merge checks.

## Placement and incremental integration

Follow the existing tool-private layout:

```text
src/topp/OpenDIA/
  CMakeLists.txt
  OpenDIA.cpp
  include/
  source/
  test/
```

Move the existing `src/topp/OpenDIA.cpp` without behavior changes first. Extract orchestration into a private static/OBJECT library linked by the tool and its tests, registered with `openms_add_tool_class_test`. New contracts are unstable v0, with no installed headers or pyOpenMS/API promise. Keep dependency direction tool -> private engine -> OpenMS core. Promote generic code only for a real second consumer. Fix genuine shared-core bugs directly in core.

Raw tiles, streamed libraries, shared extraction and compact outputs enter the current engine behind options, with identical/tolerance-equal result checks. Only retrieval and family inference require new scientific algorithms. WP-05 adapts existing mass/RT/IM calibration so data handling does not wait for WP-06.

## Units and contracts

The execution unit is an acquisition-local raw tile. The inference unit is a local competing peptidoform family. Neither a tile boundary nor a worker defines statistical confidence.

| v0 record | Required semantics |
|---|---|
| Raw tile/region | Exact RT, m/z, frame/scan and isolation geometry, optional mobility axis, bounded owning storage and leased views |
| Hypothesis store | Stable family/precursor IDs, sparse fragments, target/decoy identity, prediction and library provenance, paged cold metadata |
| Extraction query | Run/lane, acquisition selection, RT/mz/IM bounds and units, tolerance/edge/aggregation policy, isotope and calibration identity |
| Evidence block | Sparse RT/mobility support, exact coordinates and missingness masks, shared immutable evidence |
| Candidate family | Competing site forms, charge/adduct constraints, local region, shared/discriminating fragment incidence, alternatives/nulls and completeness state |
| Score batch | Separate identification/localization/quantification results, ambiguity and provenance, independent of tile lifetime |

Mobility and families are represented from the start. Freeze family/score semantics only after WP-07's reviewed null/decoy experiment. Initial contracts can merge before that review.

## Scientific invariants

- Keep MS1/MS2 axes and frame/scan isolation geometry. FAIMS voltages are discrete strata. Missing acquisition is not measured zero. Non-IM inputs use an absent mobility axis. Lookup bins never replace measurement coordinates.
- Share only identical complete extraction queries. Keep raw access for apex/isotope/precursor evidence. Freeze calibration per pass and invalidate dependent caches when it changes.
- Define tolerance, edge inclusion and aggregation mathematically. Fix shared tophat first-centroid omission and beyond-last double counting in a separate develop PR, with IM/non-IM edge tests. Do not introduce a mode that preserves these bugs.
- Preserve the current picker's complete assay domain or equivalent state for compatibility adapters. For new retrieval, define completeness and boundary reconciliation, use RT interval overlap and uncertainty-aware mobility, and retain a recovery path without requiring detectable MS1 or strong anchors. Apex ownership alone cannot restore truncated evidence.
- Score shared and site-discriminating evidence jointly without treating shared observations as independent copies. Preserve unresolved/coexisting forms. Avoid dense forms-by-fragments-by-RT-by-mobility tensors and precomputed pairwise UIS expansion. Review pruning/approximation as a scientific change.
- Keep required candidates, rejected targets, decoys, alternatives and distinguishing evidence until confidence is finalized. Validate target/decoy selection and held-out/entrapment calibration. Do not concatenate tile-level q-values. Separate identification, localization and protein confidence.
- New family science is single-run in v1. Existing multi-input and matrix export stay available unchanged through the compatibility path. Cross-run alignment, experiment-level FDR and merging new artifacts need a later explicit model and retained-population contract.

## Resource and output policy

Use one byte budget, bounded queues and a fixed stage acquisition order, with cancellation and backpressure. Count decoder/index/evidence/model/writer workspace and leave runtime/native headroom. External RSS/cgroup measurements, including children, establish whether the proposed 4 GiB/2 GiB process targets are met. Split, spill or fail oversized work explicitly, never silently truncate hypotheses. WP-02 hardening runs alongside component work.

| v1 output | Contents and limits |
|---|---|
| `results` (default) | Compact identities, localization/ambiguity, quantities, confidence and provenance. No full-library/raw copies or dense traces. Not an arbitrary rescoring input. |
| `rescore` (opt-in) | Full candidate/target/decoy/alternative population and evidence needed by explicitly named rescoring operations, including required calibration/model/library references. Test supported operations and reject unsupported ones. |

Use a schema version and document required dependencies. Retention does not change current-run scientific results. Derive retained columns/populations from actual consumers. Audit `FEATURE_TRANSITION` fan-out and repeated metadata before choosing a storage format. A Parquet substitution alone is not a size improvement, and compact results cannot be compared as capability-equivalent to full OSW evidence.

For v1, check free space, write temporary artifacts on the destination filesystem, publish by rename only after successful finalization, and report failed writes explicitly. Multi-file outputs need a final completion marker. Cleanup removes only owned intermediates, never inputs or retained dependencies. Measure actual final and peak disk use, including sidecars and publication overlap. Reuse library/index caches across compatible runs and report per-run output multiplied by run count. Per-class/per-volume quota frameworks, diagnostic profiles and richer capability manifests are deferred.

## Validation ownership

WP-00 pins upstream/CI build and data provenance and agrees quality/runtime/size thresholds. Use the ProteoBench Astral run from #10642 plus a diaPASEF dataset, not the tiny historical fixture, to prioritize work. Include preparation, calibration and rescoring/inference/export in timing.

Component WPs close on local tests and their incremental compatibility checks. WP-09 owns the public raw-to-export integration, actual reopening/rescoring and failure/cancellation tests. WP-10 owns representative cold/warm timing, complete process memory, final/peak disk, scientific non-inferiority and the default-adoption decision. Preserve candidate identities/sets and tolerance-bounded scores across tile sizes and worker counts. The final scientific slice must exercise mobility and peptidoform competition together, even though the first runnable increment need not.
