# WP-00: Acceptance baseline and datasets

Repository context: `doc/design/opendia-laptop-first/README.md` and `doc/design/opendia-laptop-first/architecture.md`.

- **Status:** proposed, not implemented.
- **Owner:** unassigned. Suggested expertise: benchmarking and scientific validation.
- **Reviewer:** unassigned, scientific and build/release reviewers required.
- **Dependencies:** none. This enables WP-01 and supplies every package's acceptance inputs.
- **Suggested issue title:** `[OpenDIA][WP-00] Establish scientific, resource and output-size baselines`.
- **Parent epic:** OpenDIA laptop-first search, issue not created.
- **Issue:** not created.
- **Repository source:** `doc/design/opendia-laptop-first/work-packages/00-acceptance-baseline.md`.

## Goal

Make future correctness, memory and speed claims reproducible using a coherent local build and representative public workflows. Establish what constitutes acceptable scientific output before tuning for speed.

## Scope and non-goals

Own the dataset manifest, benchmark protocol, reference outputs and scientific thresholds. Do not design the new detection/scoring model or claim that small legacy tests validate it. Do not commit proprietary acquisitions or large raw files to the source tree.

## Code placement

Keep OpenDIA-specific test scenarios and assertions in registered `src/tests/topp/` and `src/tests/class_tests/openms/` tests, not in generic core production code. When testing new workflow code, use the proposed paired `src/openms/include/OpenMS/ANALYSIS/OPENDIA/` and `src/openms/source/ANALYSIS/OPENDIA/` interfaces. Reusable primitive tests belong with their core classes. Shared measurement helpers must not require OpenDIA-specific types. No new source directories are created by this planning issue.

## Output-size contract

Baseline the real OSW producer, not just a database's compressed download size. Inventory table row counts, indexes, repeated library metadata, transition/family fan-out, logical/allocated file bytes and sidecars. Measure the whole run's peak workspace, reusable caches and final outputs separately, including WAL/journals and staging/compaction overlap. Define separate final/scratch/cache quotas and a combined per-volume limit, with quota-exhaustion cases in the baseline manifest.

Agree numerical byte ceilings and reduction goals for `results`, `rescore` and selected `diagnostic` profiles with WP-08/WP-10 before acceptance. Compare capability-matched artifacts: compact results are not equivalent to a complete rescoring database. Record bytes per reported result/family and retained evidence row, output read/write time, and artifact capabilities. Use fixed raw data with increasing absent hypotheses and site alternatives to expose unnecessary replication. Profiles must preserve current-run scientific results, not make files smaller by skipping inference inputs.

## Deliverables and owned interfaces

- Versioned manifest: input checksums, license/access, acquisition mode, library composition, modifications, known truth, parameter files, expected outputs and build identity.
- A reproducible out-of-tree build and executable/shared-library provenance record. Resolve the missing configured Thermo/.NET native dependency and stale-binary failures described in the evidence report before relying on those workflows.
- Small distributable CTest fixtures plus a separately provisioned representative benchmark suite.
- A result schema covering commands, outputs, exit status, stage/total time, process/cgroup peak memory, disk use, candidates, identification quality, localization error and quantitative error.
- Agreed numerical tolerances, non-inferiority margins and speed acceptance policy, signed off before evaluating the new implementation.
- Per-profile footprint baselines, byte ceilings/reduction goals, peak per-filesystem occupancy and a capability matrix covering result reading, localization, rescoring and optional diagnostics.

## Dataset coverage

| Input | Purpose |
|---|---|
| Conventional non-IM DIA with nonzero signals | Extraction semantics, regression and absent mobility |
| Real mobility/PASEF acquisition | Geometry, missing scans and overlapping RT/IM components |
| FAIMS acquisition | Discrete strata without invalid continuous-axis merging |
| Known-site peptidoform mixtures | Discriminating ions, ambiguous and coexisting forms, localization error |
| Entrapment/held-out hypotheses | Retrieval recall and confidence calibration |
| Fixed raw run with scaled predicted libraries | Preparation/index cost, absent-hypothesis growth and candidate fan-out |
| Dense regions and enlarged modification families | Worst-case memory, overflow and approximation diagnostics |

At least one end-to-end dataset must combine mobility with peptidoform competition. Separate IM and IPF fixtures are insufficient. If truth data is unavailable, record the blocked quality gate rather than substituting row-count equality.

## PR-sized implementation steps

1. Restore a coherent build, rerun relevant existing TOPP tests in a clean output directory, and archive failures with executable provenance.
2. Add the small-fixture manifest and nonzero numerical reference cases, including the known first-centroid issue.
3. Provision representative data externally with checksums and retrieval instructions. Define truth annotations and search spaces.
4. Add benchmark orchestration and machine-readable result collection. Agree quality, speed and per-profile output-size gates with scientific owners, including real OSW footprint and finalization peaks.

## Acceptance checks

- [ ] Real OpenDIA/OpenSWATH/OpenSwathInfer producers and export consumers run from the same build. Never count comparison success if its producer failed or reused stale output.
- [ ] Re-running the baseline reproduces IDs and tolerance-bounded scores/quantities, or documents nondeterminism explicitly.
- [ ] Both optional MS1/MS2 IPF evidence modes are identified in provenance. Scientific evidence cannot be disabled silently for a speed comparison.
- [ ] Every explicit architecture requirement maps to an input, command and expected observable result. Unavailable truth or unsupported formats are marked blocked, not passed.
- [ ] Cold and warm runs distinguish library/index preparation and filesystem/cache state. Repeat runs and report distributions rather than a single fastest timing.
- [ ] Real OSW footprint and new per-profile totals include all files/sidecars and peak disk overlap. Numeric size goals and retained capabilities are fixed before declaring an improvement.
- [ ] Result reopening/export and supported rescore operations verify the capabilities used in size comparisons. No smaller-file claim relies on undisclosed evidence removal.
- [ ] Test and helper placement respects OpenDIA versus reusable-core ownership, with no core dependency on OpenDIA policy.

## Risks and handoff

Dataset access, licensing and missing ground truth may block acceptance, not implementation of unrelated contracts. Give each owner a small legal fixture and the representative protocol. Deliver the completed manifest to WP-10; subsequent packages append observed evidence rather than rewriting the baseline after seeing results.
