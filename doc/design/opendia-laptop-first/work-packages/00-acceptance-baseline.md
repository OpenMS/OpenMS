# WP-00: Acceptance baseline and datasets

[Plan index](../README.md) | [Architecture](../architecture.md) | [Prior evidence](../evidence.md)

- **Status:** proposed, not implemented.
- **Owner:** unassigned. Suggested expertise: benchmarking and scientific validation.
- **Reviewer:** unassigned, scientific and build/release reviewers required.
- **Dependencies:** none. This enables [WP-01](01-contracts.md) and supplies every package's acceptance inputs.

## Goal

Make future correctness, memory and speed claims reproducible using a coherent local build and representative public workflows. Establish what constitutes acceptable scientific output before tuning for speed.

## Scope and non-goals

Own the dataset manifest, benchmark protocol, reference outputs and scientific thresholds. Do not design the new detection/scoring model or claim that small legacy tests validate it. Do not commit proprietary acquisitions or large raw files to the source tree.

## Deliverables and owned interfaces

- Versioned manifest: input checksums, license/access, acquisition mode, library composition, modifications, known truth, parameter files, expected outputs and build identity.
- A reproducible out-of-tree build and executable/shared-library provenance record. Resolve the missing configured Thermo/.NET native dependency and stale-binary failures described in the evidence report before relying on those workflows.
- Small distributable CTest fixtures plus a separately provisioned representative benchmark suite.
- A result schema covering commands, outputs, exit status, stage/total time, process/cgroup peak memory, disk use, candidates, identification quality, localization error and quantitative error.
- Agreed numerical tolerances, non-inferiority margins and speed acceptance policy, signed off before evaluating the new implementation.

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
4. Add benchmark orchestration and machine-readable result collection. Agree quality and speed gates with scientific owners.

## Acceptance checks

- Real OpenDIA/OpenSWATH/OpenSwathInfer producers and export consumers run from the same build. Never count comparison success if its producer failed or reused stale output.
- Re-running the baseline reproduces IDs and tolerance-bounded scores/quantities, or documents nondeterminism explicitly.
- Both optional MS1/MS2 IPF evidence modes are identified in provenance. Scientific evidence cannot be disabled silently for a speed comparison.
- Every explicit architecture requirement maps to an input, command and expected observable result. Unavailable truth or unsupported formats are marked blocked, not passed.
- Cold and warm runs distinguish library/index preparation and filesystem/cache state. Repeat runs and report distributions rather than a single fastest timing.

## Risks and handoff

Dataset access, licensing and missing ground truth may block acceptance, not implementation of unrelated contracts. Give each owner a small legal fixture and the representative protocol. Deliver the completed manifest to [WP-10](10-performance-release.md); subsequent packages append observed evidence rather than rewriting the baseline after seeing results.
