# WP-10: Laptop performance and release gates

[Plan index](../README.md) | [Architecture](../architecture.md)

- **Status:** proposed, not implemented.
- **Owner:** unassigned. Suggested expertise: CPU performance, benchmarking and QA/release.
- **Reviewer:** unassigned, with scientific and release-owner approval.
- **Dependencies:** [WP-09](09-cli-integration.md), using [WP-00](00-acceptance-baseline.md) data and pre-agreed gates. Benchmark infrastructure can start earlier.

## Goal

Demonstrate that the complete new public workflow is genuinely useful on a resource-limited laptop, without hiding preparation costs, reducing scientific scope or relying on unit-test proxies.

## Scope and non-goals

Own acceptance execution, profiling-guided tuning and the adoption recommendation. Do not claim a universal speedup, replace defaults before quality gates pass, or silently alter retrieval/scoring settings to win benchmarks. Algorithmic changes found necessary return to their owning package for scientific review.

## Deliverables and owned interfaces

- Reproducible benchmark recipes for declared CPU, RAM, storage, OS, power mode, thread counts, executable/dependency versions and dataset checksums.
- Per-requirement results with commands, actual outputs and pass/fail/blocked status, plus raw logs and summary tables.
- Stage and total wall/CPU time, RSS/cgroup peaks, page faults/swap, raw rereads, queue high-water marks, temporary/final disk use and search coverage.
- Quality comparison: pre-FDR recall, sequence/peptidoform/protein confidence, localization error, unresolved/coexisting forms and quantitative error.
- A release decision recording passed gates, unsupported cases, known limits and whether the opt-in workflow is ready for wider use or default adoption.

## Benchmark matrix

| Axis | Required settings |
|---|---|
| Process memory | Proposed 2 GiB and 4 GiB settings on a CPU-only 16 GiB laptop, with external measurement |
| Parallelism/locality | One-thread reference, 2 and 4 workers, multiple tile/queue sizes |
| Cache state | Cold library/index preparation plus run; warm reuse reported separately |
| Acquisition | Non-IM, mobility/PASEF, FAIMS, and a combined mobility/peptidoform truth workflow |
| Library growth | Fixed acquisition with increasing absent hypotheses, fragments and localization alternatives |
| Difficult evidence | Weak/no MS1, weak distinguishing fragments, overlap, uncertain calibration and tile-boundary peaks |
| Resource failure | Oversized region/family, slow writer, low disk, low budget, interrupted processing and corrupt checkpoint |

Repeat and interleave baseline/new runs where practical. Report cache-state preparation rather than calling a warm filesystem run cold. Include model loading, prediction, index construction, confidence estimation and export in first-use cost.

## PR-sized implementation steps

1. Finalize numeric speed and quality thresholds with WP-00 before viewing candidate-engine performance. Register small end-to-end regression recipes.
2. Run the full public workflow matrix and profile the dominant measured costs. Assign targeted optimizations to package owners, one hypothesis per change.
3. Re-run matched-policy comparisons after each optimization, checking full outputs and scientific quality as well as timing/memory.
4. Validate installed package/public consumers and user documentation, then publish the evidence-backed release recommendation.

## Acceptance checks

- The end-to-end process meets the agreed 2/4 GiB limits without swapping. Accounting totals alone are not accepted as RSS evidence. Record process-tree/cgroup behavior and native allocations.
- Total cold and warm runtime passes the pre-agreed speed gate on representative data at matched hypotheses and scientific quality. A tiny extraction-only pool benchmark is not sufficient.
- Known-site/mixture results pass localization and quantitative quality margins independently of sequence identification. Targets, decoys and entrapment retain valid global confidence behavior.
- Candidate identities/sets and tolerance-bounded scores/quantities are stable across tile and worker configurations. Any approximation is declared and evaluated separately.
- Scaling and adversarial tests remain bounded or fail explicitly without complete-looking partial output. Disk spill, cancellation and restart do not change successful final results.
- Real installed CLI/export workflows and any changed public C++/Python interfaces pass their relevant consumers. Unit/header/package tests supplement these paths.

## Risks and handoff

Laptop power/thermal variation, dataset bias and cache state can dominate timings. Report repetitions and uncertainty, not only best cases. If no improvement or a quality regression remains, retain the opt-in implementation and report the failed gate rather than changing success criteria retrospectively. Hand the release owner the benchmark manifest, per-requirement evidence, support matrix and an explicit default-adoption decision.
