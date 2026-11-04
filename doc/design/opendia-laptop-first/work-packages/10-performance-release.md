# WP-10: Laptop performance and release gates

- **Track:** output/integration.
- **Status:** separate native-fixture engineering speed and memory gates passed; representative scientific/scale/release acceptance remains open. **Owner / reviewer:** unassigned.
- **Parent epic:** [#10661](https://github.com/OpenMS/OpenMS/issues/10661). **Issue:** [#10671](https://github.com/OpenMS/OpenMS/issues/10671).
- **Shared policies:** [architecture](../architecture.md) and [plan](../README.md).

## Objective

Decide adoption from representative end-to-end quality, runtime, memory and disk evidence, not isolated component speedups.

## Current measured increment, 4 November 2026

The [actual-apex scheduler](../scoring-lifetime-20261104.md#measured-native-fixture-outcome) passes two predeclared counterbalanced native TimSim pairs: median wall 116.395 to 40.040 s, unchanged scientific outputs, scoring decodes about 86.6% lower, with 17.12% higher median process peak RSS. The subsequent [memory follow-through](../scoring-lifetime-20261104.md#scoring-context-memory-follow-through), measured against the retained apex scheduler in a separate session, passes its prospective gate: median RSS 1,688,332 to 1,505,442 KiB (10.83% lower), lower RSS in each pair, median wall 38.880 to 38.795 s and all scientific checks passing. This establishes no new runtime speedup or combined earlier-parent result. It closes useful engineering gates, not this package's release gate. Astral/native Thermo, large predicted libraries, cold/warm complete process/child memory/disk and competitive scientific quality remain explicitly unverified. No default adoption or package closure follows from this fixture.

## Scope

- Run pinned upstream/CI develop plus #10642 baselines against the opt-in workflow on ProteoBench Astral and diaPASEF.
- Measure repeated cold preparation and warm reuse separately, including prediction, indexing, rescoring and export.
- Exercise proposed 2/4 GiB settings and 2/4 workers on a declared CPU-only laptop, recording actual process/cgroup peaks and swap.
- Compare results and rescore at matched retained populations, with table cardinality, final bytes and actual peak disk including temporary/sidecar overlap.
- Evaluate single-run family science separately from unchanged existing multi-input/matrix regressions. Do not infer experiment-wide FDR validity.
- Keep default adoption blocked on failed or unmeasured gates, without retrospectively relaxing thresholds.
- Compare pinned available DIA-NN and Spectronaut versions under the [quality-matched comparison protocol](../performance-strategy.md#comparable-to-dia-nn-and-spectronaut-means-measured-science). Match search/library/MBR and inference scope, externally check error rates and assess quantitative accuracy/completeness, not just ID counts or equal nominal q-values. Unavailable licensed comparators remain explicitly unrun without blocking local OpenDIA prototyping.
- Evaluate behavior-preserving engine rewrites by reference output equivalence, and explicitly new scientific scorers by prospectively declared non-inferiority margins and independent truth/error-control checks. Include native decoding, optional indexing, prediction/training, scoring, rescoring and export costs. No numerical scientific margin or comparator parity is approved or established by this plan update.

## Dependencies

- **Acceptance execution / completion:** completed WP-00 baseline plus WP-09 integration. WP-10 agrees prospective candidate thresholds before acceptance runs. The current baseline need not meet rewrite resource targets.
- Benchmark tooling may be prepared earlier, but neither tiny fixtures nor component-only timing closes these gates.

## Local merge criteria

The [WP-00 owner/reviewer decision matrix](00-acceptance-baseline.md#ownerreviewer-decision-matrix)
records proposed future-candidate gates: behavior-preserving refactors retain
exact accepted sets and frozen-policy numeric results without added missing or
nonfinite cells; matched same-host median wall time, final results plus rescore
bytes, complete-process/cgroup peak memory and actual peak disk do not exceed
their matched baseline. These are proposals pending owner/reviewer agreement,
not accepted numerical gates or measured improvements. Agree repetitions, cold/
warm configuration, artifact capabilities and retained populations before
candidate evaluation. No new numerical truth-accuracy/FDP cutoff is specified;
report observed metrics and keep missing known-site localization truth blocked.
Selected merged-source workstation verification is in progress, not release
acceptance. Complete adversarial review before any PR.

- [ ] A machine-readable matrix records exact commands, source/input hashes, hardware, repetitions and pass/fail/blocked results for Astral and diaPASEF.
- [ ] Cold/warm runtime, memory, swap, final bytes and actual peak disk meet prospective candidate gates agreed by WP-10 using the WP-00 baseline, with rescoring and output costs separately reported.
- [ ] Known-site/mixture and held-out/entrapment checks meet separate identification, localization and quantitative-error margins at matched search scope.
- [ ] Tile/worker/budget variations preserve identities and tolerance-equal scores/quantities, while oversized, low-disk and interruption cases fail safely.
- [ ] Public CLI reopen/export/rescore and unchanged existing multi-input/matrix regressions pass from the tested build.
- [ ] The adoption record enumerates all failed/unmeasured gates and retains opt-in status unless the agreed release criteria pass.

**Integration gates:** WP-09 supplies public workflow evidence, and WP-10 owns the final representative quality/resource/performance and adoption decision.
