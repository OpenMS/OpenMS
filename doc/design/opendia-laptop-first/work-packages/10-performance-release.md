# WP-10: Laptop performance and release gates

- **Track:** output/integration.
- **Status:** proposed. **Owner / reviewer:** unassigned.
- **Parent epic:** [#10661](https://github.com/OpenMS/OpenMS/issues/10661). **Issue:** [#10671](https://github.com/OpenMS/OpenMS/issues/10671).
- **Shared policies:** [architecture](../architecture.md) and [plan](../README.md).

## Objective

Decide adoption from representative end-to-end quality, runtime, memory and disk evidence, not isolated component speedups.

## Scope

- Run pinned upstream/CI develop plus #10642 baselines against the opt-in workflow on ProteoBench Astral and diaPASEF.
- Measure repeated cold preparation and warm reuse separately, including prediction, indexing, rescoring and export.
- Exercise proposed 2/4 GiB settings and 2/4 workers on a declared CPU-only laptop, recording actual process/cgroup peaks and swap.
- Compare results and rescore at matched retained populations, with table cardinality, final bytes and actual peak disk including temporary/sidecar overlap.
- Evaluate single-run family science separately from unchanged existing multi-input/matrix regressions. Do not infer experiment-wide FDR validity.
- Keep default adoption blocked on failed or unmeasured gates, without retrospectively relaxing thresholds.

## Dependencies

- **Acceptance execution / completion:** WP-00 and WP-09.
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
- [ ] Cold/warm runtime, memory, swap, final bytes and actual peak disk meet the predeclared WP-00 gates, with rescoring and output costs separately reported.
- [ ] Known-site/mixture and held-out/entrapment checks meet separate identification, localization and quantitative-error margins at matched search scope.
- [ ] Tile/worker/budget variations preserve identities and tolerance-equal scores/quantities, while oversized, low-disk and interruption cases fail safely.
- [ ] Public CLI reopen/export/rescore and unchanged existing multi-input/matrix regressions pass from the tested build.
- [ ] The adoption record enumerates all failed/unmeasured gates and retains opt-in status unless the agreed release criteria pass.

**Integration gates:** WP-09 supplies public workflow evidence, and WP-10 owns the final representative quality/resource/performance and adoption decision.
