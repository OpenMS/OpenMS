# WP-00: Acceptance baseline and datasets

- **Track:** baseline.
- **Status:** proposed. **Owner / reviewer:** unassigned.
- **Parent epic:** [#10661](https://github.com/OpenMS/OpenMS/issues/10661). **Issue:** [#10662](https://github.com/OpenMS/OpenMS/issues/10662).
- **Shared policies:** [architecture](../architecture.md) and [plan](../README.md).

## Objective

Establish reproducible scientific, runtime and storage comparisons without blocking independent implementation.

## Scope

- Pin an upstream/CI `develop` revision and the applicable #10642 revision/state, commands, dependencies and input checksums.
- Use ProteoBench Astral plus diaPASEF as representative benchmarks, with known-site data for localization validation.
- Record cold preparation and warm reuse, stage/total runtime, peak memory, actual peak disk and output row cardinalities.
- Inventory OSW `FEATURE_TRANSITION` fanout and rescoring cost for WP-08, alongside identification, localization and quantitative error.
- The local 600pt fixture is a regression aid, not the baseline. Repairing local `libnethost` is not a prerequisite.

## Dependencies

- **Start / local merge:** none.
- Supply pinned inputs and thresholds to WP-10. WP-01 and the WP-07 early experiment can proceed independently.

## Local merge criteria

- [ ] A manifest records exact source/CI and #10642 provenance, checksums, parameters and executable commands for both representative datasets.
- [ ] Fresh successful producer logs and reference outputs are attached for Astral and diaPASEF, with no stale-output reuse.
- [ ] A repeat baseline run reproduces identities and scores/quantities within recorded numerical tolerances.
- [ ] The measurement table includes stage runtime and memory, process/cgroup peak memory, rescoring share, table row counts, final bytes and observed peak disk including sidecars/staging.
- [ ] Numerical quality, runtime and results/rescore size gates are recorded before candidate-engine evaluation, with unavailable truth explicitly blocked.

**Integration gates:** WP-09 owns public workflow checks, and WP-10 owns representative quality/resource/performance acceptance.
