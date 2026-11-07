> **Archived 7 November 2026.** This laptop-first WP-00 to WP-10 plan is closed and superseded by the [OpenDIA rewrite master plan](../../../README.md). It is kept unchanged for provenance. Its native-input policy, read-only inputs and gate discipline carry over.

# WP-09: Early opt-in CLI and final integration

- **Track:** output/integration.
- **Status:** existing-path integration checks pass for the current engineering increment; combined-engine acceptance remains open. **Owner / reviewer:** unassigned.
- **Parent epic:** [#10661](https://github.com/OpenMS/OpenMS/issues/10661). **Issue:** [#10675](https://github.com/OpenMS/OpenMS/issues/10675).
- **Shared policies:** [architecture](../architecture.md) and [plan](../README.md).

## Objective

Keep a runnable opt-in path throughout development, then validate the combined bounded single-run family workflow.

## Current integration evidence, 4 November 2026

The [shared-evidence, actual-apex and scoring-context lifetime increments](../scoring-lifetime-20261104.md) pass all 178 affected CTests after building prerequisites, one/eight-thread scorer checks, actual native ownership/pressure integration and full native scientific-output comparisons. Core/MRM/featureXML/Parquet scheduling is unchanged. Final committed memory-change source and executable/three libraries match the measured candidate. The context-lifetime commit is scientifically and locally validated, but the current implementation branch is **not a standalone PR branch against `develop`**: its local `develop` comparison contains 266 files of stacked OpenDIA work. The commit is ready for scoped review once the owner selects a dependency-aware PR base or an agreed split/rebase. This is existing-path integration evidence, not combined sparse-family execution or a new installed-consumer replay. Whole-job admission, native Thermo science, reopen/rescore and complete failure/publication gates remain open. No PR or merge has been opened.

## Scope

- Start with v0 contracts and existing-engine adapters, adding WP-03/04/05/08 independently rather than waiting for a replacement engine.
- Wire the tool-private entry point and shared private implementation library, with registered private tests and public CLI regression commands.
- Include IM and families in interfaces from day one, but do not make completed new IM/family science a prerequisite for the first runnable engine.
- Add final retrieval/family science after its experiments, with explicit single-run scope and no experiment-wide FDR claim.
- Preserve existing defaults, multi-input and matrix behavior. Expose only results/rescore retention for v1.
- Integrate byte admission, free-space preflight, temporary-file rename, explicit failure and cleanup without a new quota or capability framework.

## Dependencies

- **Start / local skeleton merge:** WP-01 only.
- **Final integration completion:** WP-02 through WP-08. Early adapter merges do not depend on this completion gate.
- WP-10 receives reproducible public commands and measured artifacts after final integration.

## Local merge criteria

- [ ] An out-of-tree tool/private-test build succeeds, and an opt-in existing-engine CLI smoke test produces fresh nonempty output before new science is ready.
- [ ] Adapter switches reproduce existing-engine identities and tolerance-equal scores/quantities. New single-run science rejects multi-run requests explicitly, while compatibility-path defaults and multi-input/matrix regressions remain unchanged.
- [ ] Final-slice CLI tests cover non-IM, diaPASEF, discrete FAIMS, missing acquisition and known-site/ambiguous family cases with expected scientific fields.
- [ ] Public reopen/rescore commands reproduce the declared full-population single-run results, and results-only inputs reject unsupported rescoring.
- [ ] Slow-writer, oversized-work, disk-write and cancellation tests terminate explicitly without complete-looking partial artifacts or lost alternatives.
- [ ] Final-slice runs record commands, hashes, tolerances, process/cgroup memory and actual peak disk across worker/tile settings for WP-10.

**Integration gates:** WP-09 owns final public end-to-end integration, and WP-10 owns representative quality/resource/performance and adoption gates.
