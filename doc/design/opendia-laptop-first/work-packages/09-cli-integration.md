# WP-09: Early opt-in CLI and final integration

- **Track:** output/integration.
- **Status:** proposed. **Owner / reviewer:** unassigned.
- **Parent epic:** [#10661](https://github.com/OpenMS/OpenMS/issues/10661). **Issue:** [#10675](https://github.com/OpenMS/OpenMS/issues/10675).
- **Shared policies:** [architecture](../architecture.md) and [plan](../README.md).

## Objective

Keep a runnable opt-in path throughout development, then validate the combined bounded single-run family workflow.

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
