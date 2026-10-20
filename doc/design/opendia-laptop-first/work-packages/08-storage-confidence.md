# WP-08: Compact storage and single-run confidence

- **Track:** output/integration.
- **Status:** proposed. **Owner / reviewer:** unassigned.
- **Parent epic:** [#10661](https://github.com/OpenMS/OpenMS/issues/10661). **Issue:** [#10670](https://github.com/OpenMS/OpenMS/issues/10670).
- **Shared policies:** [architecture](../architecture.md) and [plan](../README.md).

## Objective

Reduce existing output fanout and rescoring cost early, preserving the complete population needed for confidence.

## Scope

- Prioritize existing `FEATURE_TRANSITION` fanout and the reported roughly 40% rescoring runtime share, remeasuring both on pinned inputs.
- Plug storage/rescoring adapters opt-in into the existing engine before new family science. Choose format from measured row cardinality, bytes and read/write cost, not format replacement alone.
- v1 retains only default `results` and optional `rescore`, with schema versions. Defer diagnostic profiles, general quota frameworks and capability manifests beyond schema version.
- Preserve full target/decoy/alternative populations until single-run confidence completes, with distinct identification/localization/protein confidence and no tilewise q-value concatenation.
- Check free space upfront, write temporary files then rename on success, report write failures and clean task-owned temporaries. Measure actual peak disk including sidecars and staging.
- Keep existing multi-input/matrix behavior unchanged. Defer new multi-run inference and experiment-wide FDR claims until explicit cross-run rescore populations, competition units and retained fields are defined.

## Dependencies

- **Start / local merge:** WP-01 only. Neither WP-02 hardening nor WP-07 completion blocks existing-engine adapters.
- Final family integration consumes WP-06/07 records through WP-09, not as a prerequisite for local writer/rescorer closure.

## Local merge criteria

- [ ] A pinned existing-engine comparison records `FEATURE_TRANSITION` counts, rescoring time, bytes and format-choice rationale at matched scientific scope.
- [ ] Adapter on/off and results/rescore runs yield identical identities and tolerance-equal scores, confidence and quantities on the existing-engine fixture.
- [ ] Reopened rescore artifacts contain exact required target/decoy/alternative IDs, including rejected candidates, and reproduce the declared rescoring operation.
- [ ] Single-run disk-backed confidence matches an in-memory oracle across input partitions/order, without substituting identification confidence for localization.
- [ ] Injected low-space, write, rename and interruption failures return failure, preserve prior final output and clean owned temporaries without touching inputs.
- [ ] A scale run reports actual peak disk plus final bytes, and existing multi-input/matrix regression outputs remain unchanged.

**Integration gates:** WP-09 owns final family/output workflow checks, and WP-10 owns representative quality/resource/performance acceptance.
