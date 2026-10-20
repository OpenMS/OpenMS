# OpenDIA: laptop-first implementation plan

**Status:** proposed, updated 20 October 2026. [Epic #10661](https://github.com/OpenMS/OpenMS/issues/10661). Ownership is tracked in GitHub.

Build a CPU-first DIA engine for large predicted libraries on a 16 GiB laptop. Target 4 GiB process memory and a 2 GiB low-memory mode, with compact output and preserved scientific quality. These are acceptance targets, not demonstrated capabilities.

## Approach and v1 scope

Replace data handling **one component at a time in the existing pipeline**, behind an option, checking identical or tolerance-equal results. Build new algorithms only for candidate retrieval and family scoring. Mobility and peptidoform families belong in contracts from day one, but do not block the first runnable increment.

- Tool-private code in `src/topp/OpenDIA/`, following DIAuditor/OpenSwathInfer. Move the current entry point first, then extract orchestration into a private library shared by tool and tests. Promote code to core only for a real second consumer. Shared core bugs are fixed in core directly.
- One memory byte budget, bounded queues and fixed stage order. Harden admission in parallel rather than blocking data handling.
- v1 output: compact `results` and optional `rescore`. Keep the statistical population until confidence is finalized. Reduce `FEATURE_TRANSITION` fan-out and measure rescoring early.
- Disk safety: check free space, write temporary output then rename on success, fail explicitly and clean up owned intermediates. Measure final and peak disk use. Defer quota frameworks, diagnostic profiles and elaborate capability manifests.
- **New family-search science is single-run in v1.** Preserve existing multi-input and matrix-export behavior in the compatibility path. Cross-run alignment, experiment-wide FDR and merging new per-run artifacts are deferred, not implied by “global confidence”. WP-08 must derive future retained evidence from those operations before claiming support.

## Five tracks, existing WP IDs

Numbers identify packages, not a serial implementation order. Local merge checks live in each WP. Whole-workflow gates live in WP-09/10.

| Track | Work packages | Start / integration dependencies |
|---|---|---|
| Baseline | [WP-00](work-packages/00-acceptance-baseline.md) | Start now. Pin develop + #10642, ProteoBench Astral and one diaPASEF dataset. |
| Data handling | [WP-01](work-packages/01-contracts.md), [WP-02](work-packages/02-resource-budget.md), [WP-03](work-packages/03-raw-tiles.md), [WP-05](work-packages/05-shared-extraction.md) | WP-01 starts independently of baseline. Others start on v0 contracts, without waiting for hardened admission. WP-05 owns the existing-calibration adapter. |
| Library preparation | [WP-04](work-packages/04-hypothesis-store.md) | Start on WP-01 contracts. Integrate streamed preparation into the current engine. |
| Science | [WP-06](work-packages/06-candidate-retrieval.md), [WP-07](work-packages/07-family-scoring.md) | Start WP-07 null/decoy experiment now. Retrieval integrates WP-03/04/05. Final family scoring integrates extraction/retrieval. |
| Output, integration and benchmarks | [WP-08](work-packages/08-storage-confidence.md), [WP-09](work-packages/09-cli-integration.md), [WP-10](work-packages/10-performance-release.md) | WP-08 and WP-09 start on v0 contracts, not after science. WP-09 completes integration, then WP-10 assesses release. |

See [issue map](issue-posting.md) for live tickets and native dependency relationships.

## Sequence

1. **Baseline and early experiments:** profile the representative current pipeline, establish tool-private v0 contracts, and test the localization null/decoy model against exhaustive and IPF references. Fix shared extractor edge bugs separately on develop.
2. **Incremental wins:** plug raw tiles, streamed preparation, shared extraction and compact output into the existing engine. Prioritize by measured time/RSS, including rescoring/inference/export. Each component can merge with its own tests and compatibility comparison.
3. **New science:** review the localization experiment before freezing family/score contracts. Integrate high-recall retrieval and joint family scoring with mobility, ambiguity and calibrated confidence.
4. **Release decision:** exercise the public raw-to-export path on representative data, then decide whether to change defaults. A runnable early increment is not acceptance of the final scientific slice.

## Shared acceptance rules

- WP-00 pins provenance, numerical tolerances, scientific margins and performance/size targets before comparisons. Local scratch experiments are historical evidence only.
- Every implementation PR builds locally out of tree and runs relevant tests. Packages close on their local merge checks, without waiting for all end-to-end gates.
- Scheduling/storage changes preserve results. Science changes require measured retrieval recall, separate identification/localization error, ambiguity and quantitative quality. Never truncate hypotheses to meet a resource target.
- WP-09 owns integration, real artifact reopening/rescoring and failure paths. WP-10 owns representative cold/warm timing, whole-process memory, peak disk/output size and release compatibility. No default switch before those gates pass.

## Decisions still requiring review

- WP-07 with WP-01/08: family partition, localization nulls/decoys and confidence model, before freezing contracts.
- WP-00/03/08/10: initial acquisition coverage, storage format and numerical quality/runtime/size gates, based on representative measurements.
- Multi-run extension: experiment-level inference and its retained population must be designed before extending the new single-run path.

[Shared architecture and invariants](architecture.md) · [Historical measurements](evidence.md)
