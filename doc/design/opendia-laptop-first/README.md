# OpenDIA: laptop-first implementation plan

**Status:** proposed, updated 24 October 2026. [Epic #10661](https://github.com/OpenMS/OpenMS/issues/10661). Ownership is tracked in GitHub.

## Current verification status

Full fresh workstation A/B output comparisons passed at selected clean merged source
`3e47bc7e7468656d0d53237f11fc4f049d278efe`, including #10698 (merged
23 October 2026), with independently pinned harness `5a804f5e7c92f5c7e541050123a66f2978d23368`.
TimSim target-only, mixed and all-three Astral comparisons passed all 12 files
with zero changed numeric cells or scientific identities. Raw rows per repeat
are 11,484 / 6,783 / 6,274 and standard accepted counts 6,000 / 5,919 / 4,343.
Independent provenance review, comprehensive resource/stage acceptance, and owner-approved candidate gates remain open. The separate profiling branch builds and passes its profile fixture plus all 29 OpenDIA TOPP checks. The branch is for refactoring work and is not planned for integration. Four instrumented cold/warm followups passed frozen-output equivalence and produced stage summaries. Standalone TimSim rescoring processed 16,441 rows and wrote 13,421,129 apparent bytes; Astral rescoring was refused because the frozen SQLite source retains WAL/SHM sidecars. A local evidence archive passed clean-extraction verification, but durable public publication remains open. These separate-binary results are not same-binary repeat or candidate-acceptance evidence.

PXD004573's 579 published strings are authenticated, but raw retrieval, site
mapping, purity and blank interpretation remain unresolved. See WP-00 for the
requirement-specific map. Earlier failed repeats remain historical evidence; the
22 October same-binary repeat checks passed on a premerge build, not on this
selected source. The corrected Parquet scores are not interchangeable with
affected historical scores.

The declared workstation protocol uses 17 OpenDIA threads, with independently
captured engine/harness provenance and actual subprocess/affinity settings.
Historical four-thread laptop timings are not a matched baseline. The
[resource measurement contract](https://github.com/singjc/OpenMS/blob/9abc264f4c95cb550d67490030dcac649750060c/tools/scripts/opendia_benchmark/README.md#dedicated-workload-accounting-and-limitations)
distinguishes verified final dedicated-scope kernel peaks, diagnostic samples,
inner/outer timing and sampled disk lower bounds. Four successful instrumented
followups additionally provide stage timing and aligned sampled-resource summaries,
but these are descriptive observations, not exact peaks or acceptance evidence.

[WP-00](work-packages/00-acceptance-baseline.md) records the frozen existing
repeat-comparator policy, the exact pending evidence and the owner/reviewer
decision matrix. Proposed candidate no-regression/runtime/footprint gates are
not approved thresholds. Scientific accuracy margins and known-site localization
truth remain open. Freeze reviewed gates before candidate evaluation, complete
adversarial review before any PR, and keep public GitHub updates pending evidence.

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
| Baseline | [WP-00](work-packages/00-acceptance-baseline.md) | Start now. Pin develop + #10642 and the established TimSim diaPASEF fixture and exploratory Astral inputs below. |
| Data handling | [WP-01](work-packages/01-contracts.md), [WP-02](work-packages/02-resource-budget.md), [WP-03](work-packages/03-raw-tiles.md), [WP-05](work-packages/05-shared-extraction.md) | WP-01 starts independently of baseline. Others start on v0 contracts, without waiting for hardened admission. WP-05 owns the existing-calibration adapter. |
| Library preparation | [WP-04](work-packages/04-hypothesis-store.md) | Start on WP-01 contracts. Integrate streamed preparation into the current engine. |
| Science | [WP-06](work-packages/06-candidate-retrieval.md), [WP-07](work-packages/07-family-scoring.md) | Start WP-07 null/decoy experiment now. Retrieval integrates WP-03/04/05. Final family scoring integrates extraction/retrieval. |
| Output, integration and benchmarks | [WP-08](work-packages/08-storage-confidence.md), [WP-09](work-packages/09-cli-integration.md), [WP-10](work-packages/10-performance-release.md) | WP-08 and WP-09 start on v0 contracts, not after science. WP-09 completes integration, then WP-10 assesses release. |

See [issue map](issue-posting.md) for live tickets and native dependency relationships.

## Established benchmark inputs and tooling

- **TimSim `validation_3x3`:** six diaPASEF inputs, three controls and three treatments, with target-only and target-plus-entrapment transition-library searches and packaged identification/abundance truth.
- **Astral:** three HeLa 50 ng, 15 min technical replicates. The current FASTA is a randomly sampled 100k-peptide database. This is exploratory, not official full ProteoBench global precursor validation. Current exports lack global precursor q-values and the exact modified-peptide pair mapping needed for paired FDP. Global peptide/protein scores do not fill those gaps.
- [Benchmark scripts](https://github.com/singjc/OpenMS/blob/9abc264f4c95cb550d67490030dcac649750060c/tools/scripts/opendia_benchmark/README.md) establish build, runners, metrics and report generation. System-resource monitoring records machine, RAM, logical threads and available RAM/idle CPU during each search. Fresh full TimSim and three-file Astral runs completed successfully on 21 October 2026 at source `d3b00d382efc`, with a regenerated report/PDF. See [WP-00 #10662](https://github.com/OpenMS/OpenMS/issues/10662) for results. This is a measured current-pipeline baseline, not completion of acceptance gates.

Neither dataset supplies known-site localization truth. WP-00 must select separate known-site data before localization acceptance.

Reproducibility packaging and strict repeat-comparison tools are now established.
Historically, the 21 October independent TimSim target repeat preserved all 6,000 accepted precursor
identities but failed frozen score/quantity and raw-identity checks. The mixed
TimSim and three-file Astral repeats also completed successfully as searches,
but failed strict scientific comparisons. Their standard-filter accepted identity
totals changed from 5,949 to 5,935 and 4,337 to 4,291 respectively. Final evidence
and provenance gaps are tracked in [WP-00](work-packages/00-acceptance-baseline.md).
A verified bundle is not proof that the scientific repeat gate passes.

## Sequence

1. **Baseline and early experiments:** profile the representative current pipeline, establish tool-private v0 contracts, and test the localization null/decoy model against exhaustive and IPF references. Shared extractor edge bugs were fixed separately on develop in #10673; retain their mathematical edge checks.
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
