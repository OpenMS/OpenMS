# Shared evidence and scoring lifetime, 4 November 2026

[Master plan](README.md) · [WP-02](work-packages/02-resource-budget.md) · [WP-05](work-packages/05-shared-extraction.md) · [WP-09](work-packages/09-cli-integration.md) · [WP-10](work-packages/10-performance-release.md)

## Delivered engineering increment

The current implementation is stacked on `feature/opendia-wp05-readonly-extraction`, not merged to `develop`.

- `fce959ce726d68b03bfca380ec546da7022277dc` shares calibrated immutable scan evidence between extraction and spectral/IM scoring. Mutable public access remains detached. Scoring copies fell from 9,432 to zero, but the two wall-time pairs were mixed. That increment alone did not demonstrate a useful speedup.
- `3c13c1593065c2c23c95879914fa1892cab53df2` separates scorer preparation, one-candidate scoring and finalization, reusing the same numerical body. The private active SQLite OSW path scores actual picked candidate apices in stable RT order within the existing assay batch. Complete chromatograms, candidates, neighbor spectra, scoring formulas, ranking/ties and original assay output order remain intact. Core, MRM, featureXML and Parquet scheduling remain unchanged.

### Measured native fixture outcome

Two counterbalanced full-workflow pairs, identical raw/PQP/options, eight threads, 256 MiB raw cache, 32 MiB maximum scan and fixed seeds:

| Pair order | Parent wall | Actual-apex wall | Parent peak RSS | Actual-apex peak RSS |
|---|---:|---:|---:|---:|
| Parent then candidate | 116.13 s | 40.95 s | 1,429,324 KiB | 1,796,440 KiB |
| Candidate then parent | 116.66 s | 39.13 s | 1,524,652 KiB | 1,663,348 KiB |
| Median | 116.395 s | 40.040 s | 1,476,988 KiB | 1,729,894 KiB |

This is a **2.91x native-fixture speedup**, with both pairs improving and the prospective speed gate passing. Median RSS is **17.12% higher**. Scoring decodes fall 6,480/6,439 to 860/875, about 86.6%. Total decodes fall 11,231/11,242 to 5,669/5,686, about 49.5%. Total calibration copies remain 460 in every run. This is reduced rereading, not once-only access.

Scientific acceptance over all four runs passes: export 4,699 rows and 55 scientific columns, complete OSW evidence/score tables within established tolerances, exact precursor/peptide/protein matrices with 942/942/251 rows, and identical rank-one target precursor sets at 1% (944) and 5% (1,017) FDR. Nonzero IM is present. Only generated run/file/feature identities are normalized. Separate calibration-coefficient auditing is not claimed.

All **176 existing selected OpenDIA/OpenSwathWorkflow/MRMFeatureFinder/OpenSwathAnalyzer CTests** pass after rebuilding their prerequisites. Actual native five-lane integration and pressure/lease checks also pass, as do one/eight-thread scorer checks. Final committed executable and three shared libraries byte-match the measured artifacts. Installed-consumer replay and clang-format availability are not claimed by this increment.

The frozen raw is native TimSim `OpenSwathTimSim_control_01.d.zip`, SHA256 `4d7ab11d27132d5295a4bb12083bec40b9db67b6dcded0a4011c81f62afcaa10`. The PQP SHA256 is `3b1ac20cb2ecba9fec0d445c2c96bcbd74d1e463dfb78ffea0fb8a5c2b20a03e`, with 4,000 target/decoy precursor entries and 24,000 transitions. Local detailed commands, comparison results and provenance are retained in `apex-scoring-20261104`; this is not a published portable benchmark bundle.

## Scoring-context memory follow-through

The parent actual-apex implementation prepared every group's S/N estimates, detection/identification subsets and scoring scratch before the first task, then retained every context until all tasks finished. Existing assay-batch sizing limits count/lifetime, not total process bytes.

Commit `b539859c3bac006ee8046be85ce7da855621aad8` implements on-demand context preparation and prompt release at the last consumer, preserving global candidate RT order and original output order. Non-IM deferred REMAINDER scoring retains its original raw-access/finalization order. It does not truncate hypotheses or complete traces, add a cache tier, enlarge caps or introduce a generic resource framework. Empty groups are validated independently; failed groups release their scratch without suppressing healthy groups. Output buffers preserve input-group ordering even when execution order differs.

Before candidate measurements, the declared acceptance was: all scientific outputs unchanged, at least 5% lower median process peak RSS with lower RSS in each of two counterbalanced pairs, and no more than 5% median wall regression against the retained actual-apex implementation. The gate **passes**, without changing inputs, options, cache limits or seeds:

| Pair order | Apex parent wall | Shorter-lifetime wall | Apex parent peak RSS | Shorter-lifetime peak RSS |
|---|---:|---:|---:|---:|
| Parent then candidate | 36.46 s | 37.04 s | 1,605,184 KiB | 1,552,444 KiB |
| Candidate then parent | 41.30 s | 40.55 s | 1,771,480 KiB | 1,458,440 KiB |
| Median | 38.880 s | 38.795 s | 1,688,332 KiB | 1,505,442 KiB |

Median process peak RSS is **10.83% lower**, with both pairs lower. Wall time is effectively unchanged (0.22% lower median), not a new speedup claim. The previous 2.91x comparison used a different parent and measurement session; do not multiply improvements or infer a combined earlier-parent result.

Scoring decodes are 893/899 in the parent and 873/910 in the candidate. Total decodes are 5,693/5,726 versus 5,687/5,774. Scoring copies remain zero and total calibration copies remain 460 in every run. These counters do not establish once-only decoding or a new decode reduction from this lifetime change.

### Requirement-to-check results for the retained memory change

| Requirement / changed output | Whole-result check and observed result |
|---|---|
| Prepare only demanded scoring scratch and release at last use | Exact production scorer regression counts peak live contexts: disjoint groups need one, overlapping groups two, versus all three including an empty group in the parent. Deferred raw work retains contexts until its original finalization phase. |
| Preserve complete traces, candidates, ranking/ties and output order | Ordinary versus batch scorer checks exact features, chromatograms and scores for overlapping/disjoint groups, reverse output/execution order, aliased output maps, detection subsets, UIS, MS1, IM-enabled scoring and report limits -1/0/1/2/10. All pass. |
| Preserve request ordering and independent errors | Exact candidate/deferred REMAINDER request sequences pass. Empty and failed groups, invalid nonfinite groups and whole-batch nonfinite fallback pass. |
| Preserve scientific public outputs | All four real CLI runs pass: 4,699 exported rows and 55 scientific columns exact after generated run/file/feature identity normalization, all eight OSW feature/score tables within relative 1e-6 / absolute 1e-8 tolerance, exact 942/942/251-row precursor/peptide/protein matrices, identical target ID sets at 1% (944) and 5% (1,017), and nonzero IM. No fallback. |
| Preserve neighboring workflows and native ownership | Local Release build and all 178 affected CTests pass, including core scorer and native tests. One/eight-thread scorer checks pass. Actual native integration passes five lanes, four acquisitions/four decodes, 14 raw hits, lease/pressure and physical-frame reuse checks. |
| Demonstrate useful process-memory reduction without losing speed | The two prospective counterbalanced pairs above pass the RSS and wall gate. No extra pairs or parameter search. |
| Validate the actual retained result, not an earlier patch | Final committed source matches the frozen complete patch; executable and three libraries byte-match measured artifacts; exact raw/PQP hashes and all commands verified before and after runs. Implementation working tree is clean. |

Detailed commands, acceptance, comparisons, source/artifact audit and logs remain in the local `scoring-memory-20261104` evidence bundle. They are not a published portable benchmark package. Overlapping candidate spans or deferred raw work can still retain the whole existing assay batch. Complete picked groups and output buffers remain materialized. This is shorter lifetime, **not a strict byte cap**, a guarantee of extract-once processing, or proof of large predicted-library scaling.

## Roadmap and PR boundary

- The earlier R1 indexed sparse extraction experiment was reverted after its integrated path was slower. It is **not** a completed retained redesign milestone. R2 family scoring, R3 fused native execution and R4 scale/scientific acceptance remain open.
- This increment reuses the current scorer as an engineering bridge. It neither freezes the future scientific model nor substitutes a full legacy score suite plus another model for the proposed sparse engine.
- WP-02 owns pre-materialization whole-job byte admission. Shorter context lifetime alone cannot close that gate.
- WP-03/05 retain native evidence geometry, calibration and ownership obligations. Extraction, calibration, cross-lane and cache-eviction rereads still exist.
- WP-09 owns integration and existing-path compatibility. WP-10 owns Astral/native Thermo, large predicted-library, cold/warm, complete process/child memory, disk and competitive-quality acceptance. One TimSim fixture cannot close those gates.
- Keep implementation and planning branches separate. This document is not authorization to open a PR, merge into `develop`, close work packages, or claim release readiness. Notify the owner when a scoped implementation PR is ready; the combined engine remains subject to the existing scientific/resource acceptance decision.
