# OpenDIA performance exploration and engine ownership

**Owner direction, 27 October 2026.** Optimize for competitive end-to-end identification and quantification, not preservation of the current accessor/scorer architecture. This supersedes the incremental-only next-step restrictions in the earlier WP-03/05 notes. This is an experiment and implementation plan, not evidence that any proposed engine or model is faster.

## What changes, and what does not

- New raw access, indexing, scheduling, extraction, candidate generation and scoring implementations belong in `src/topp/OpenDIA/include`, `source` and `test`. Entire chunks may be replaced. Existing OpenMS/OpenSWATH accessors, spectrum containers and scorer entry points are reference/compatibility adapters, not mandatory internal interfaces.
- Reuse correct low-level vendor decoding/calibration where it helps. A rewrite does not require reimplementing a compression codec or adding a dependency. Keep shared OpenSWATH users working and fix genuine generic defects in core. Do not place OpenDIA-only architecture in the shared library.
- Bruker and Thermo support stay enabled. Read original vendor data directly. Do not require or silently convert it to indexed mzML. A small optional derived index is an experiment only if its build, storage, invalidation and cold-run costs are included. Do not hide a full converted copy behind the word index.
- Preserve exact acquisition coordinates, calibration provenance, discrete FAIMS lanes, missingness, complete statistical populations and byte ownership. These are correctness requirements, not reasons to retain an inefficient algorithm.
- Behavior-preserving engineering must retain reference identities/evidence and tolerance-equal outputs. A new scientific scorer is a different evaluation lane: it may change scores and accepted sets, but must meet independent error-control and quantitative non-inferiority gates. Exact agreement with old OpenSWATH scores is not a requirement for a scientifically different model.

## Current evidence is not an optimality claim

The tested private adapter reduces memory on one historical comparison at substantial runtime cost. The latest known-size admission observation reduces scoring decodes and decoded bytes, but takes 312.19 s versus the previous 291.96 s; scoring mutable copies remain 9,432. These are groundwork results, not a competitive fast engine.

The scoring phase owns most observed raw decodes in the earlier attribution run. That does **not** prove the scoring formulas or final classifier consume most exclusive CPU/wall time. Read timing includes waiting, nested decode and copy timing is not a disjoint critical-path decomposition, and worker times overlap. Replacing the classifier cannot be assumed to eliminate vendor decoding, evidence extraction or feature computation.

## Immediate experiment: native regions through real scoring

The first concrete storage/query slice is implemented and checked in [native column-region contract](native-region-contract.md). Its native producer is currently a compatibility bridge, not a frame-native backend or an integrated scoring engine. A separate [calibration getter correction](work-packages/05-shared-extraction.md#native-region-contract-and-calibration-hotspot) shows observed full native wall 309.34 to 231.23 s, user CPU 1,124.07 to 235.10 s, with all twelve outputs preserved. RSS rises; the comparison is instrumented and unbalanced. Keep that distinction when reporting measurements.

Do this bounded experiment before another general accessor migration or model-training campaign. WP-03 owns access, WP-05 extraction and the fused scoring handoff, WP-07 scorer alternatives, and WP-10 the comparison decision.

1. **CPU diagnostic completed; remaining attribution is scoped.** The actual unchanged full-library native Bruker profile identified 75.61% self user CPU in shared calibration getters, reduced to 2.57% by `c64e647997`. Native sorting is now the largest self symbol at 25.62%. This is not an allocation or off-CPU lock profile, nor an exclusive stage wall decomposition. Do not repeat the baseline campaign. Use actual calibrated requests and complete wider-IM evidence when integrating the region consumer, and measure cold source/index/workspace costs there. Thermo's direct contract checks pass, but full Thermo science/performance is not inferred from Bruker.
2. **Three bounded access/execution candidates, the same request set:**

   | Candidate | Concrete implementation hypothesis | Decisive observation |
   | --- | --- | --- |
   | A: current adapter | Existing source/index, scan LRU, detached scoring copies | Control, including observed decodes, copies, actual outputs and full runtime |
   | B: native-region access | OpenDIA-owned frame/block and acquisition index, bounded decoded blocks, leased column slices; batch overlapping RT/isolation/IM requests before m/z lookup | Fewer redundant frame decodes and allocations with equal exact-query evidence and lower cold-inclusive cost |
   | C: acquisition-local execution | Decode a bounded region once, fuse extraction and batched numerical scoring for all relevant hypotheses, then retain sparse evidence for revisits | Lower full workflow cost without truncated candidate windows, mobility evidence or statistical populations |

   Reuse existing Bruker/OpenTIMS and Thermo bridge capabilities first. Investigate lower-level frame/batch access rather than forcing a round trip through `MSSpectrum`, `SpectrumPtr` or `ISpectrumAccess`. Source-specific decoders may differ behind one OpenDIA region contract. A source capability limit must be measured, not worked around by mandatory format conversion. Do not install a new backend or write all three engines before testing the smallest native-region-to-real-score slice.
3. **Measure useful random access, not just call latency.** Preserve the same physical query support and score outputs. Compare real-order requests with an equivalent acquisition-grouped order, reconciling output identity and reduction order. Record region/query throughput and latency distribution, decoder calls and materialized bytes, allocation/copy work, lock contention, total wall time, process/cgroup peak memory, temporary/final disk and cold index cost. Raw file reads and decoded materialized bytes are distinct measures. A memory-mapped file or a warm OS cache alone is not evidence of fewer decompressions.
4. **Build one real vertical slice.** Put the winning prototype in OpenDIA, driving the real extraction/scoring consumer without mutable legacy facades. The focused native regression must cover calibration, exact and wider-IM support, interval/tile boundaries, empty/missing data, clone/lifetime/cancellation and final release. Build affected real targets. Run one unchanged full native Bruker comparison over all twelve outputs, then validate the chosen path on native Thermo with its stated acquisition scope. Tiny reader checks do not establish full Thermo scientific acceptance.
5. **Stop or promote from evidence.** Reject a candidate if preparation, model or retained-evidence cost erases its end-to-end benefit. If only repeated use benefits, report index break-even runs and cold/warm results separately. Do not hide a memory/disk increase. Promote one candidate to representative repeated quality/resource comparisons, not an unlimited matrix or another publication/readback loop.

This is not permission to weaken unsafe-admission/cancellation guards, enlarge caps silently, truncate search scope or mark an unrun candidate accepted. The automatic low-cap failure remains a separate real regression that the selected execution strategy must address.

## Scoring exploration is allowed, including replacement

Scoring is not frozen to the existing handcrafted feature suite. Evaluate two levels separately:

- **Engineering control:** fuse and batch existing DIA/IM integrations, reuse calibrated evidence/isotope computations, use contiguous sparse data and vectorizable kernels where the profile justifies them. This isolates dataflow benefits from a new statistical model and supplies an explainable reference.
- **Learned alternative:** first compare the current confidence scorer against one compact CPU model on the same retained candidates and evidence. A small MLP is a hypothesis, not a selection in advance. If handcrafted feature generation is the real cost, test one evidence-native model that can actually replace those expensive features, rather than running the whole old pipeline and adding a network afterwards. Large transformers, GPU dependence and training from scratch are not prerequisites.

Charge and time feature/evidence preparation, fitting, calibration, inference and retained model/data storage. Keep training, tuning and final evaluation disjoint at peptide/family and run/instrument levels appropriate to the claim. Same-spectrum/shared-evidence leakage and near-identical peptidoforms require explicit grouping. Freeze model/version/seeds, target-decoy construction and split identities before the held-out comparison. Check both identification and localization calibration, including unresolved/coexisting forms. Train on targets/decoys only through a justified selection procedure; target labels are not ground truth. Never use held-out truth or an entrapment panel to tune the final model and then report that panel as independent validation.

Do not choose a neural network because another tool uses one. Evaluate accuracy, calibration, CPU cost, memory, training data/provenance and redistribution licenses together. A model can improve identification while leaving quantification or runtime worse.

## Comparable to DIA-NN and Spectronaut means measured science

WP-10 must compare pinned available versions/configurations of OpenDIA, DIA-NN and Spectronaut on the same declared native-input/hardware/search scopes. Tools or licensed runtime unavailable locally are explicit unrun comparator cells, not a reason to stop the OpenDIA prototype or invent superiority. Do not copy proprietary implementation or redistribute restricted model weights.

- Separate empirical-library, predicted/library-free, single-run and match-between-runs workflows. Match FASTA/PTMs, candidate scope, library/prediction policy, inference units and run/global confidence filters as far as supported; publish unavoidable differences.
- Compare precursor/peptide/protein identifications at matched externally checked error levels, not raw ID totals or a common nominal q-value alone. Map protein grouping consistently. Report entrapment/FDP with the correct search-space ratio, known-truth sensitivity, completeness and identification reproducibility. Current sampled-FASTA Astral diagnostics are not full official ProteoBench validation.
- Compare quantitative ratio bias/compression, precision/CV, missingness and differential-abundance false positives/power on independent mixture/abundance truth. Separate intersection-based precision from coverage; selecting only shared easy IDs must not hide lost sensitivity. Known-site localization remains its own truth gate.
- Include complete cold/warm time, preparation/model/prediction costs, rescoring/export, peak whole-process/child memory and disk. Keep the CPU laptop target primary; GPU results require a separately declared hardware/cost comparison.
- Predeclare numerical non-inferiority margins, repeat policy and primary endpoints before selecting a winner. Existing future-candidate proposals are not approved numerical science margins. Explore the quality/runtime/memory Pareto frontier without retrospectively weakening failed gates. No claim of parity or superiority has been established yet.

## Primary evidence consulted, not OpenDIA benchmark results

- AlphaTims, *Molecular & Cellular Proteomics* (2021), "Indexing Trapped Ion Mobility Spectrometry-TOF Data for Fast and Easy Accession and Visualization", PMCID PMC8526765: a concrete sparse/indexed raw-coordinate access design. Its import/index and memory costs must not be assumed to fit OpenDIA's budget.
- OpenTIMS, *Journal of Proteome Research* (2021), DOI `10.1021/acs.jproteome.0c00962`: native frame access and comparative decoding evidence, not a guarantee for our query workload.
- DIA-NN, *Nature Methods* (2020), DOI `10.1038/s41592-019-0638-x`, and the official DIA-NN README consulted 27 October 2026: neural classification, interference correction and current cross-validation guidance. The original paper is not the specification of current releases.
- DIA-BERT, *Nature Communications* (2025), DOI `10.1038/s41467-025-58866-4`: learned identification/quantification is a legitimate alternative, but large training sets and model cost cannot be ignored.
- AlphaDIA, *Nature Biotechnology*, online 21 October 2025, DOI `10.1038/s41587-025-02791-w`: feature-free signal handling, transfer learning, entrapment and quantitative benchmarking. Its reported ion-mobility runtime limitations also caution against equating newer ML with faster analysis.

These sources motivate bounded experiments, not new dependencies, ready-made compatibility or proof that a published accuracy gain transfers to OpenDIA.
