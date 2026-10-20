# WP-05: Exact-query shared extraction and immutable local evidence

- **Suggested issue title:** [OpenDIA][WP-05] Share exact extraction without trace-output explosion
- **Parent epic:** not created
- **Issue:** not created
- **Repository source:** `doc/design/opendia-laptop-first/work-packages/05-shared-extraction.md`

Status: proposed/not implemented
Owner: unassigned
Reviewer: unassigned
Suggested expertise: numerical extraction kernels, ion mobility, lifetime-safe C++ and regression tests.
Repository context: `doc/design/opendia-laptop-first/README.md` and `doc/design/opendia-laptop-first/architecture.md`.

## Goal

Extract equivalent evidence once per complete query while retaining local RT/mobility structure.
Share immutable evidence across families without treating repeated references as independent signals.
Establish mathematical extraction correctness separately from legacy-output compatibility.

## Dependencies

- Requires WP-02 budgets, WP-03 raw tiles
  and WP-04 hypothesis store, using their WP-01 contract dependencies.
- Supplies WP-06 retrieval and WP-07 scoring.
- Integrate public acceptance through WP-09 and WP-00.

## Scope and non-goals

- Implement byte-bounded interval batches, shared extraction kernels and immutable EvidenceBlocks.
- TraceBlock may be a chromatographic view of an EvidenceBlock, not the canonical evidence contract.
- Retain compact RT-by-mobility patches or sparse support, not dense forms-by-fragments arrays.
- Support shared fragments and lazy site-discriminating/loss-channel requests with the same API.
- Do not implement candidate pruning, localization inference or per-tile q-values.
- Do not silently repair legacy extraction while presenting the change as scheduling-only parity.

## Deliverables and interface ownership

- Own query canonicalization, exact equality checks, extraction cache and evidence-buffer lifetimes.
- Query identity includes raw run/content version, lane/FAIMS stratum and acquisition geometry.
- Include MS level, exact RT/scan interval and axis, m/z bounds, tolerance units and edge policy.
- Include isotope/fragment channel, mobility bounds/units, aggregation and missingness policy.
- Include calibration, extraction-policy and numeric-representation versions in cache identity.
- Similar masses, overlapping windows or shared precursor identity alone never authorize sharing.
- Hashes accelerate lookup but full canonical-key equality must resolve collisions.
- Freeze calibration within a pass and invalidate all affected queries after revision.
- Publish contiguous values, ragged offsets, exact axes, validity masks and optional mass errors.
- Keep separated mobility components and raw-region references usable by later scoring.
- Consumers receive leased immutable views, with bounded re-extraction if evidence is evicted.
- Batch by raw locality and expected bytes, including non-monotone interval starts/ends.
- Apply actual PASEF frame/scan geometry and keep discrete FAIMS strata separate.
- Preserve missing acquisition versus measured zero, with independent MS1/MS2 timing.
- Reconcile tile cores/halos without duplicate points or truncated long-interval hypotheses.
- In compatibility/reference mode, preserve complete legacy picker input or equivalent state
  for whole-XIC totals and overlap suppression. Greenfield candidate semantics belong to WP-06.
- Supply counters for unique queries, references, centroid visits, bytes and cache invalidations.

## Normative extraction boundary

- Version the sum/max rule, ppm/Da width convention, isotope spacing and endpoint inclusion.
- Proposed reference compatibility sum policy uses strict open m/z bounds: `lower < mz < upper`.
- This proposal requires WP-01 approval and is not a fixed new greenfield extraction policy.
- Explicitly specify RT/IM boundary inclusion separately in WP-01 before implementing the kernel.
- Hand-computed unit cases must cover first/last centroid, singleton/empty scans and duplicate m/z.
- Test exactly-on-boundary and just-inside/outside values in both non-IM and IM paths.
- Add a first-centroid regression where a multi-step left walk must include the first valid point.
- The earlier nonzero fixture suggests a legacy omission, not an approved normative behavior.
- Isolate its regression and reviewed fix/policy decision before accepting migration comparisons.
- Keep in-process numeric tolerances separate from lossy exported chromatogram encoding tolerances.

## Code placement

- Provisional until implementation: reusable extraction kernels belong in core `PROCESSING`/`IONMOBILITY`, with matching header/source paths.
- OpenDIA query orchestration, canonicalization policy and evidence-cache policy belong in NEW proposed `src/openms/include/OpenMS/ANALYSIS/OPENDIA/` and `src/openms/source/ANALYSIS/OPENDIA/`, neither currently present. This plan creates no directories or code.
- Core must not depend on OPENDIA. Existing shared `OPENSWATH` components stay shared, and `src/topp/OpenDIA.cpp` remains a thin CLI.

## Output-size contract

- Own evidence accounting by unique complete query, stored support points, shared axes, validity masks and consumer-reference edges. Encoded payload grows with unique evidence plus sparse references, not consumers multiplied by trace length or a dense family/fragment/RT/IM cube.
- Persist each shared evidence block once per representation and reference it by stable ID; count dictionary, axis and sparse-incidence bytes separately. Charge recomputed/overlapping representations explicitly rather than concealing duplicate storage.
- Extracted evidence is working inference scratch or an explicitly reusable cache, not a default deliverable. Budget scratch, final deliverables and caches separately plus per-volume combined peaks, counting WAL/journals, indices, atomic staging and compaction old+new copies.
- Default results-only output contains compact identification/localization/quantification summaries and ambiguity provenance, never whole-library traces/transition cubes. Full-OSW compatibility is an opt-in export only.
- Optional rescore profile preserves declared prefilter target/decoy/alternative inputs required by supported rescoring, including necessary evidence or durable references. Diagnostics are opt-in by explicit query/region/ID selection and byte capped before serialization.
- Keep full inference population and required evidence available until global confidence completes, including legacy whole-XIC state where claimed. Cache eviction may re-extract, not silently drop candidates or prune evidence to satisfy quota.
- Fail explicitly when mandatory output cannot fit, without complete-looking truncated results. Compression or Parquet alone is insufficient, and CLI budget names remain proposals only.

## Small PR slices

1. Specify canonical queries and normative boundary fixtures, including the first-centroid regression.
2. Implement bounded scalar reference extraction and immutable evidence leases for non-IM and IM.
3. Add batched shared-query kernels, sparse mobility support and calibration-aware invalidation.
4. Add stitching/picker handoff and public end-to-end comparisons with a reviewed compatibility policy.

## Acceptance checks

- [ ] Review header/source placement and build dependencies: reusable core has no OPENDIA dependency, shared OPENSWATH stays shared, and the CLI contains only thin orchestration.

- [ ] Multiply consumers of fixed exact queries and verify payload/support-point counts stay fixed while only reference edges grow, including serialized bytes.
- [ ] Verify results-only writes no whole-library trace payload, selected diagnostics stay byte capped, and declared rescore inputs remain complete.
- [ ] Exhaust evidence/export quota with pinned consumers and staging/compaction overlap; verify explicit failure rather than changed inference populations or truncated successful output.

- [ ] Unit: normative sums/maxima and tolerance conversions pass independently of legacy output.
- [ ] Unit: vary every query-key field and force hash collisions to prove no false evidence reuse.
- [ ] Unit: shared queries return identical evidence and release buffers only after the last consumer.
- [ ] Unit: non-monotone intervals, overlapping PASEF geometry, missing scans and FAIMS remain correct.
- [ ] Unit: cancellation, oversized patches and stalled consumers obey WP-02 without partial success.
- [ ] Integration: run the real public OpenDIA pipeline from raw/prepared inputs through final export.
- Use a coherent build once WP-09 is integrated, with nonzero MS1 and MS2 positive fixtures.
- Include real IM overlap, known-site mixtures, unresolved forms and weak discriminating fragments.
- Compare sharing on/off, tile sizes and 2/4 workers for evidence, candidates, scores and quantities.
- Validate boundary cases against declared mathematical policy, not blind buggy legacy equality.
- Require fresh producer success before output checks, including ordinary DIA and PASEF fixtures.
- Measure unique-query reduction and full pipeline time/RSS/disk, not an extraction-only speed claim.
- Proposed 2/4 GiB settings are acceptance targets, not established boundedness or speedup results.
- Save commands, policy/build IDs, tolerances and fresh output checks in the WP-00 ledger.

## Failure risks

- An incomplete query key shares incompatible evidence and can bias localization or interference scores.
- In compatibility/reference mode, halos alone do not reproduce legacy whole-input picker decisions.
- Greenfield boundary/candidate behavior must satisfy the WP-06 contract, not implicit legacy parity.
- Cache pinning and dense mobility patches can exceed budgets unless charged before production.
- Shared evidence must not be counted repeatedly as independent support by downstream inference.

## Handoff

Give WP-06/WP-07 canonical query examples, normative fixtures and lifetime/ownership guarantees.
Provide the compatibility decision and measured sharing counters to WP-09/WP-10.
Keep scientific and resource acceptance pending until the integrated public pipeline passes.
