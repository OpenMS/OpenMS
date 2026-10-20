# WP-04: Bounded hypothesis preparation and indexed family store

- **Suggested issue title:** [OpenDIA][WP-04] Bound hypothesis preparation and family indexes
- **Parent epic:** not created
- **Issue:** not created
- **Repository source:** `doc/design/opendia-laptop-first/work-packages/04-hypothesis-store.md`

Status: proposed/not implemented
Owner: unassigned
Reviewer: unassigned
Suggested expertise: peptide chemistry, digestion/decoys, prediction and external-memory indexes.
Repository context: `doc/design/opendia-laptop-first/README.md` and `doc/design/opendia-laptop-first/architecture.md`.

## Goal

Prepare and query hypotheses without materializing a whole predicted library or localization graph.
Bound the entire path from FASTA or imported library through digestion, decoys and prediction.
Preserve search coverage, stable identities and provenance independently of batching and storage.

## Dependencies

- Requires WP-01 contracts and WP-02 budgets.
- Supplies WP-05 extraction and WP-06 retrieval.
- Coordinate evidence identity consumers with WP-08.
- Public acceptance closes through WP-09 using WP-00.

## Scope and non-goals

- Stream imported PQP/library records and FASTA digestion with explicit enzyme/search settings.
- Bound missed-cleavage expansion, modification enumeration, deduplication and protein mappings.
- Generate decoys with equivalent supported search/preparation rules and recorded random seeds.
- Build searchable disk-backed family/precursor/fragment indexes using bounded sort/merge passes.
- Include model loading, batched predictions, index builds, temporary files and cold-start time.
- Do not implement candidate gating, localization confidence, global FDR or output export here.
- Do not claim lazy prediction if building the retrieval index already requires full predictions.

## Deliverables and interface ownership

- Own preparation iterators, durable hypothesis-store implementation and paged query adapters.
- Use WP-01 stable IDs for backbone/composition families, forms, charges, fragments and decoys.
- Initially partition families by relevant charge/adduct and precursor constraints.
- Preserve site alternatives, loss-channel annotations and target/decoy provenance losslessly.
- Retain sequence/protein/source metadata out of hot loops but available through bounded lookup.
- Expose paged precursor m/z, RT and optional mobility interval-overlap queries.
- Store uncertainty and missing predictions, allowing broad/recovery queries rather than false exclusion.
- Queries return bounded pages and continuation tokens, not an in-memory full-store directory.
- Freeze run calibration separately from the reusable index and version run-specific query mappings.
- Generate shared fragments once and discriminating fragments lazily for requested local hypotheses.
- Preserve a sparse fragment-to-hypothesis incidence relation without pairwise UIS expansion.
- Stream huge families or spill partial work while keeping complete alternative/null context.
- No hidden top-k, beam or modification truncation is allowed as a memory optimization.
- Any explicit approximation records policy, excluded scope and incomplete-search status for evaluation.
- Persist preparation manifests: input hashes, chemistry settings, seeds, model/version and schema.
- Make restart/reuse validate manifests and reject incompatible or partially published indexes.
- Reserve all preparation buffers/workspaces via WP-02 and enforce temporary/final disk quotas.
- Publish import/digestion/decoy/prediction/index timings and peak resident/temporary/final bytes.

## Code placement

- Provisional until implementation: reusable digestion, modification and decoy chemistry belongs in core `CHEMISTRY`, with matching header/source paths.
- OpenDIA family index, preparation orchestration and provenance policy belong in NEW proposed `src/openms/include/OpenMS/ANALYSIS/OPENDIA/` and `src/openms/source/ANALYSIS/OPENDIA/`, neither currently present. This plan creates no directories or code.
- Core must not depend on OPENDIA. Existing shared `OPENSWATH` components stay shared, and `src/topp/OpenDIA.cpp` remains a thin CLI.

## Output-size contract

- Own a preparation ledger for distinct families/forms/precursors/fragments, dictionary entries, protein mappings, sparse incidence edges and index entries, with encoded bytes per record plus variable metadata. Store shared fragments once, never pairwise alternative/UIS expansions or scan-by-library products.
- Prepared stores are reusable caches, not default final deliverables; sort runs and prediction intermediates are working inference scratch. Budget these classes separately plus per-volume combined peaks, including WAL/journals, indices, atomic staging and compaction old+new copies.
- Final results-only output contains compact identification/localization/quantification summaries and ambiguity provenance using IDs/dictionaries and shared coordinate/evidence references, not a copied predicted library or whole-library traces/transition cubes. Full-OSW compatibility is opt-in export only.
- Optional rescore profile preserves declared prefilter target/decoy/alternative inputs required by supported rescoring. Resolve referenced chemistry/provenance through a charged self-contained subset or declared durable dependency, not an assumed surviving temporary store.
- Diagnostics are opt-in, selected and byte capped. Compression or Parquet alone is insufficient: report prepared-store and final-output bytes/cardinalities separately.
- Preserve every supported enumerated hypothesis and the full inference population until global confidence completes. No silent top-k, dropping or evidence pruning for quota; if mandatory storage/output cannot fit, fail explicitly without complete-looking truncated results.
- CLI budget names are proposals only; index construction must reserve peak merge/publish space before work, not infer feasibility from final compressed size.

## Small PR slices

1. Add deterministic family IDs, paged imported-library access and provenance round-trip tests.
2. Stream digestion/modifications/decoys with bounded deduplication and complete coverage checks.
3. Add budgeted prediction batches, lazy fragments and external-memory index construction.
4. Add restart validation and public cold/reuse pipeline tests across growing search spaces.

## Acceptance checks

- [ ] Review header/source placement and build dependencies: reusable core has no OPENDIA dependency, shared OPENSWATH stays shared, and the CLI contains only thin orchestration.

- [ ] Grow absent hypotheses and alternatives independently: report exact store row/edge counts and bytes, with no pairwise UIS or library-by-scan expansion.
- [ ] Verify results-only exports do not copy the prepared library, while mandatory ambiguity links and declared rescore populations/provenance remain recoverable.
- [ ] Exhaust quota during external merge and atomic publication, measuring old/new plus journal/index bytes and rejecting incomplete stores/results as complete.

- [ ] Unit: small exhaustive digestion/modification/decoy truth sets match streamed enumeration exactly.
- [ ] Unit: IDs, protein mappings and target/decoy associations remain stable across page/worker sizes.
- [ ] Unit: huge family pagination retains all alternatives and does not replicate shared evidence.
- [ ] Unit: absent RT/IM predictions, broad uncertainty and calibration changes preserve recovery access.
- [ ] Unit: prediction workspace estimates, sort spill, disk-full and cancellation never escape budgets.
- [ ] Unit: incompatible manifests and interrupted builds cannot be reused as complete stores.
- [ ] Integration: execute a fresh public OpenDIA preparation-to-export workflow once WP-09 is wired.
- Exercise both imported libraries and library-free FASTA, with targets, decoys and known-site mixtures.
- Increase absent hypotheses and site alternatives on a fixed real mobility-resolved run.
- Compare prepared populations/IDs and downstream candidates/scores/quantities across page sizes.
- Verify target/decoy selection symmetry and retain unresolved localization alternatives downstream.
- Measure cold preparation and warm reuse separately, including model startup and disk construction.
- Measure full process/cgroup peaks at proposed 2/4 GiB and 2/4 workers, not extraction alone.
- Force insufficient memory/disk and require explicit incomplete/failure status, not a smaller search.
- Record manifests, commands, coherent build IDs and nonempty outputs in the WP-00 ledger.
- Numeric resource/performance targets remain proposals until these public integrated checks pass.

## Failure risks

- Full prediction needed for retrieval can dominate cold preparation despite lazy scoring predictions.
- Combinatorial modifications and decoys may exceed feasible disk/time even with bounded RAM.
- Batch-dependent deduplication or decoy construction can silently change the tested population.
- A second whole-library index in memory would invalidate the streaming contract.

## Handoff

Provide WP-05/WP-06 query examples, exhaustive mini-fixtures, ID rules and restart manifests.
Give WP-02 measured workspace minima and WP-08 provenance/coverage fields required for confidence.
Do not close acceptance on indexed lookup speed without bounded upstream preparation evidence.
