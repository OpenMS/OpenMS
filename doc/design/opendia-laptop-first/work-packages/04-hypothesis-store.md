# WP-04: Bounded hypothesis preparation and indexed family store

Status: proposed/not implemented
Owner: unassigned
Reviewer: unassigned
Suggested expertise: peptide chemistry, digestion/decoys, prediction and external-memory indexes.
Index: [Laptop-first design](../README.md)

## Goal

Prepare and query hypotheses without materializing a whole predicted library or localization graph.
Bound the entire path from FASTA or imported library through digestion, decoys and prediction.
Preserve search coverage, stable identities and provenance independently of batching and storage.

## Dependencies

- Requires [WP-01 contracts](01-contracts.md) and [WP-02 budgets](02-resource-budget.md).
- Supplies [WP-05 extraction](05-shared-extraction.md) and [WP-06 retrieval](06-candidate-retrieval.md).
- Coordinate evidence identity consumers with [WP-08](08-storage-confidence.md).
- Public acceptance closes through [WP-09](09-cli-integration.md) using [WP-00](00-acceptance-baseline.md).

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

## Small PR slices

1. Add deterministic family IDs, paged imported-library access and provenance round-trip tests.
2. Stream digestion/modifications/decoys with bounded deduplication and complete coverage checks.
3. Add budgeted prediction batches, lazy fragments and external-memory index construction.
4. Add restart validation and public cold/reuse pipeline tests across growing search spaces.

## Acceptance checks

- Unit: small exhaustive digestion/modification/decoy truth sets match streamed enumeration exactly.
- Unit: IDs, protein mappings and target/decoy associations remain stable across page/worker sizes.
- Unit: huge family pagination retains all alternatives and does not replicate shared evidence.
- Unit: absent RT/IM predictions, broad uncertainty and calibration changes preserve recovery access.
- Unit: prediction workspace estimates, sort spill, disk-full and cancellation never escape budgets.
- Unit: incompatible manifests and interrupted builds cannot be reused as complete stores.
- Integration: execute a fresh public OpenDIA preparation-to-export workflow once WP-09 is wired.
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
