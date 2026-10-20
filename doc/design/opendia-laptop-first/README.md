# Laptop-first OpenDIA implementation plan

- **Status:** proposed architecture and assignable backlog, not an implemented engine.
- **Updated:** 20 October 2026.
- **Owner:** unassigned. Package owners should be assigned below before implementation.
- **Suggested epic title:** `[OpenDIA] Laptop-first, mobility-native search with bounded memory and output size`.
- **Epic issue:** not created. The user will post this plan and its work packages to OpenMS.
- **Planning branch:** `rewrite/opendia`, separate from `develop`. Publishing the branch and issues is a separate action.

## Start here

- [Architecture and scientific invariants](architecture.md): the larger plan and proposed resource contract.
- [Measured baseline and limitations](evidence.md): what was actually exercised, and what remains unverified.
- [Issue posting guide](issue-posting.md): copy/paste workflow, stable package IDs and dependency tracking. No issues have been posted by this task.
- The work packages below are implementation tickets. Each includes scope, dependencies, deliverables, PR-sized steps, and acceptance checks.

This is the repository handoff for the earlier scratch design. It supersedes that design's migration-first priorities: **ion mobility and peptidoform inference belong in the first vertical slice**, not in later extensions. Existing algorithms are references or optional adapters, not architectural constraints.

## Outcome

A CPU-first DIA search engine that is fast on a resource-limited laptop even with large predicted libraries. Execute against **byte-budgeted acquisition-local raw tiles** and infer **local competing peptidoform families**, sharing sparse RT/mobility evidence. Do not retain all library chromatograms or a whole-run feature graph in memory.

Initial engineering targets are a 16 GiB laptop, a 4 GiB process-memory setting, and a 2 GiB low-memory setting. These are proposed acceptance targets, not demonstrated capabilities. End-to-end speed and scientific quality must be measured together.

**Disk and output size are first-class requirements too.** Large SQLite OSW outputs must not be replaced by equally large columnar files, duplicated caches or unbounded temporary storage. The default result artifact should scale with reported scientific results and their necessary provenance, not a materialized library-by-transition-by-time-by-mobility product.

## Output-size policy

- Default **results** profile: compact identification, localization, quantification and ambiguity summaries, with stable IDs, confidence values and provenance. No full-library copy, raw-data copy, full XICs, mobilograms or all-transition diagnostic tables by default.
- Optional **rescore** profile: preserve the explicitly supported rescoring inputs, including required rejected candidates, targets/decoys and alternatives. Results-only files must not advertise arbitrary rescoring capability.
- Optional **diagnostic** profile: selected traces, mobility patches or detailed evidence, with explicit selectors and byte limits. Full OSW compatibility export is opt-in and size-accounted, not an additional automatic output.
- Retention profiles do not change the current run's inference. Keep the full required statistical population in bounded working storage until global confidence is finalized. Never silently drop hypotheses or distinguishing ions to meet an output quota.
- Budget **final artifacts**, **working scratch**, and **reusable caches** separately, with a combined peak limit per filesystem. Count sidecars, indexes, WAL/journals, atomic staging and compaction copies. Fail explicitly if mandatory output cannot fit.
- WP-00/WP-08/WP-10 must agree per-profile byte targets and reduction goals against a capability-matched OSW baseline before acceptance. Report total bytes, peak disk, bytes per result/family/evidence row, and read/write cost. No numeric size reduction is claimed yet.

## Code ownership and placement

OpenDIA-specific workflow, family-search policy, run orchestration, retention profiles and result-schema policy belong in the proposed OpenDIA subtree:

```text
src/openms/include/OpenMS/ANALYSIS/OPENDIA/   # OpenDIA-specific interfaces
src/openms/source/ANALYSIS/OPENDIA/           # OpenDIA-specific implementation
src/topp/OpenDIA.cpp                         # thin CLI adapter only
```

`OPENDIA` follows the existing uppercase analysis-directory convention. These directories do not yet exist in this checkout and are to be introduced during implementation, not by this documentation update. Reusable primitives stay in the appropriate core OpenMS header/source trees: for example `DATASTRUCTURES`, `SYSTEM`, `FORMAT`, `PROCESSING`, `IONMOBILITY`, `CHEMISTRY` and `MATH`. Existing reusable `ANALYSIS/OPENSWATH` code stays shared.

Dependency direction is **TOPP CLI -> OpenDIA -> reusable OpenMS core**, never core -> OpenDIA. Reuse means a generic contract with meaningful use outside this workflow, not simply moving OpenDIA policy into core. New modules need registration in `src/openms/includes.cmake` and their own source/header `sources.cmake` lists, plus installed-header and consumer checks. Do not bulk-move existing APIs or fork shared algorithms as part of this plan.

## Assignable work packages

All packages are **proposed / unassigned**. Numbers identify packages, not a requirement to implement them serially.

| ID | Work package | Suggested owner expertise | Depends on | Concrete handoff |
|---|---|---|---|---|
| WP-00 | [Acceptance baseline and datasets](work-packages/00-acceptance-baseline.md) | Benchmarking + proteomics validation | None | Pinned build, fixtures, baseline manifest and quality gates |
| WP-01 | [Contracts and scientific semantics](work-packages/01-contracts.md) | C++ architecture + scientific lead | 00 | Versioned records, lifetimes, query and boundary semantics |
| WP-02 | [Global resource budget](work-packages/02-resource-budget.md) | C++ runtime/concurrency | 01 | Reservations, bounded queues, spill policy, telemetry |
| WP-03 | [Mobility-native raw tiles](work-packages/03-raw-tiles.md) | Raw formats + acquisition geometry | 01, 02 | Bounded raw access and exact-coordinate tiles |
| WP-04 | [Disk-backed hypothesis store](work-packages/04-hypothesis-store.md) | Library preparation + indexing | 01, 02 | Streamed preparation, family index and bounded prediction |
| WP-05 | [Shared sparse extraction](work-packages/05-shared-extraction.md) | Extraction kernels + numerical testing | 02, 03, 04 | Query sharing and sparse evidence patches |
| WP-06 | [Calibration, candidate retrieval and reconciliation](work-packages/06-candidate-retrieval.md) | Calibration/detection + recall evaluation | 03, 04, 05 | Calibration snapshots, high-recall families, recovery and boundaries |
| WP-07 | [Peptidoform-family scoring](work-packages/07-family-scoring.md) | Localization + quantitative proteomics | 05, 06 | Joint local evidence, ambiguous forms and quantities |
| WP-08 | [Compact outputs, evidence storage and global confidence](work-packages/08-storage-confidence.md) | Storage + statistical inference | 01, 02; final integration with 07 | Size-bounded result/rescore/diagnostic profiles and global inference |
| WP-09 | [Public CLI vertical slice](work-packages/09-cli-integration.md) | TOPP/OpenDIA integration | 03 through 08 | Opt-in end-to-end IM + peptidoform workflow |
| WP-10 | [Performance and release gates](work-packages/10-performance-release.md) | Performance + QA/release | 09 | Representative laptop validation and adoption decision |

## Parallelization and milestones

1. **M0: trustworthy baseline and contract.** WP-00 and WP-01 establish scientific policies, representative data, stable identity and build provenance. Benchmark tooling can continue evolving afterwards.
2. **M1: bounded data plane.** WP-02 establishes accounting. WP-03 and WP-04 can then proceed in parallel. WP-08 schema work starts in parallel against WP-01 contracts.
3. **M2: first scientific slice.** WP-05, WP-06 and WP-07 deliver extraction, retrieval and family scoring. WP-08 completes inference/export integration. The slice includes real mobility and localization cases.
4. **M3: public, measured workflow.** WP-09 integrates the actual CLI. WP-10 runs the agreed laptop, correctness, quality, overflow and packaging gates. Only then decide whether to replace existing defaults.

WP-09's integration skeleton can begin after contracts are agreed, but its completion depends on the upstream implementations. Synthetic producer/consumer tests unblock development, not scientific acceptance.

## Assignment and completion rules

- Put the responsible person and reviewer in the package's `Owner` / `Reviewer` fields. Keep its dependency list and status current.
- Use the package's suggested issue title and checkbox acceptance criteria. Replace stable WP dependency IDs with actual issue references only after posting; do not invent issue numbers or assume project labels exist.
- Start with the package's first PR-sized slice. Do not reinterpret a work package as a single large pull request.
- Agree contract changes with both producer and consumer owners. Record any change to scientific behavior separately from a transport optimization.
- For every acceptance check record the command, pinned input/build, output artifact and observed result. A proposed check is not a pass.
- Build locally out of tree and run relevant tests before opening a PR. Follow repository header registration, bindings and package rules when public APIs change.
- Do not introduce a new default engine until WP-10 passes. Do not claim reduced scientific scope as a performance improvement.
- Use the real public CLI and installed package where applicable. Unit tests and synthetic fixtures supplement, but do not replace, the acceptance workflow.

## Open decisions that need owners

| Decision | Accountable package | Must be resolved before |
|---|---|---|
| First supported acquisition formats and raw-cache format | 03 with 00 | M1 completion |
| Candidate completeness domain and boundary reconciliation | 01 with 06/07 | M2 quality comparison |
| Family partition and localization/null model | 07 with 01/08 | Localization acceptance |
| Exact versus approximate retrieval/model fitting | 06/08 with 00 | Any performance claim |
| Canonical evidence schema and supported export backends | 08 with 09 | M3 |
| Per-profile output byte targets, OSW comparison and retained rescoring capability | 00 with 08/10 | Before schema/performance acceptance |
| OpenDIA/core placement and one-way dependency review | 01 with every package owner | Each implementation PR |
| Numerical tolerances and scientific non-inferiority margins | 00 with scientific reviewers | Before benchmarking the new engine |
| Quantitative speed target and supported laptop/storage class | 10 with 00 | Before performance acceptance |

## Definition of done for the overall plan

The public workflow, including preparation, decoding, search, scoring, confidence estimation and export, passes the requirement matrix in [architecture.md](architecture.md). It demonstrates bounded whole-process memory, bounded peak disk and compact final artifacts, useful end-to-end speed, preserved identification/localization/quantification quality, and safe overload behavior on representative data. OpenDIA policy and reusable core code respect the ownership boundary above. Until then this is a design and implementation backlog, not a performance promise.

## Epic acceptance checklist

- [ ] WP-00 through WP-10 are posted, linked to this epic and assigned with explicit dependencies.
- [ ] Public raw-to-export workflows meet the agreed memory, runtime and per-profile output-size gates.
- [ ] Compact outputs retain current-run identification, localization, ambiguity and quantification semantics; optional rescoring capabilities are verified separately.
- [ ] Peak disk accounting covers scratch, caches, final artifacts and finalization overlap, with explicit failure and cleanup tests.
- [ ] OpenDIA-specific code lives under `ANALYSIS/OPENDIA`, reusable core has no dependency on it, and the TOPP entry point remains thin.
- [ ] Real-data, scientific-quality, installed-consumer and compatibility evidence is attached before default adoption.
