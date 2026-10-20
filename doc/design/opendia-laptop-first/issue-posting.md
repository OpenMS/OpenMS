# Posting the plan as OpenMS issues

[Master plan](README.md) | [Architecture](architecture.md)

**Status:** issue-ready Markdown only. No GitHub issues, project items, labels or assignments have been created by this documentation task.

The documents live on the local `rewrite/opendia` branch, not `develop`. Push/publish that branch when desired before adding remote document permalinks. This task has not pushed it.

## Posting workflow

1. Create one epic using the title and standalone body below. The checked-in OpenMS Feature Request template has an **Enhancement Description** field; the body can be pasted there.
2. Create one issue per work-package file, using its **Suggested issue title** and the file contents as its body. The packages contain owner/reviewer fields, dependency IDs, code placement, output-size obligations, PR-sized steps and unchecked acceptance criteria.
3. Replace `Parent epic` and `Issue: not created` with the actual references after creation. Replace dependencies such as `WP-02` with their real issue references, keeping the WP IDs for cross-document traceability. Add the issues to the intended OpenMS project and assign owners/reviewers there.
4. Update the epic checklist and repository metadata with the real references. Never infer issue numbers, project fields, available labels or assignment permissions from this plan.

Package bodies intentionally use readable WP IDs and repository paths rather than relative Markdown links, because such links would resolve incorrectly when pasted into an issue. The master index retains clickable local links for repository browsing. After these documents are available on the remote, optionally add their actual permalinks to the issues. A local commit is not evidence that a remote link is already live.

Keep the repository design and issue decisions synchronized. Record scope changes and scientific/output-capability decisions in both rather than letting an issue silently supersede an invariant. Split PR-sized steps into sub-issues if useful, but retain the package's acceptance gate and dependency owner.

## Suggested epic title

`[OpenDIA] Laptop-first, mobility-native search with bounded memory and output size`

## Standalone epic body to paste

```markdown
## Goal

Build a CPU-first OpenDIA engine for resource-limited laptops and large predicted libraries. Use byte-budgeted acquisition-local tiles for execution and local competing peptidoform families for inference. Ion mobility and localization belong in the first vertical slice.

Initial process-memory targets are 4 GiB and a 2 GiB low-memory mode on a 16 GiB laptop. These are proposed targets, not implemented guarantees. Benchmark total cold/warm runtime at matched search space and scientific quality.

## Output-size requirements

- Default results are compact identification/localization/quantification summaries with ambiguity, confidence and provenance, not whole-library traces or transition cubes.
- Optional rescore retention preserves the complete population and evidence needed by explicitly supported rescoring operations. Results-only artifacts must not claim arbitrary rescoring capability.
- Optional diagnostics are selected and byte-capped. Full OSW compatibility export is opt-in and size-accounted.
- Keep required candidates/decoys/alternatives through global confidence. Retention settings must not silently change the current run's science to meet a byte limit.
- Budget final artifacts, temporary workspace and reusable caches separately, with a combined peak limit per filesystem, including sidecars, indexes, journals and staging/compaction overlap.
- Agree numerical per-profile size ceilings and capability-matched reduction goals against OSW before acceptance. Compression or a backend switch alone is not success.

## Code placement

Introduce OpenDIA-specific interfaces under src/openms/include/OpenMS/ANALYSIS/OPENDIA/ and implementation under src/openms/source/ANALYSIS/OPENDIA/. These are proposed new directories using existing uppercase analysis-module conventions.

Keep generic containers, I/O, resource accounting, chemistry, numerical and mobility kernels in the appropriate reusable OpenMS core modules. Dependency direction is CLI -> OpenDIA -> core, never the reverse. src/topp/OpenDIA.cpp remains a thin adapter. Preserve existing reusable OPENSWATH APIs and register/test new public headers normally.

## Child issues

Replace the placeholders with real issue references after posting:
- [ ] WP-00: acceptance baseline, data and output-footprint targets (issue not created)
- [ ] WP-01: contracts and scientific semantics (issue not created)
- [ ] WP-02: global memory/disk admission (issue not created)
- [ ] WP-03: mobility-native raw tiles (issue not created)
- [ ] WP-04: bounded hypothesis preparation/store (issue not created)
- [ ] WP-05: shared sparse extraction (issue not created)
- [ ] WP-06: calibration, candidate retrieval and reconciliation (issue not created)
- [ ] WP-07: peptidoform-family scoring (issue not created)
- [ ] WP-08: compact outputs, retained evidence and global confidence (issue not created)
- [ ] WP-09: public CLI integration and output profiles (issue not created)
- [ ] WP-10: laptop runtime, output-size and release gates (issue not created)

## Acceptance

- [ ] Representative public raw-to-export workflows meet agreed memory, runtime and per-profile output-size gates.
- [ ] Identification, localization, unresolved/coexisting forms and quantification pass separate scientific-quality checks, including joint mobility/localization data.
- [ ] Reopening, export and supported rescoring verify each retained artifact's declared capabilities.
- [ ] Peak disk and failure tests include output finalization, disk exhaustion, cancellation and incomplete-artifact protection.
- [ ] OpenDIA/core ownership and installed consumers pass their checks, with no new default until release gates pass.

## Repository references and current evidence

Master: doc/design/opendia-laptop-first/README.md
Planning branch: rewrite/opendia (publish before adding remote permalinks)
Architecture: doc/design/opendia-laptop-first/architecture.md
Evidence: doc/design/opendia-laptop-first/evidence.md
Tasks: doc/design/opendia-laptop-first/work-packages/

The proposed engine and compact-output profiles are not implemented. Existing small public extractor/IPF tests are baseline evidence only. No output-size reduction or large-library whole-engine speed/memory guarantee has been demonstrated. Known build/test blockers and remaining acceptance work are recorded in the evidence report and WP-00.
```

## Dependency checklist for project setup

- WP-00 establishes the baseline and feeds WP-01.
- WP-02 follows WP-01; WP-03 and WP-04 then proceed in parallel.
- WP-05 requires WP-02/03/04; WP-06 requires WP-03/04/05; WP-07 requires WP-05/06.
- WP-08 schema/writer work starts after WP-01/02 and integrates WP-07 before completion.
- WP-09 integrates WP-03 through WP-08. Its CLI skeleton may start earlier against reviewed contracts.
- WP-10 completes after WP-09, although benchmark preparation starts with WP-00.

## Evidence to attach when closing an issue

Include a coherent local build identity, real public commands and pinned inputs, observed outputs and failure paths, per-requirement checks, scientific/profile compatibility results, memory and total-disk measurements where applicable, and source-placement review. A checked task list alone does not establish implementation success. Attach or link large evidence artifacts through the project's approved storage rather than uploading raw acquisitions into issue bodies.
