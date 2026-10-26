# WP-01: Contracts and scientific semantics

- **Track:** data handling.
- **Status:** locally complete and final adversarial review passed 25 October 2026. **Integration deliberately deferred by the owner on 26 October 2026.** Continue dependent packages on stacked branches/PRs until the combined refactor is demonstrably correct, fast and memory-efficient. Earlier merge-ready notes below describe local correctness only and do not authorize a `develop` merge. Issue remains open. Scientific freeze remains a separate WP-07 milestone. **Owner / reviewer:** ownership is tracked in GitHub.
- **Parent epic:** [#10661](https://github.com/OpenMS/OpenMS/issues/10661). **Issue:** [#10663](https://github.com/OpenMS/OpenMS/issues/10663).
- **Shared policies:** [architecture](../architecture.md) and [plan](../README.md).

## Objective

Provide small, explicitly unstable v0 interfaces that unblock independent adapters before scientific contracts are frozen.

## Scope

- Define raw tiles, paged hypotheses, complete extraction keys, sparse evidence, candidate families and score records.
- Include IM coordinates, absent mobility, family alternatives and nulls from day one, without requiring finished IM/family science for the first runnable engine.
- Keep OpenDIA code tool-private under `src/topp/OpenDIA/{OpenDIA.cpp,include,source,test,CMakeLists.txt}`.
- Use one private static/OBJECT library shared by the tool and tests, without installed OpenDIA headers or a new core API.
- Specify ownership, units, stable identities, calibration snapshot identity and schema-version rejection, not a general capability-manifest framework.

## Implementation roadmap

Deliver small, reviewable commits in the order below. All four increments are
implemented on `feature/opendia-wp01-private-engine`, based on clean upstream
`3e47bc7e7468656d0d53237f11fc4f049d278efe`, with standalone, full native
in-tree and installed-SDK workflow validation complete. Keep implementation
commits separate from this scratch-only planning branch.
Do not include the optional WP-00 profiler.

### 1. Private engine boundary (implemented and locally validated)

- Relocate the current tool into `src/topp/OpenDIA/` without changing its CLI,
  parameter defaults, scientific algorithms or output policy.
- Put the existing orchestration behind one private STATIC/OBJECT target. Both
  the thin executable and tool-private tests must use that same implementation.
- Retain standalone tool builds against the OpenMS CMake package. Do not install
  or export an OpenDIA library/header or add a core/pyOpenMS API.
- Check the shared entry point, parameter/INI compatibility and invalid CLI paths.
  Build out of tree and rerun the existing public OpenDIA TOPP producers,
  SQLite/Parquet comparisons, predicted mode and failure-path tests.

### 2. Coordinates, identities and ownership (implemented and locally validated)

- Add minimal v0 raw tile, hypothesis-page and sparse evidence value records.
  Document units and exact coordinate references, stable run/precursor/family IDs,
  calibration identity and explicit schema version.
- Represent absent mobility, measured zero versus missing evidence, and discrete
  FAIMS strata without substituting lookup bins for measurement coordinates.
- Start with bounded owning containers and explicit checked views/leases. Avoid
  raw pointers that can outlive their owner and avoid a general storage framework.
- Tests transfer records between an actual producer and consumer, retain owned
  data across producer destruction, and reject invalid offsets, dimensions,
  schema versions and invalid views. Round trips here need not invent a wire format.

### 3. Complete extraction query identity (implemented and locally validated)

- Define a value key covering run/lane and acquisition selection, MS level and
  isolation geometry, RT/mz/IM bounds and units, discrete FAIMS selection,
  tolerance, edge/aggregation policy, isotope policy and calibration snapshot.
- Use exact equality on canonical scientific values. A hash is an index hint,
  not proof that two queries can share evidence. Reject nonfinite/invalid inputs.
- Vary every field independently in tests. Force hash collisions and prove that
  unequal queries do not share. Calibration changes must invalidate equivalence.
- Do not implement a cache or extraction algorithm in this increment.

### 4. Families, scores and v0 handoff (implemented and locally validated)

- Add minimal candidate-family and score-batch records with stable references,
  unresolved/coexisting alternatives and null hypotheses. Keep identification,
  localization and quantity distinct. State candidate completeness explicitly.
- Retained scores/evidence must not depend on a destroyed raw tile or local worker.
  Do not freeze new family partitioning, null/decoy science or confidence semantics.
- Exercise the complete v0 producer/consumer path with non-IM, IM and discrete
  FAIMS examples, absent versus zero support, alternatives and malformed records.
- Keep the [contract narrative](../contracts-v0.md) with field units, ownership,
  validators and schema rejection on this planning branch, not in the implementation
  branch or its pull request. Mark scientific freeze pending WP-07 experiment review.
  Hand v0 to WP-02/03/04/05/08 and early WP-09 without waiting for that freeze.

## Observed starter implementation, 25 October 2026

Implementation commit [e7e67276ce](https://github.com/OpenMS/OpenMS/commit/e7e67276ce099a92e08780b9f7235863d0f9b891)
is pushed on the separate implementation branch. The production-source move is
99% identical: one private-header include and `main` renamed to `OpenDIA::run`.
A thin main and boundary test link the same `OpenDIA_impl` static archive. No
scientific changes, instrumentation, installed headers or new public API were added.

Observed local checks, independently repeated by the coordinator:

- Standalone Debug tool/test build against the matching real OpenMS build-tree
  CLI package passed. The shared-entry test passed help and invalid-option calls
  sequentially in one process.
- All **28 existing OpenDIA TOPP regressions passed** with the new executable,
  including real SQLite/Parquet production and identity checks, predicted-library
  materialization/end-to-end, preconditions and failure/cleanup paths. Replay
  retained the existing generated CTest definitions, changing only the executable
  and writable fixture-output paths. Frozen core/runtime dependencies were read-only.
- Complete generated INI XML agrees with the baseline after sorting only the
  protease allowed-value list. Raw byte comparison differed in that list's order,
  not its values, defaults or any other attributes. No broader normalization was used.
- Link rules and archive symbols confirm both consumers share `OpenDIA::run`.
  Standalone install produced an empty manifest. Git whitespace checks passed.

**Deeper native and installed-consumer validation:** a fresh root OpenMS Debug
build compiled the real core, CLI layer, production tool and automatically
registered private test, without importing a substitute OpenMS core. All **32
selected native CTests passed**: the 28 public workflows, shared-entry test,
INI generation/section check and real tool-registry readability. Generated test
backtraces identify `openms_add_tool_class_test`, and runtime library bindings
point to the freshly built native libraries.

This exposed and fixed one missed integration edge: the installed-consumer fixture
still referenced the removed `src/topp/OpenDIA.cpp`. With the original declaration,
a real configure against the newly installed SDK failed with `Cannot find source
file: ../../topp/OpenDIA.cpp`. Fix [094bf09ef8](https://github.com/OpenMS/OpenMS/commit/094bf09ef8)
adds the production tool subdirectory instead. The same configuration then built
and its existing `TestExternalOpenDIA` INI test passed. The full real SDK fixture
passed **12 setup/configure/build/run CTests**, including **all six nested consumer
tests** for CLI, Qt, OSWPQ, OpenDIA, core code and PeptDeep. Installed headers and
CMake exports contain no private OpenDIA interface or archive.

The native configuration used ONNX enabled and GUI, opentims and Thermo RAW disabled.
Models were staged into the new tree. Existing Eigen/HiGHS and ONNX dependencies
were reused read-only, with ONNX runtime still resolved from that dependency
location. This verifies installed-SDK consumption, not a self-contained relocated
runtime bundle. The first native build invocation hit the tool's 600-second cap,
and its incremental continuation completed successfully.

**Still outside this check:** other platforms and optional RAW/opentims build
configurations. `clang-format` is not available locally. No implementation PR has
been opened. This mechanical move did not deliver coordinates/ownership. The
separate step-2 increment below now does. Step 3 below adds extraction keys.
Families and scores remain the final increment.

## Observed owning-record increment, 25 October 2026

Implementation [919681d30b](https://github.com/OpenMS/OpenMS/commit/919681d30b)
adds concrete `RawTile`, CSR `HypothesisPage` and sparse `EvidenceBlock` records
to the same private `OpenDIA_impl` archive. Follow-up
[6d4801f95a](https://github.com/OpenMS/OpenMS/commit/6d4801f95a) preserves exact
native acquisition IDs and adds real-loader interoperability coverage. Both are
on `feature/opendia-wp01-private-engine`, separate from this planning branch.
The planning-only [contract document](../contracts-v0.md)
specifies units, caller identity namespaces, immutable owning copy leases,
count limits, sparse missingness and schema rejection. Scientific freeze remains
pending WP-07 review.

Observed checks over the complete step-2 result:

- Actual native ClassTest checks pass copies, moves, assignment, producer
  destruction, detached value independence, empty records, exact limits and
  invalid access. Production validators reject unknown schemas/enums, nonfinite
  values, malformed/overflowing CSR offsets, reversed geometry and exceeded caps.
  Explicit missing evidence, measured zero, absent mobility, both continuous IM
  units, zero FAIMS and negative baseline-corrected intensity remain distinct.
- `OpenDIARecordsInterop_test` loads existing mzML files using the real public
  `MzMLFile`, destroys the producer experiment, and compares **every retained raw
  point** against an independent load. Observed counts are **513 non-IM points
  from 114 spectra**, **32,750 IM points from 150 spectra**, and **44,354 discrete
  FAIMS points from 19 spectra**. RT/mz/intensity, caller frame/scan IDs, source
  native ID, MS level, isolation edges, continuous IM units/values and FAIMS
  presence/voltage compare exactly. The source RT conversion to seconds is kept.
- Sparse evidence retains **257 / 16,375 / 22,177** selected observations after
  its originating tile and hypothesis page are destroyed, with exact identities,
  coordinates and values. This tests transfer/ownership, not scientific matching.
- The actual `TransitionTSVFile` loads the workflow-17 library. A fixed caller
  identity registry, independent of page position, preserves **six hypotheses
  and 18 transitions**, their mz/intensity, target/decoy state and provenance.
  Normalized library RT and unit-unspecified IM are deliberately absent from
  physical prediction fields, not silently relabelled as seconds or known IM.
  Original metadata remains in the library store.
- **34/34 selected native CTests pass**, including all 28 existing public OpenDIA
  workflows, entry/INI/registry checks and both new record tests. Standalone Debug
  SDK build/tests pass **3/3**. A separate final `-DNDEBUG` build passes both
  record tests **2/2**, so rejection and test checks do not disappear with asserts.
- The real installed-SDK fixture again passes **12/12** setup/configure/build/run
  CTests and **6/6** nested consumers. No private OpenDIA headers, archive or CMake
  interface are installed/exported. Native archive/link inspection confirms the
  tests use the production validators and real class-test helper, not copied code.
- The first real-loader run correctly rejected the assumption that workflow-17
  mzML was non-IM. The non-IM case now uses existing workflow-1, and the acquisition
  coverage assertion remains active. Legacy library label-group and FAIMS-file
  metadata warnings are emitted by the real loaders, not hidden or treated as a
  changed scientific result. Git whitespace checks pass, no new compiler warnings.

**Limits of this increment:** count caps are not byte accounting or run-level
memory admission. The by-value record layout is not a packed performance design.
The real-loader test is a boundary fixture, not a production WP-03/04/05 adapter,
streaming implementation or end-to-end scientific extractor. No existing workflow
algorithm/default changed. Other platforms, optional RAW/opentims configurations
and self-contained runtime relocation remain untested. `clang-format` is unavailable.
No implementation PR is open. Step 3 below now adds complete extraction-query
identity. Step 4 still owns alternatives/nulls, scores and the complete v0 handoff.

## Observed extraction-query increment, 25 October 2026

Implementation [b71ce465a7](https://github.com/OpenMS/OpenMS/commit/b71ce465a7)
adds one checked `ExtractionQuery` value in the existing private implementation.
Follow-up [cfbb419390](https://github.com/OpenMS/OpenMS/commit/cfbb419390)
adds real-loader and actual public-extractor interoperability checks.
Its detached specification covers caller context/calibration, opaque acquisition
selection, MS level, isolation geometry, RT/mz/IM intervals and edge policies,
mobility units, discrete FAIMS, tolerance, aggregation and isotope policy.
Only acquisition-set order/duplicates are canonicalized. Physical values are
finite and exact, with signed zero equal and adjacent doubles distinct.
Schema 0 is the only supported version. The contract specifies a fixed tophat
profile, full-width tolerance units, absent/all versus present-empty selection,
and producer responsibility for canonical physical bounds.

Observed checks over the complete result, including all earlier WP-01 changes:

- The production key test retains **33 distinct payloads under one constant
  hash**, covering **32 independent valid field mutations**, including changed
  calibration identity. All payload lookups remain distinct. Reordered/duplicate
  acquisition IDs and signed zero compare equal, without overwriting the original
  payload. Missing axes/FAIMS differ from present zero. Caller/getter mutations
  cannot change a retained key.
- Constructors reject unknown schema/enums, nonfinite values on every interval
  endpoint, FAIMS and tolerance, reversed intervals, negative tolerance, zero
  MS level and ambiguous isotope policy. All four edge policies and zero-width
  requests remain explicit. These checks also pass with `-DNDEBUG`.
- Actual mzML-loaded non-IM, continuous IM and discrete FAIMS records preserve
  query coordinates, units, acquisition IDs and context. Changing only the
  calibration snapshot produces an unequal key. Existing full raw/evidence/TSV
  ownership checks still pass with the counts reported for step 2.
- A query round trip through the actual public `ChromatogramExtractorAlgorithm`
  preserves **all 59 time/intensity samples exactly**, with measured nonzero
  signal, on its existing mzML fixture. This exercises one supported all-RT,
  all-acquisition, non-IM SUM request, not a new production adapter.
- **35 selected native CTests pass**, including all **28 public OpenDIA TOPP
  workflows**, all four private tests, INI and registry checks. The actual
  installed-SDK fixture passes **12/12**, for **47/47 combined**, plus **6/6**
  nested consumers. Standalone SDK tests pass **4/4** and separate final NDEBUG
  record/query tests pass **3/3**. No private headers/archive/interface are
  installed or exported. All use the same production `OpenDIA_impl` target.

There is no cache, production hasher, new extraction algorithm, dependency,
scientific/default change or performance claim. Equality conservatively retains
the whole request, including tolerance even when physical bounds coincide.
Consumers must honor every field or reject unsupported requests. Production
sharing/calibration adapters belong to WP-05. Optional platforms/RAW configurations
and self-contained runtime relocation remain outside this Linux check, and
`clang-format` is unavailable. Step 4 below now completes the local v0 handoff.
WP-01 remains open for implementation review/integration, not for scientific freeze.

## Observed family/score increment and complete v0 handoff, 25 October 2026

Implementation [3a8b2a64db](https://github.com/OpenMS/OpenMS/commit/3a8b2a64db9ad9b311a327f0f3d60eb5deeea307)
adds concrete immutable `FamilyBatch` and `ScoreBatch` records to the same private
archive. Follow-up [49deb7c5c5](https://github.com/OpenMS/OpenMS/commit/49deb7c5c58e03b92445c73e22ff0a37e28a3cfd)
checks the complete real-loader handoff and genuine public CLI scored outputs.
No scientific algorithm, production adapter, public API or dependency was added.

- Family records retain stable member/alternative IDs, per-member target/decoy,
  optional signed charge/adduct, checked local query region, sparse shared versus
  discriminating fragment incidence, provenance and explicit candidate completeness.
  Singleton, competing, coexisting and empty-member null alternatives survive
  without selecting a winner. Empty alternatives and explicit null differ.
- Score records retain independent optional identification/localization raw scores
  and explicitly unit-tagged quantity, plus ambiguity, method provenance and checked
  family/alternative/evidence references. Missing is not zero. Raw scores are not
  probabilities. Immutable copy leases keep families/evidence alive independently
  of raw tiles, pages, producers and parent-score lifetime.
- Production native/standalone/NDEBUG checks preserve zero and full uint64 IDs,
  negative and positive charges, absent constraints, all completeness declarations,
  zero/missing scores/support and detached-value independence. They reject invalid
  schemas/enums, duplicate identities, dangling references, context/calibration
  mismatch, nonfinite scores, empty present units, exceeded count limits, bounds
  errors and moved-from leases. Permitted repeated evidence references are retained.
- The complete mzML/library fixture independently verifies all **513 non-IM,
  32,750 IM and 44,354 FAIMS raw points**. After producer destruction each mode
  retains **four unresolved alternatives/scores**, including null and coexistence,
  and **259 / 16,377 / 22,179 evidence rows**, respectively. Those totals include
  two explicitly declared missing/zero support rows in addition to the earlier
  sparse measured observations. These are boundary examples, not inferred families.
- Two registered public TOPP checks consume actual SQLite and Parquet OpenDIA
  `results.tsv` outputs through the real `CsvFile` API. Each independently re-read
  export retains **four rows, including two decoys**, with exact run/feature/precursor
  IDs, charge, RT/mz, genuine `d_score` and integrated intensity. Localization stays
  absent, candidate completeness UNKNOWN, and added null alternatives unscored.
  A filtered export is not treated as the complete scientific candidate population.
- The final whole-result replay on **25 October 2026, 13:46:17–13:47:10 UTC** rebuilt
  all affected native targets and passed **38 native plus 12 actual installed-SDK
  CTests, 50/50 combined**, **6/6 nested consumers**, **5/5 standalone SDK tests**
  and **4/4 explicit NDEBUG tests**. Both genuine public-export handoffs also pass
  under NDEBUG. All **28 existing** public OpenDIA workflows and the **two new**
  handoff workflows passed, with no skipped/not-run tests. Earlier 33-payload
  forced-collision checks and the 59-sample actual extractor replay still pass.
- Installed headers/archives/exports still contain no private OpenDIA interface.
  Native class-test helper/production archive and actual native/installed runtime
  bindings are verified. No new compiler warnings, whitespace errors or changes to
  frozen baseline, third-party or vcpkg trees.

All **five local merge criteria** are now satisfied on the implementation branch.
WP-01 remains open for upstream integration after the final adversarial review below, with no implementation PR opened.
Downstream WP-02/03/04/05/08 and early WP-09 can consume the documented **unstable v0**
now. Start the next implementation increment with WP-02 resource accounting and
admission, not more record scaffolding. Count caps are not byte budgeting, packed
streaming, runtime admission or a memory/performance improvement. Production raw,
hypothesis, extraction and family/scoring adapters remain their downstream work.
The Linux configuration still excludes GUI/opentims/Thermo, other platforms and
self-contained runtime relocation. `clang-format` remains unavailable. **Scientific
family/null/decoy/confidence freeze stays pending WP-07 experiment review**, not
an extra local v0 closure requirement.

## Final adversarial review (25 October 2026)

Reviewed the entire WP-01 change set from baseline
`3e47bc7e7468656d0d53237f11fc4f049d278efe` through
[7dcbcea35e](https://github.com/OpenMS/OpenMS/commit/7dcbcea35e85c10e3c6291c6a7134668e83a2207),
including every production validator, ownership/query/family/score contract,
existing caller, native/standalone target and installed-consumer boundary.
**No blocking production defect was found.** The orchestration relocation remains
an include/entry-point change, not a scientific algorithm rewrite.

The concrete test gap was multi-family namespace and cumulative-cap coverage.
The test-only commit above adds 17 executed assertions: non-positional family IDs
including UINT64_MAX, family-local alternative IDs, cross-family member/evidence
rejection, all context axes, repeated evidence references and exact cumulative
family/score cap boundaries. No production abstraction or dependency was added.

The predeclared whole-result acceptance loop was then run over the final commit:

| Requirement / boundary | Observed final result |
|---|---|
| Real native tool, regressions and installed SDK | Fresh build and **50/50 CTests** passed on 25 October 2026, **14:26:19–14:26:45 UTC**: 38 native plus 12 installed-consumer fixture steps, including all 28 existing TOPP workflows and both genuine score handoffs. Six nested installed consumers passed again independently. |
| Standalone SDK and release-style validation | **5/5 standalone** and **4/4 explicit NDEBUG** tests passed. Both actual SQLite/Parquet exports also passed the NDEBUG handoff. |
| New adversarial family/reference/cap checks | The strengthened test passed in native, standalone, NDEBUG and sanitizer builds. A family cannot borrow a member or evidence row from another family. Cumulative family bounds accept 22 and reject 21, and score references accept 3 and reject 2. |
| Owning raw-to-score handoff | Real non-IM/IM/FAIMS loaders retained exactly **513 / 32,750 / 44,354 raw points** and **259 / 16,377 / 22,179 evidence rows**, with four unresolved alternatives/scores after producer destruction. Missing/zero, null and coexistence semantics remain explicit. |
| Real extraction and scored public outputs | The real public extractor retained **59 samples** exactly. Each independently consumed SQLite/Parquet CLI export preserved **four genuine rows, two decoys**, exact d_score/intensity and absent localization. Query tests retained **33 distinct payloads** under one forced hash. |
| Instrumented actual production archive | A separate GNU 13 ASan/UBSan build instrumented the real private implementation and tests, not copied source or a mock. **5/5 tests** plus both genuine CLI output handoffs passed with no sanitizer diagnostics or suppressions. OpenMS core/external shared dependencies were not transitively instrumented, so this is not whole-stack sanitizer coverage. |
| Unchanged CLI parameters and capabilities | A fresh comparison initially exposed different SDK build dates and enabled optional RAW readers. Rebuilding the actual relocated tool against the **same frozen baseline SDK** resolved that comparison mismatch. Every INI option/default/attribute then matched, sorting only unordered enzyme allowed-value order. Canonical SHA-256: `3c080bc948989187266f91181a89acdd43f3c5a5b70baf62d0983d9ded600fbd`. Capability/provenance differences were not normalized away. |
| Private packaging and real linkage | Native/installed runtime bindings and native production-archive linkage passed. No private OpenDIA header/archive/target leaked into installed headers or CMake exports. Header-only include compilation and whitespace checks passed. |

These checks close the local review loop for all five merge criteria below.
The implementation branch passed local correctness review. Integration is deferred under the 26 October stacked-validation policy, not authorized by these tests.
No Windows/macOS, GUI/opentims/Thermo native configuration, fully self-contained
runtime relocation or whole-stack sanitizer result is claimed. `clang-format`
was unavailable. This review does not approve new scientific acceptance margins,
prove byte-admission/resource targets, freeze contracts or reopen WP-00.
WP-02 accounting/admission is the next implementation increment.

## Requirement-to-check map

| Requirement | Acceptance check and current observation |
|---|---|
| One real private implementation, no public API | Native target graph and link commands show executable/test linking the same target. Build-tree and standalone tool builds work. Installation/export files contain no OpenDIA interface. |
| Behavior-preserving first move | All 28 existing TOPP workflow tests pass with the relocated binary. Against the identical baseline SDK, every parameter/INI option/default/attribute matches after sorting only unordered enzyme allowed-value order. Different optional-reader configurations are not compared as equivalent. |
| Coordinate and missingness semantics | Step-2 production records and full step-4 mzML/TSV family/score handoffs pass with exact coordinates/IDs, missing/zero support and all alternatives retained. |
| No false query sharing | Step-3 production tests retain 33 distinct payloads under one constant hash after 32 field mutations, including calibration. Real-loader identity checks and actual 59-sample public extraction replay pass. Production cache integration remains WP-05. |
| Checked ownership/schema/bounds | Actual production validators and native/standalone/NDEBUG/sanitizer tests pass declared rejection, ownership, family-local reference and cumulative-cap checks. Real-loaded evidence survives producer/tile/page destruction. |
| Family ambiguity preserved | Step-4 production and real-loader/public CLI checks preserve declared singleton/competing/coexisting/null alternatives, completeness and independent identification/localization/quantity. Freeze remains WP-07, not a local v0 blocker. |

## Dependencies

- **Start / local merge:** none, including no WP-00 dependency.
- Initial v0 contracts can close before science review. Contract freeze is a later milestone, only after the independent WP-07 early localization experiment is reviewed.
- WP-02/03/04/05/08 and early WP-09 consume v0 without waiting for that freeze.

## Local merge criteria

- [x] An out-of-tree build links the tool and private tests against the same static/OBJECT implementation, with no installed OpenDIA interface.
- [x] Round trips preserve IDs, exact coordinate references, missing versus zero, discrete FAIMS strata and unresolved family alternatives.
- [x] Query-key tests vary every scientific field, including calibration identity, and reject false sharing under forced hash collisions.
- [x] Tests reject malformed offsets, unknown schema versions and invalid ownership/lifetime use.
- [x] A producer/consumer fixture passes using v0 records, and the contract document marks freeze pending WP-07 experiment review rather than blocking local closure.

The complete raw/hypothesis/query/family/evidence/score producer/consumer handoff
and unstable contract document now pass, as do both genuine public CLI output
handoffs. All five criteria and the final adversarial review are locally verified.
Upstream integration remains, not further WP-01 record implementation or the
separate WP-07 scientific freeze.

## Next WP-01 action

Open an implementation-only pull request from `feature/opendia-wp01-private-engine`
against `develop`, obtain upstream review and passing CI, then integrate and close
#10663 / move its project item to Done. No PR or merge has been performed yet.
Keep this master plan and the [v0 contract narrative](../contracts-v0.md) on the
planning branch, outside that PR. No additional v0 record scaffolding is needed.
WP-02 resource accounting/admission is the next development package and can use v0
without waiting for the separate WP-07 scientific freeze.

**Integration gates:** WP-09 owns public workflow checks, and WP-10 owns representative quality/resource/performance acceptance.
