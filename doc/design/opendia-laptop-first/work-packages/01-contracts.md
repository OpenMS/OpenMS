# WP-01: Contracts and scientific semantics

- **Track:** data handling.
- **Status:** in progress, private boundary, owning records and extraction-query identity implemented and validated 25 October 2026. Families/scores and full v0 handoff remain. **Owner / reviewer:** ownership is tracked in GitHub.
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

Deliver small, reviewable commits in the order below. The first three increments are
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

### 4. Families, scores and v0 handoff (pending)

- Add minimal candidate-family and score-batch records with stable references,
  unresolved/coexisting alternatives and null hypotheses. Keep identification,
  localization and quantity distinct. State candidate completeness explicitly.
- Retained scores/evidence must not depend on a destroyed raw tile or local worker.
  Do not freeze new family partitioning, null/decoy science or confidence semantics.
- Exercise the complete v0 producer/consumer path with non-IM, IM and discrete
  FAIMS examples, absent versus zero support, alternatives and malformed records.
- Publish a short contract document with field units, ownership, validators and
  schema rejection. Mark scientific freeze pending the WP-07 experiment review.
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
The tool-local [contract document](https://github.com/OpenMS/OpenMS/blob/feature/opendia-wp01-private-engine/src/topp/OpenDIA/CONTRACTS.md)
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
`clang-format` is unavailable. **Next: step 4, families, scores and v0 handoff.**
WP-01 remains open and scientific freeze remains pending WP-07 review.

## Requirement-to-check map

| Requirement | Acceptance check and current observation |
|---|---|
| One real private implementation, no public API | Native target graph and link commands show executable/test linking the same target. Build-tree and standalone tool builds work. Installation/export files contain no OpenDIA interface. |
| Behavior-preserving first move | Existing TOPP workflow tests pass with the relocated binary. Parameter/INI output matches the baseline apart from declared executable/build provenance fields and unordered protease allowed-value ordering. |
| Coordinate and missingness semantics | Step-2 production record tests and actual mzML/TSV producer/consumer checks pass, with exact counts and fields reported above. Family alternatives remain step 4. |
| No false query sharing | Step-3 production tests retain 33 distinct payloads under one constant hash after 32 field mutations, including calibration. Real-loader identity checks and actual 59-sample public extraction replay pass. Production cache integration remains WP-05. |
| Checked ownership/schema/bounds | Actual production validators and native/standalone/NDEBUG tests pass all declared rejection and ownership checks. Real-loaded evidence survives producer/tile/page destruction. |
| Family ambiguity preserved | Pending step 4: round trips must retain all declared alternatives/nulls and separate identification/localization/quantification fields. Freeze remains pending WP-07, not a local closure blocker. |

## Dependencies

- **Start / local merge:** none, including no WP-00 dependency.
- Initial v0 contracts can close before science review. Contract freeze is a later milestone, only after the independent WP-07 early localization experiment is reviewed.
- WP-02/03/04/05/08 and early WP-09 consume v0 without waiting for that freeze.

## Local merge criteria

- [x] An out-of-tree build links the tool and private tests against the same static/OBJECT implementation, with no installed OpenDIA interface.
- [ ] Round trips preserve IDs, exact coordinate references, missing versus zero, discrete FAIMS strata and unresolved family alternatives.
- [x] Query-key tests vary every scientific field, including calibration identity, and reject false sharing under forced hash collisions.
- [x] Tests reject malformed offsets, unknown schema versions and invalid ownership/lifetime use.
- [ ] A producer/consumer fixture passes using v0 records, and the contract document marks freeze pending WP-07 experiment review rather than blocking local closure.

The raw/hypothesis/evidence producer/consumer path and unstable contract document
now pass, as do complete extraction keys. The final fixture/handoff checkbox
remains open for family alternatives/nulls and scores in step 4, not for WP-07 freeze.

**Integration gates:** WP-09 owns public workflow checks, and WP-10 owns representative quality/resource/performance acceptance.
