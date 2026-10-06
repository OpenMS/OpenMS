# IdentificationData implementation validation

Measured on 5 October 2026 against OpenMS base commit `fe7b9cc` (Linux x86-64).
The owning model, native format, adapters and Python bindings are implemented.

## Direct OMS SQLite storage (2026-10-06)

OMS schema 6 now stores owning identifications directly in typed SQL tables.
It no longer writes temporary Parquet files or embeds a native archive. This
replaces the unreleased schema-6 experiment in place; released schemas 1–5 keep
the compatibility reader. Quantitative maps and their stable ID links remain
supported. Runs share the checked ordered score schema and numeric score columns.
Inference, custom enzymes/modifications, group arrays and metadata remain queryable.
Metadata names are interned, and `ID_MetadataValue` exposes a readable SQL view.

Final local validation:

- GCC 13.3 Release core builds with strict JSON conversions disabled.
- All **21 focused C++ suites pass**, including released OMS files, feature and
  consensus maps, scoreless catalogs and transactional I/O failures.
- The full Python suite reports **6,620 passed, 29 skipped, 10 xfailed, 7 xpassed**,
  with no failures; this includes **24 owning-model regressions**.
- All five changed C++ translation units pass Clang 18.1.3 syntax compilation
  with OpenMP enabled. Native Windows/macOS execution remains a CI responsibility.
- Added checks cover direct SQL score filtering and editing, corrupted foreign
  keys, full-width unsigned record IDs, vector order, empty observations/sources,
  absent versus empty parent catalogs, complete inference values, modifications,
  group arrays, exact NaN payloads and signed zero. The rich inference fixture
  compares native Parquet and direct SQL round trips.

Release comparison: **one million PSMs, 1,000 runs, 200 protein hits per run**.
Medians of three fresh-process writes and warm-cache reads; no competing local
builds/tests during the final measurements. Both variants use the same owning
model and benchmark executable; the baseline library is from `45f09524`.

| OMS implementation | Write (s) | Load (s) | File (MB) | Load peak RSS (MiB) |
|---|---:|---:|---:|---:|
| Previous embedded Parquet | 4.43 | 4.50 | 14.96 | 1358 |
| Direct SQLite | 18.44 | 10.43 | 494.98 | 1260 |

The direct representation is slower and larger than compressed Parquet; this is
not a bulk-I/O speed improvement. It makes the data available to ordinary SQL,
removes extraction/temporary-file I/O, and modestly reduces full-load peak RSS.
These synthetic data are highly repetitive and favor Parquet compression.
All runs preserve the same digest over sequence, observation, charge, score, RT
and m/z. Times exclude generation, conversion and digest verification. Write
process RSS includes generation/conversion and is not isolated writer memory.

The first SQL layout was 805.91 MB. SQLite page accounting showed that repeated
metadata names and overlapping metadata indexes dominated it. Name interning and
reuse of the uniqueness index reduced the final file to 494.98 MB. The final
5.4 million metadata rows and their uniqueness index still occupy 275.22 MB,
about 56% of the file. Each is a typed SQL row; the remaining work is ordinary
SQLite insertion/indexing overhead, not Parquet encoding. These measurements do
not establish cold disk performance, arbitrary query latency or billion-PSM scale.

Raw observations, source blob IDs, commands and table/index sizes are in
[`oms-sqlite.json`](../../tools/benchmarks/identification_data/results/oms-sqlite.json).

## CI portability and Python imports (2026-10-06)

The CI failures at `179705666f4adc0e1eb8ec4156a6314b8180aaa3` had four causes:

- Native manifest parsing relied on implicit JSON conversions, disabled by the
  CI dependency configuration. String and boolean reads now use explicit typed
  extraction.
- MSVC does not guarantee nonthrowing moves for the trees inside optional adducts
  and processing metadata. Match edits now stage replacement vectors on platforms
  with throwing payload moves, then commit through nonthrowing vector swaps.
  Processing metadata has unique ownership and commits through a pointer swap.
  Platforms with nonthrowing match moves retain the in-place fast path.
- Clang 18.1.3 crashed while instantiating dependent `requires` expressions in
  recursive feature lambdas. Explicit feature-type checks retain subordinate
  traversal and compile with that compiler.
- The format Python module eagerly converted the default `ProgressLogger::NONE`
  argument before the misc module registered its enum. The optional logging
  argument now defaults to `None` and resolves to `NONE` when called.

Validation on Linux x86-64, GCC 13.3, C++23, Release `-O3 -DNDEBUG`,
Arrow/Parquet 25 and `JSON_USE_IMPLICIT_CONVERSIONS=0`:

- The complete core library rebuilt and all **21 focused C++ suites passed**.
- Both original adapter/converter translation units reproduced Clang 18.1.3's
  frontend crash (exit 139); both fixed units pass syntax compilation.
- A fault-injection harness exercised every allocation failure before successful
  filtering, replacement, transformation and processing-metadata replacement.
  All 34 failure points on the fast path and 75 on a forced staged-vector path
  preserved the original run. The latter exercises the fallback on Linux; it is
  not a native Windows build.
- All actual Python extension sources were compiled and linked using the
  generated Release build commands, with nanobind 3.1 and Python 3.12. The original
  format module reproduced `ImportError: std::bad_cast` in a fresh interpreter;
  the fixed module and complete package import successfully.
- All **23 owning-model Python regressions passed**, including fresh-process
  import and FileHandler roundtrips with omitted, `None` and explicit enum logging
  arguments. The full Python suite then completed with **6,619 passed, 29 skipped,
  10 expected failures and 7 unexpected passes**, and no failures.
- New C++ regressions cover filtering mixed optional adducts, replacement in both
  directions, preserved IDs/scores/selections and independent processing metadata
  across copies and moves. Changed C++ lines were formatted and whitespace checks
  pass.

The native layout and score contract are unchanged. Native Windows/macOS builds
and wheel packaging remain subject to the new CI run; this local Python build
uses linked mode rather than the wheels' stable-ABI split mode.

## Clang/OpenMP follow-up (2026-10-06)

The next CI run at `0e8037c7066092a11c63d201de1b10dbc796730e` progressed past
the previous Clang crashes and exposed a structured-binding lambda capture in
the released-OMS compatibility reader. Clang 18.1.3 rejects this capture when
`-fopenmp` is enabled. The earlier isolated Clang checks did not enable OpenMP.
An explicit value capture of the observation ID removes that unsupported capture
without changing which empty observations are retained.

The original error reproduces locally with C++23 and `-fopenmp`; the fixed reader
compiles successfully. All **46 non-deleted C++ translation units changed by this
PR pass Clang 18.1.3 syntax compilation with OpenMP enabled**: 43 use their generated
build commands, and the external-consumer source and two benchmark sources use
the same configured consumer include paths and compiler flags. Deleted sources
are excluded. This covers the native codecs, migrated consumers, bindings,
changed class tests and TOPP sources.

The Release OMS reader was rebuilt and the core library relinked. All **21 focused
C++ suites and 23 owning-model Python regressions pass** against that library.
The full Python-suite result above belongs to the preceding CI-fix commit.

## Build and correctness

- Release `libOpenMS`: GCC 13.3, C++23, `-O3 -DNDEBUG -g1`, Arrow/Parquet 25.
- All **11 focused CTest suites passed**: the owning model, adapters, inference,
  feature/consensus workflows, native format, native inference codec, FileHandler
  integration, legacy model, legacy converter, existing FileHandler and OMS.
- All **15 Python regressions passed**, using the actual binding header in a
  standalone nanobind module linked to the Release library. The actual
  `bind_metadata.cpp` and `bind_format.cpp` translation units also compiled.
  The complete pyOpenMS package was not built.
- The five modified TOPP translation units also passed syntax compilation:
  AccurateMassSearch, IDFileConverter, IDMerger, MapAlignerIdentification and
  NucleicAcidSearchEngine. Their end-to-end pipelines were not run.
- Checks cover callback rollback, copied ownership, stable IDs and allocation
  counters, independent inference provenance, typed metadata/units/empty values,
  signed-zero and NaN bit preservation, custom modifications, optional adducts,
  missing tables, malformed descriptors, unsupported text and failed publication.
- The benchmark independently checks counts and score checksums. PyArrow inspection
  confirmed the physical file counts, inline score columns, multiple row groups,
  and absence of a stored per-query match count.

The complete library was built first. A repeatedly damaged Ninja dependency cache
in this workspace caused unnecessary rebuilds; the final changed translation units
and test executables were rebuilt and linked using their CMake-generated commands.
The recovery script and final logs are included with the implementation bundle.

## Release measurements

Machine: AMD EPYC 9V74 host, 8-CPU cgroup quota and 8 GiB memory limit.
No competing builds ran during these measurements. `OMP_NUM_THREADS=1`.
Each scan/load/filter ran in a fresh process. Linux `wait4` supplied that child's
peak RSS. Table times are measured operation times; RSS includes the entire process.
Files were recently written, so these are warm-cache measurements, not cold-storage
or remote-object-store measurements. They are individual runs, not statistical
confidence intervals.

The synthetic dataset has one candidate per query, two numeric score columns,
repeated molecular strings, parent evidence and two metadata fields. Compression
ratios and inference complexity must not be generalized to all biological data.

Each entry below is **seconds / peak MiB**. Filtering reads all candidates and
writes the retained half to a new dataset.

| Input | Score scan | Full payload scan | Filter + write | Owning load |
| --- | ---: | ---: | ---: | ---: |
| 100,000 PSMs, 1 run | 0.086 / 90.0 | 0.323 / 113.5 | 0.373 / 142.7 | 0.442 / 192.3 |
| 1,000,000 PSMs, 1 run | 0.947 / 105.0 | 3.185 / 128.9 | 4.141 / 163.2 | 4.915 / 1,023.2 |
| 1,000,000 PSMs, 1,000 runs | 1.947 / 73.8 | 5.152 / 74.9 | 8.464 / 82.1 | 6.695 / 895.5 |

A second score scan took 0.824 seconds for one million PSMs in one run and
1.943 seconds for the same count across 1,000 runs.

Constructing one million owned PSMs took 1.000 seconds for one run and
0.830 seconds for 1,000 runs. Writing these complete datasets took 3.074 and
6.440 seconds respectively. These write phases include validation, descriptor
collection, encoding and compression.

| Dataset | Total bytes | Parquet files | Query/match row groups per first run | Manifest bytes |
| --- | ---: | ---: | ---: | ---: |
| 100,000 PSMs, 1 run | 1,423,980 | 2 | 2 / 2 | 2,843 |
| 1,000,000 PSMs, 1 run | 14,007,168 | 2 | 16 / 16 | 2,847 |
| 1,000,000 PSMs, 1,000 runs | 34,533,772 | 2,000 | 1 / 1 | 2,738,772 |

## Real search-result fixture

`src/tests/topp/IDFilter_missed_cleavages_input.idXML` is an MS-GF+ search fixture
with 614 queries/matches, one run, one source and one imported protein result.
It is a small correctness-oriented example, not a representative large study.

| Operation | Seconds |
| --- | ---: |
| Parse idXML | 0.0912 |
| Import owning values | 0.00834 |
| Write native dataset | 0.0435 |
| Score scan, fresh process | 0.00422 |
| Full scan, fresh process | 0.0186 |
| Owning load, fresh process | 0.0300 |
| Filter half + write | 0.0284 |

The input occupied 1,890,671 bytes. The native output occupied 515,911 bytes in
nine Parquet files plus its manifest: three run tables and six inference tables.

## Pooled inference

A separate synthetic benchmark uses 10,000 distinct peptidoforms and 2,000
proteins, with consistent parent mappings across runs. It checks both membership
and assignment counts. Peak memory includes construction and inference.

| Inputs | Inference seconds | Process peak MiB |
| --- | ---: | ---: |
| 10,000 PSMs, 1 run | 0.0744 | 82.1 |
| 100,000 PSMs, 10 runs | 0.840 | 303.1 |

This measures the in-memory BasicProteinInference bridge for this graph. Ambiguous
peptide/protein graphs can have different costs. It is not an external-memory
inference algorithm.

## What dominates and what scales

Full payload decoding costs substantially more than reading numeric scores.
Owning loading adds object allocation and validation and has memory proportional
to the loaded data. Many short runs increase file-open, schema and row-group
work: in this example, filtering 1,000 small runs took longer than loading them.
These are phase measurements, not CPU-sampling evidence attributing time to an
individual function.

The 10x increase from 100,000 to one million PSMs increases scan/filter memory
modestly while owning load grows from 192 to 1,023 MiB. Sequential reading and
filtering do not build a PSM-sized lookup map. Strict uniqueness checks and random
lookup indexes have explicit additional memory costs.

A billion PSMs were **not tested**. The format permits processing them in bounded
batches or loading manageable runs. Its manifest/descriptor footprint, Arrow
buffers and largest individual payload remain additional memory costs. Global
protein inference and fully owning loading still need memory proportional to
their working data. The file format does not make those operations memory bounded.

## Migration

The old reference graph has been removed. Feature and consensus maps, OMS I/O,
filtering/FDR, RT alignment, accurate-mass search, NASE, IDFileConverter, IDMerger,
MapAlignerIdentification and their bindings now use the owning model. Existing
PeptideIdentification/ProteinIdentification workflows still use explicit adapters.

Release validation of this migration (5 October 2026):

- Rebuilt the library consumers and five changed TOPP executables with the generated
  Release compile/link commands. Final incremental builds include every changed C++
  implementation and regression test.
- All **21 focused C++ suites passed**, covering the owning model, adapters,
  inference, native I/O, released OMS compatibility, transactional failures,
  feature/consensus measurements and links, RNA digestion/conversion, filtering,
  pooled FDR, RT alignment, mzTab-M and existing FileHandler.
- All **20 Python regressions passed** in the standalone nanobind harness.
  The actual metadata, format, kernel, chemistry, processing and analysis binding
  translation units compiled. The full pyOpenMS package was not built.
- End-to-end invocations passed for AccurateMassSearch (native OMS annotation and
  mzTab-M), fresh NASE search, saved-digestion reuse, OMS IDMerger and OMS RT
  alignment. OMS results converted to idXML successfully. NASE fresh, reused and
  converted outputs agree on molecular identity, charge, target/decoy and score.
  This is a functional smoke check, not a pass of the entire TOPP golden-file suite.
- Native and OMS smoke writes/reads preserve a common content digest for **1,000
  PSMs across 10 runs**. No new million-PSM timing or memory claim is made.
- Removed production graph headers/classes and checked that no source or binding
  still references LegacyIdentificationData or its reference-update machinery.
  Changed C++ lines are formatted and `git diff --check` passes.

OMS schema 6 preserves typed quantitative metadata, units and string-list bytes;
feature hull comparisons canonicalize the internal cache and retain exact persisted
geometry. Empty optional native tables need no physical file. Explicit scoreless
sequence catalogs do not weaken the primary-score contract for search PSMs.
Compatibility converter APIs warn when older output formats cannot retain complete
score/provenance definitions; the adapter itself still defaults to strict export.

Other operating systems, the full Python package and billion-PSM execution remain
untested. The benchmark sections below describe earlier implementation commits and
are historical; their old OMS timings do not describe owning schema 6.

## Common primary score contract (2026-10-05)

Added checked dataset-wide primary score compatibility, an atomic score-selection
API, legacy import rejection, and native descriptor validation including streaming.
Release library and affected tests rebuilt. All 11 focused C++ suites pass; all 16
standalone Python binding tests pass (including the new contract API). This is a
focused binding harness, not a complete pyOpenMS package build.

Regression coverage includes incompatible score direction/provenance, different
local column orders, incomplete supplementary scores, atomic failure, replacement,
copy independence, mixed legacy input, streaming descriptor rejection, and inference
explicitly consuming a supplementary score while preserving a common primary score.
The existing IdXMLFile_whole fixture mixes MOWSE directions and is now tested as a
rejection case; its compatible subset retains roundtrip coverage.

## Shared Parquet tables and writer views (2026-10-05)

The per-run/per-result files were replaced with shared physical tables and bounded
row groups. Manifest slices and partition columns preserve ownership and direct
run access. The writer accesses payloads/scores by reference; readers reuse decoded
row groups across adjacent slices. The format is unreleased and no compatibility
reader for the discarded layout is retained.

Release build and all 11 focused C++ suites pass, including new shared-run and
shared-inference roundtrips/filtering tests. All 16 standalone Python binding tests
pass. The matched million-PSM comparison passes all content digests. For 1,000 runs,
write time falls from 16.86 s to 5.36 s, full load from 11.77 s to 6.33 s, disk size
from 79.65 MB to 15.48 MB, and file count from 9,001 to 10. Writing is now comparable
to existing PSM Parquet; full loading remains slower. Earlier timings above describe
the prior format/workloads. See `identification_data_legacy_comparison.md` and the
benchmark's `results/` directory for methods, limitations and raw observations.

## Inline protein-group membership (2026-10-05)

Group membership now uses an ordered Arrow list in the group row. The separate
member table, join keys, member ordinals and redundant member counts are removed.
Aliases without protein hits, optional qualified identities, duplicates, ordering,
empty groups, scores and typed arrays roundtrip unchanged. The whole group row is
checked against max_record_bytes during writing, owning loading and filtering.

The Release library and affected inference suite were rebuilt. All 11 focused C++
suites and 16 Python binding tests pass. A 1,000-PSM / 10-run roundtrip retains the
content digest; independent PyArrow inspection confirms the nested member schema
and nine output files. Previous million-PSM timings describe the shared-file layout
immediately before this membership simplification; no new speedup is claimed.

## Run-level inference input provenance (2026-10-05)

Inference inputs now retain only the run identifier/UUID, optional score definition
and selection description. The per-PSM input vector, membership-known flag, member
count and input_members table are removed from the C++ model, Python API and native
format. Match-to-parent assignments remain independent, including explicit empty
parent lists and retained references to removed PSMs. Their counter and input-run
validation remains checked. Run-level provenance alone does not reconstruct the
exact candidates used by an earlier calculation.

The Release library, focused tests, benchmark and standalone Python bindings were
rebuilt. All 11 focused C++ suites and all 16 Python tests pass; the actual metadata
and format binding translation units compile. Coverage includes pooled inference,
filtering, retained assignments, absent input runs, optional scores and native
roundtrips. A fresh 1,000-PSM / 10-run write/read preserves the content digest, and
PyArrow confirms six input columns including the partition ID and eight output
files. A pooled inference smoke check validates 10 input runs and 1,000 assignments.
Previous million-PSM timing and size measurements have not been rerun or relabelled.

## Protein/group results without PSM assignments (2026-10-05)

The base model and native format no longer contain inferred PSM-to-protein
assignments. Their C++ and Python APIs, serializer table, per-match reservation
checks, input-reference lookup and bridge bookkeeping have been removed. Original
search evidence remains on matches. Protein/group results, their scores and
run-level input provenance remain supported. A future PSM-to-peptide-to-protein
inference graph and its persistence require a separate design and may use a
different format.

The Release library and focused test targets rebuild successfully. All 11 focused
C++ suites and 16 Python binding tests pass; the real metadata/format binding
translation units compile. Updated coverage checks unchanged original candidates,
protein/group filtering, run-level provenance, live ID allocation counters,
missing-table failures and native roundtrips. The benchmark rebuilds and pooled
inference runs successfully on 1,000 PSMs from 10 runs. A fresh native write/read
preserves the content digest and confirms seven files, with only inputs, proteins
and groups declared for inference. Earlier sections describe previous layouts;
no new million-PSM timing or memory improvement is claimed.

## One ordered score schema and plain filenames (2026-10-05)

Configured runs now share the complete ordered score-definition vector and primary
column. The model checks supplementary definitions/provenance and ordering as well
as the primary selection. Primary values are required; supplementary values may be
null. Empty unconfigured runs remain valid construction placeholders and preserve
that state on roundtrip. Dataset-level getScoreDefinitions() exposes the checked
schema; ScoreId handles remain run-local. Atomic primary switching rejects schema
mismatches before changing any run.

Native inspection, loading, scanning and filtering validate the contract before
reading rows. The physical writer uses a single schema per table, rejects conflicting
layouts and writes plain names such as matches.parquet. Schema fingerprints,
variant counters and numeric filename suffixes are removed. No implicit union or
column remapping occurs. Supplemental legacy values remain typed match metadata.

The Release library and focused tests were rebuilt. All 11 focused C++ suites and
16 Python tests pass, including reordered/missing/extra supplementary columns,
provenance mismatches, required primary scores, nullable supplementary values,
atomic failures, transactional malformed-manifest rejection and shared row groups.
An additional regression verifies projected scans across empty unconfigured runs.
The standalone Python extension and real metadata/format binding translation units
compile. A rebuilt benchmark writes and reads 1,000 PSMs over 10 runs with matching
content digests, seven plain filenames and one shared matches table; pooled inference
also succeeds. No new million-PSM performance claim is made.


## Release benchmark after format simplification (2026-10-05)

Reran the unchanged four-format driver against implementation commit
`4cef855e6c6fc60cc8c90a29093f31ccc1c83502`. Release benchmark targets were up to date;
no implementation or benchmark code changed. All 80 fresh-process operations
passed common-content digest checks, covering 1,000 PSMs / 2 runs, 100,000 / 1,
1,000,000 / 1 and 1,000,000 / 1,000. Every native dataset has seven plain filenames.

For one million PSMs over 1,000 runs, native write/load measured
6.18 s / 6.77 s, versus
6.01 s / 4.85 s for existing PSM Parquet.
Native peak read-process RSS is 1381 MiB and disk
size is 14.92 MB. Loading remains slower than existing PSM Parquet;
the file-count simplification does not establish an additional speed improvement.
Loads are warm-cache medians of three measured reads after warmup; writes are
single samples excluding preparation/conversion. No concurrent builds ran.

The latest tables, conversion costs and limitations are in
`identification_data_legacy_comparison.md`; raw observations and per-file sizes
are in `tools/benchmarks/identification_data/results/shared-score-schema.json`.
Earlier sections retain their original implementation-specific results.


## Copy/allocation reductions and bounded Arrow threading (2026-10-05)

Release library, benchmarks and focused test targets rebuilt. All 11 focused C++
suites pass. The standalone nanobind harness passes 18 tests, including serial and
four-worker roundtrip/scan variants. Actual metadata/format binding translation
units compile. Formatting checks cover changed C++ lines. The full pyOpenMS package
and other operating systems were not built.

New regressions exercise moved/copy-preserved payloads, dense-score shape and
finite/required-value validation, arbitrary and duplicate record IDs, shared row
groups, serial/threaded interoperability, concurrent independent pools, unchanged
Arrow global worker capacity and transactional failure while decoding a damaged
column page in parallel. Existing inference, file-handler and OMS coverage passes.
The constructor/import boundaries retain every validation guarantee; no revision
machinery or extra physical tables were introduced.

Final benchmarks cover one million PSMs over one and 1,000 runs, existing PSM
Parquet and native worker counts 1/2/4/8. All 50 process results pass the common
content digest, as do 12 small worker-count checks and 10 intermediate-stage
measurements. Serial native load is 3.61 / 4.49 s;
four-worker load is 3.48 / 4.32 s. Serial writing is
3.91 / 4.73 s and four-worker writing is
3.58 / 4.37 s. Loads are warm-cache medians of three reads;
writes are single samples excluding generation/conversion. No concurrent builds
ran. The native layout still has seven files. Full results, ranges, memory, disk
sizes, source blob IDs and limitations are in the comparison report and benchmark
results directory. Threading stays opt-in via `Options.threads`; default is 1.
