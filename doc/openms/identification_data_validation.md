# IdentificationData implementation validation

Measured on 5 October 2026 against OpenMS base commit `fe7b9cc` (Linux x86-64).
The owning model, native format, adapters and Python bindings are implemented.

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

## Migration boundary

The former reference-based implementation remains explicitly named
`LegacyIdentificationData` for consumers not yet migrated. The new public
`IdentificationData` is the owning implementation. No compatibility reader for
unreleased experimental native formats exists, and there are no revision tokens
or freshness gates.

Feature/consensus association adapters are implemented. Complete quantitative-map
persistence and converting every high-level workflow to the new API are separate
follow-up migrations. Existing idXML/OMS and PSM Parquet paths remain available;
strict conversion reports unsupported information instead of silently dropping it.
Validation here covers Linux; other platform builds were not run.

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
