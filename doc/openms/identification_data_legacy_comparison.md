# IdentificationData: Release comparison against legacy

## Optimized loading, validation and optional threading (2026-10-05)

Implemented and measured after the seven-file baseline below. The score contract,
seven-file layout and scientific record order are unchanged. The measured source
blob IDs and complete raw observations are recorded in
`tools/benchmarks/identification_data/results/optimized-threading.json`.

The implementation now:

- Moves decoded observation/match payloads and inference results into their owners,
  and moves temporary scan records into callback batches.
- Validates each fully loaded run once, reads dense score values directly and avoids
  constructing an optional-score vector for every match during validation.
- Checks uniqueness using compact ID vectors: a linear pass for sorted IDs, or
  sorting ID copies when IDs are not ordered. Every duplicate is still rejected.
- Binds projected columns once per batch and caches metadata registry IDs in the
  run-local descriptor dictionary. Registry IDs are not written to disk.
- Exposes `Options.threads` (default 1) in C++ and Python. Higher counts use a private
  CPU pool shared by the operation's tables, with parallel column decoding and
  buffered row-group encoding. The global Arrow pool is unchanged.

Both run-count cases contain one million PSMs and use the unchanged synthetic
content generator. Existing PSM Parquet was regenerated and measured alongside
native with 1, 2, 4 and 8 workers. All processes run sequentially; no builds ran
concurrently. Load values are medians of three fresh-process reads after warmup;
load order rotates. Writes are single samples. I/O timings exclude generation,
conversion and digest verification. Memory is process peak RSS in MiB; disk is MB.

| Runs | Format / CPU workers | Write (s) | Load (s) | Load range (s) | Peak read RAM (MiB) | Disk (MB) |
|---:|---|---:|---:|---:|---:|---:|
| 1 | Existing PSM Parquet (defaults) | 5.46 | 3.90 | 3.76–3.95 | 1622 | 10.79 |
| 1 | Native / 1 | 3.91 | 3.61 | 3.61–3.86 | 1233 | 17.12 |
| 1 | Native / 2 | 3.80 | 3.60 | 3.35–3.90 | 1227 | 17.12 |
| 1 | Native / 4 | 3.58 | 3.48 | 3.40–3.72 | 1243 | 17.12 |
| 1 | Native / 8 | 3.77 | 3.64 | 3.48–3.78 | 1312 | 17.12 |
| 1,000 | Existing PSM Parquet (defaults) | 5.90 | 3.94 | 3.90–4.13 | 1682 | 10.82 |
| 1,000 | Native / 1 | 4.73 | 4.49 | 4.32–4.61 | 1375 | 14.92 |
| 1,000 | Native / 2 | 4.44 | 4.55 | 4.30–4.57 | 1380 | 14.92 |
| 1,000 | Native / 4 | 4.37 | 4.32 | 4.15–4.59 | 1412 | 14.92 |
| 1,000 | Native / 8 | 4.33 | 4.13 | 4.11–4.50 | 1475 | 14.92 |

Relative to the immediately preceding native baseline, serial loading decreases
from 5.82 to 3.61 s for one run (38% less time), and from 6.77 to 4.49 s for
1,000 runs (34% less time). Serial writing measures 3.91 / 4.73 s, compared with
the previous 5.48 / 6.18 s. These are sequential before/after experiments, not
randomized paired trials.

Against the simultaneously rerun existing PSM Parquet API, native writing is
faster in these samples. Full loading is comparable for one run, but the serial
native loader still takes 14% longer for 1,000 runs. Eight workers reduce that
gap to about 5%, with higher peak memory than serial native loading.

Four workers measure 3.58 / 4.37 s writing and
3.48 / 4.32 s loading (one / 1,000 runs).
The extra threading gains are modest and not monotonic across worker counts;
several load ranges overlap. Removing repeated per-record work supplies most of
the observed loading improvement. Keep the default at one worker and expose
explicit control for callers to measure against their workloads. There is no
claim of linear scaling with cores. The threaded path keeps one decoded row group
per cached physical table, but encoding/decoding workers need additional buffers.

Conversion from the common legacy vectors adds 2.41 /
4.90 s for serial native in these runs. Whole-process write RSS and
CPU time include that conversion and input generation, so they are retained only
in the raw observations. Projected scan throughput, inference cost and cold-disk
service time are separate measurements.

An intermediate serial stage (moves, dense-score validation, one validation pass
and lookup caching, before compact ID validation) measured 4.84 / 4.81 s writing
and 4.02 / 4.33 s loading. Its raw records and exact patch against the preceding
commit are retained as `serial-optimization-stage.json` and
`serial-optimization-stage.patch` in the benchmark results directory. This measures
a bundle of changes and does not isolate the contribution of each optimization.

All 50 final benchmark subprocesses, 12 worker-count smoke subprocesses and 10
intermediate-stage subprocesses passed count/score/content digest checks. Every
final native dataset has seven files. Independent Parquet inspection verifies
query/match row counts and row-group counts. All 11 focused C++ suites and 18
standalone Python binding checks pass; the real metadata/format binding translation
units compile. This is not a complete pyOpenMS package or cross-platform build.
New coverage includes consuming imports, unordered/duplicate IDs, dense-score
validation, projected/shared row groups with four workers, concurrent operations,
serial/threaded cross-reading, unchanged global thread-pool capacity and a corrupt
encoded page that fails during parallel decoding without modifying the destination.

Reproduce the final matrix with:

```bash
python3 tools/benchmarks/identification_data/run_legacy_comparison.py \
  /absolute/path/IdentificationDataLegacyBenchmark /absolute/path/new-output \
  --formats parquet native --cases 1000000:1 1000000:1000 \
  --native-threads 1 2 4 8
```

The Release environment and workload details below still apply. This synthetic
million-PSM comparison does not establish billion-PSM feasibility or full workflow
performance. Earlier tables below are retained as historical measurements.

## Why full loading needs about 1.2 GiB per million PSMs

A diagnostic of the same Release implementation separates the owning model from
its compressed Parquet representation. On this Linux x86_64 build, the per-PSM
object sizes and occupied payloads in the one-run workload are:

| Component | Bytes per PSM | MiB for one million |
|---|---:|---:|
| Match object, including inline optional fields and vector/string headers | 376 | 358.6 |
| Observation/identification object | 120 | 114.4 |
| One parent-evidence entry | 160 | 152.6 |
| Metadata objects, occupied entries and separate string objects | 264 | 251.8 |
| Score storage, heap string capacity and spare query capacity | about 36 | 34.2 |
| **Accounted PSM storage, before the remaining overhead** | **about 956** | **911.5** |

This is an accounting lower bound, not a complete heap attribution. The remaining
roughly 321 MiB includes allocation rounding/headers, metadata spare capacity,
run/protein/inference data, libraries and Arrow/allocator buffers. Peak RSS was
1,232.8 MiB in this diagnostic. After loading returned, RSS was still 1,233.0 MiB;
explicitly trimming free glibc pages reduced it only to 1,217.3 MiB. Thus the large
footprint cannot be explained mainly by temporary decoding buffers or benchmark
verification. The diagnostic performs neither input conversion nor digest checking.
Small differences between sampled RSS and peak RSS reflect the different kernel
accounting interfaces.

Concrete design costs:

- `std::optional<AdductInfo>` occupies **120 bytes even when absent**: 114.4 MiB
  across these one million matches. This is already included in the match row
  above. An absent optional formula similarly occupies 40 bytes. Empty name,
  identifier and annotation containers also retain their inline headers.
- Parent evidence owns database/accession strings, two optional 64-bit coordinates
  and two string-valued flanking residues. That is 160 bytes before extra string
  storage, even for a peptide with single-character flanks.
- Every query owns a match vector and every match owns a score vector. Even this
  one-hit, one-score workload therefore makes many small allocations.
- The legacy adapter copies metadata wholesale: spectrum references and
  target/decoy labels coexist with their dedicated fields. Repeated metadata
  values such as the engine name also become independent C++ values.
- The generator has only 10,000 distinct peptide representations and 200 protein
  accessions in its one-run case. Parquet compresses repeated values well; the
  owning model constructs individual values for every PSM.

A full-projection streaming scan of the same one-run dataset, retaining no callback
batches, peaked at **192.2 MiB**, with one million queries and matches counted.
This scan does not materialize the whole dataset or load protein inference; it
uses the default bounded buffering and disables optional set-based global ID
uniqueness checks. It is not a billion-PSM memory guarantee. The 1,000-run full-load
diagnostic peaked at 1,370.7 MiB. That workload also contains 200 protein hits per
run rather than 200 total, so its extra memory is not solely run overhead.

The next memory-oriented changes should move uncommon molecular fields into an
optional owned extension, compact parent evidence, avoid duplicate canonical
metadata and reduce small allocations for common one-hit/one-score records. These
are design recommendations, not changes included in the I/O optimization. They
can retain simple value ownership without introducing a graph of mutable
cross-references. Run-owned dictionaries are another option for repeated strings,
provided their lifetime and editing rules remain explicit.

Raw phase measurements, object counts and the diagnostic source are archived in
`tools/benchmarks/identification_data/results/memory-footprint.json`.

## Seven-file baseline before these optimizations (2026-10-05)

Rerun on implementation commit `4cef855e6c6fc60cc8c90a29093f31ccc1c83502` after the
ordered score-schema contract and format simplifications. All four formats were
regenerated and measured again using the unchanged comparison driver.

**Native writing remains comparable to existing PSM Parquet; full owning loading
remains slower.** The current format produces seven files for both run counts.

Each case contains 1,000,000 PSMs. Load time and peak load-process RSS are medians
of three fresh-process reads after one warmup, with rotated format order. Writes
are single measurements, excluding generation, conversion and digest checks.
Disk sizes are decimal MB; resident memory is MiB. These are warm-cache API
measurements, including decoding and object construction.

| Runs | Format | Write (s) | Load (s) | Load range (s) | Peak load RAM (MiB) | Disk (MB) | Files |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | idXML | 3.98 | 11.19 | 10.92–11.34 | 857 | 584.44 | 1 |
| 1 | Existing PSM Parquet | 5.99 | 3.66 | 3.66–3.87 | 1623 | 10.79 | 4 |
| 1 | OMS / LegacyIdentificationData | 12.90 | 8.59 | 8.07–8.93 | 872 | 309.86 | 1 |
| 1 | Current native format | 5.48 | 5.82 | 5.80–6.05 | 1289 | 17.12 | 7 |
| 1,000 | idXML | 15.31 | 12.61 | 12.04–13.22 | 921 | 618.75 | 1 |
| 1,000 | Existing PSM Parquet | 6.01 | 4.85 | 4.46–4.90 | 1682 | 10.82 | 4 |
| 1,000 | OMS / LegacyIdentificationData | 14.71 | 16.85 | 16.51–17.49 | 957 | 357.54 | 1 |
| 1,000 | Current native format | 6.18 | 6.77 | 6.68–6.78 | 1381 | 14.92 | 7 |

Native full loading takes 1.59× as long as existing PSM Parquet for one run and
1.40× for 1,000 runs, with respectively 21% and 18% less peak read-process RAM.
The latest simplification reduces the prior shared-table layout
from ten files to seven; these measurements do not establish an additional speedup
over that implementation. Single-write samples cannot establish small performance
differences. The earlier measurements below remain historical observations.

Conversion from the common peptide/protein vectors costs another
3.13 s / 5.99 s for native and
5.55 s / 29.70 s for OMS (one / 1,000 runs).
Include this cost when evaluating conversion pipelines. Process peak RSS for
writing includes generation and conversion, so it is retained in the raw records
but omitted from the table.

All 80 subprocesses across the 1,000-, 100,000- and million-PSM cases completed
successfully. Counts, primary-score sums and the common-field content digest agree
across formats, conversion and every roundtrip. These checks do not compare every
metadata field. Native output has exactly `manifest.json`, `queries.parquet`,
`matches.parquet`, `parents.parquet`, `inputs.parquet`, `proteins.parquet` and
`groups.parquet` in every case. Raw observations, environment, file sizes and
summary statistics, including independently inspected Parquet row counts and
row-group counts, are archived in
`tools/benchmarks/identification_data/results/shared-score-schema.json`.

This rerun uses the same synthetic modified-peptide workload, default codecs,
Release build and environment described below, with no concurrent builds.
It measures full loading rather than projected scans, and does not establish
cold-disk throughput, inference performance or billion-PSM scalability.

## Optimization audit before implementation

The following audit motivated the implementation reported at the top of this
document. It describes the previous code; a pipeline overlapping row preparation
with encoding remains future work.

| Priority | Previous work | Proposed change at that stage |
|---:|---|---|
| 1 | `readRun` passes decoded observations and match payloads through const references; `importIdentification` and `importMatch` copy them into owned records. Scan batching also copies temporary records. | Add an internal consuming import/batch path and move owned strings, evidence, annotations and metadata. Keep copying APIs for callers that retain their inputs. |
| 2 | `addRun` validates each loaded run, and `load` then calls dataset `validate`, which validates those runs again. Every `Run::validate` calls `match.getScores()`, allocating an optional-score vector for every match. | Validate dense score storage directly and consolidate the load-time validation passes while preserving every required check and transactional publication. |
| 3 | `readMatch` repeatedly resolves columns by name and constructs score-column names. Metadata code repeatedly allocates name lists and resolves registry names. | Bind typed columns once per batch, precompute score indices and cache metadata descriptor-to-registry mappings. Preallocate bounded batch storage. |
| 4 | The Parquet reader uses Arrow's default column-threading setting, which is false in the installed Arrow 25 headers. | Measure `set_use_threads(true)` with an explicitly bounded Arrow worker pool. This parallelizes decoding, not the subsequent construction of owned OpenMS records. |
| 5 | The writer synchronously builds rows and calls `WriteTable`. | Evaluate buffered row-group writing with parallel column encoding/compression; then consider a bounded pipeline that overlaps row preparation with encoding. |

Relevant source paths are `src/openms/source/FORMAT/IdentificationDataFile.cpp`,
`IdentificationDataFileSupport.cpp` in the same directory, and
`src/openms/source/METADATA/ID/IdentificationData.cpp`. `addInferenceResult` also
copies a freshly loaded result, so a consuming overload can avoid that protein/group
copy. The owning writer already uses `QueryView`/`MatchView` references; moving data
out of the caller's records is not appropriate for a const serialization API.
Likewise, applying `std::move` to the loader's existing const-reference parameters
alone does not enable the required moves.

Arrow's writer threading option applies to **buffered row-group mode**. The current
`WriteTable` path must be adapted before relying on this setting for column
parallelism. See the installed Arrow 25 `parquet/properties.h` and
`parquet/arrow/writer.h`, and the official
[Arrow C++ format API](https://arrow.apache.org/docs/cpp/api/formats.html).
Arrow also documents a nested-executor deadlock risk when parallel file writing
and parallel column writing share the same executor.

For further parallel loading, prefer bounded work on physical row groups over
launching a task for every run: many small runs share a row group, so independent
run readers can decode it repeatedly. The current `ReadPool` contains mutable
shared cached tables and cannot simply be used concurrently by those tasks.
Owned record construction must preserve ordering and use disjoint output storage;
metadata-name registry operations also contain OpenMP critical sections. Keep
queues and decoded batches bounded and measure RSS as well as throughput.

Recommended experiment order: remove payload copies and optional-score temporaries,
consolidate validation, then compare Arrow decoding/encoding with 1, 2, 4 and 8
workers on both run-count cases. Profile phases separately to distinguish Parquet
work from owned-record construction, and preserve the seven-file layout and checked
score contract throughout.

## Historical baseline before consolidation

Measured on implementation commit `d1fcc143eaf64431da9b25bda9ba078dfd79118e` with the accompanying benchmark additions.

Read time and peak resident RAM are medians of three fresh-process loads after one warmup. Write time is a single measurement, excludes data generation and conversion, and should be treated as indicative. Disk sizes are decimal MB; resident memory is MiB.

### 1,000,000 PSMs / 1 run

| Format | Read (s) | Read range (s) | Write (s) | Peak read RAM (MiB) | Disk (MB) | Files |
|---|---:|---:|---:|---:|---:|---:|
| idXML | 11.07 | 10.94–11.34 | 4.55 | 857 | 584.44 | 1 |
| Existing PSM Parquet | 4.00 | 3.79–4.03 | 5.67 | 1622 | 10.79 | 4 |
| OMS / LegacyIdentificationData | 8.65 | 8.20–8.89 | 13.79 | 872 | 309.86 | 1 |
| New native format | 5.75 | 5.74–6.33 | 5.88 | 1278 | 17.11 | 10 |

### 1,000,000 PSMs / 1,000 runs

| Format | Read (s) | Read range (s) | Write (s) | Peak read RAM (MiB) | Disk (MB) | Files |
|---|---:|---:|---:|---:|---:|---:|
| idXML | 12.39 | 11.74–12.86 | 15.34 | 921 | 618.75 | 1 |
| Existing PSM Parquet | 4.10 | 3.93–4.78 | 6.37 | 1683 | 10.82 | 4 |
| OMS / LegacyIdentificationData | 17.37 | 17.23–17.55 | 15.11 | 958 | 357.54 | 1 |
| New native format | 11.77 | 11.64–11.78 | 16.86 | 1285 | 79.65 | 9,001 |

## Baseline interpretation

- For one run, native loading takes 44% longer than existing PSM Parquet, while using 21% less peak RAM. It is 1.92× faster than idXML and 1.50× faster than OMS.
- Across 1,000 runs, native loading takes 2.87× as long as existing PSM Parquet. Native write time rises from 5.88 s to 16.86 s.
- The native legacy import produces three run tables and six inference tables per run, plus a manifest: 9,001 files for 1,000 runs. The importer preserves one protein result per input analysis run. Existing PSM Parquet produces four files. This makes file/table overhead an optimization target; these measurements do not isolate its exact share of CPU or I/O time.
- Native disk size grows from 17.11 MB to 79.65 MB across 1,000 runs. Tiny tables and repeated run/protein metadata affect compression and file overhead. A producer without protein results would produce fewer tables.
- These are warm-cache API load measurements, not raw disk speed. They include decoding and construction of the respective native object models.
- Streaming score scans are a different operation and are deliberately excluded from this full-load comparison.

## Baseline conversion overhead

Converting the common peptide/protein vectors before writing costs 2.91 s (one run) / 6.34 s (1,000 runs) for the new model, and 6.98 s / 27.47 s for LegacyIdentificationData. Add these to write time when evaluating a conversion pipeline. Generation itself is separate. Write-process peak RAM is not presented because it includes generation and simultaneous conversion inputs.

## Workload and validation

All formats receive the same PSM content: 10,000 distinct sequences, 20% modified peptides (Oxidation(M) and Carbamidomethyl(C)), one candidate per query, primary score, PEP metadata, charge, RT, m/z, spectrum reference, target annotation, engine metadata, and protein evidence. Each run carries corresponding protein hits and search metadata. Splitting into more runs preserves global PSM content but naturally adds run and protein catalogs. This synthetic workload is richer than the earlier two-sequence benchmark; its timing and size are not directly comparable to that earlier 14 MB result.

Every conversion and roundtrip checks PSM count, primary score sum, and an order-independent digest of modified sequence, spectrum reference, charge, score, RT and m/z. All checks passed, including the 1,000-PSM smoke and 100,000-PSM intermediate cases. The digest checks common key content, not every metadata field. Verification runs outside timed I/O phases.

Both Parquet formats use their default ZSTD settings. Existing PSM Parquet produced one PSM row group for one million PSMs; native matches produced 16. This compares current API defaults, not identical compression/layout tuning. The in-memory types also differ: peptide/protein vectors for idXML and existing Parquet, a reference-based graph for OMS, and the owning model for native.

## Environment and reproduction

Linux x86_64, AMD EPYC 9V74 host, 8 CPU quota, 8 GiB memory limit; GCC 13.3, Release (`-O3 -DNDEBUG -g1`), Arrow/Parquet 25. `OMP_NUM_THREADS=1`; no concurrent builds. Load order rotates between repetitions. Fresh child-process peak RSS is collected with Linux `wait4`; filesystem caches remain warm.

Build the standalone benchmark using its CMakeLists.txt against a Release OpenMS build, then run:

```bash
python3 run_legacy_comparison.py /absolute/path/IdentificationDataLegacyBenchmark /absolute/path/new-output-directory
```

Set `LD_LIBRARY_PATH` and `OPENMS_DATA_PATH` for your build before invoking the driver. The archive includes source, the measurement driver and raw observations. The driver regenerates every dataset and checks each roundtrip. Timings are machine- and workload-specific; billion-PSM scalability and workflow/inference speed are not established by this benchmark.

## Historical native writer phase profile

Temporary instrumented Release build, one run per case, same one-million-PSM
workload. These are wall-clock phase measurements, not a statistical CPU sampling
profile. Instrumentation was kept outside the production code.

| Phase | 1 run (s) | 1,000 runs (s) |
|---|---:|---:|
| Full store | 5.738 | 18.345 |
| Model validation | 0.825 | 0.290 |
| Metadata dictionary collection | 0.627 | 0.706 |
| Query/match row loop, including automatic flushes | 3.631 | 2.928 |
| Query/match final flush and close | 0.032 | 5.218 |
| Parent tables | 0.001 | 1.339 |
| Inference serialization and descriptor insertion | 0.004 | 6.355 |

Nested timings (already included above; do not add to the phase totals):

| Operation across all tables | 1 run (s) | 1,000 runs (s) |
|---|---:|---:|
| Parquet WriteTable calls | 1.328 | 7.444 |
| Parquet Close calls | 0.004 | 1.996 |
| Arrow batch ValidateFull | 0.081 | 0.261 |

Unlisted time includes writer construction, descriptors/manifest, provenance checks,
and directory publication. Parquet WriteTable combines encoding, compression and
output; the instrumentation does not separate disk service time from CPU work.
With 1,000 small runs, query/match buffers flush on close rather than in the row loop,
so interpreting the loop alone would be misleading.

The many-run penalty is concentrated in writing/finalizing thousands of small
tables. The existing writer aggregates runs into four tables. Native store also
performs full model validation, a metadata dictionary traversal, temporary owning
QueryRecord/MatchRecord copies, score-vector materialization, metadata encoding and
record-size estimation. For one large run, row preparation and validation are
substantial; the profile does not isolate the individual cost of those allocations.
Consolidating tables across runs should be prioritized for the many-run case;
reference-based serialization views and avoiding repeated metadata/score work are
separate targets for the per-row path. Preserve validation guarantees when optimizing.

## Shared-file optimization (2026-10-05)

At this stage, the writer shared physical files and row groups across compatible runs/results, wrote through references instead of temporary owning records, and read slices through a bounded file/row-group cache. Complete supplementary score layouts could differ and were grouped into compatible match files. The common primary PSM score contract was checked; the stricter current contract described above supersedes this arrangement.

Same one-million-PSM workload and settings as above; new and existing Parquet were rerun without concurrent builds. Load values are medians of three fresh-process measurements after warmup. Write values are single measurements, excluding conversion.

| Runs | Format | Write (s) | Load (s) | Load range (s) | Peak load RAM (MiB) | Disk (MB) | Files |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | Existing PSM Parquet | 5.41 | 3.83 | 3.70–3.90 | 1623 | 10.79 | 4 |
| 1 | Optimized native | 5.35 | 6.12 | 5.88–6.15 | 1280 | 17.13 | 10 |
| 1,000 | Existing PSM Parquet | 5.60 | 3.83 | 3.83–3.88 | 1683 | 10.82 | 4 |
| 1,000 | Optimized native | 5.36 | 6.33 | 6.32–6.45 | 1383 | 15.48 | 10 |

For 1,000 runs, native writing improves from 16.86 s to 5.36 s (3.15×), full loading from 11.77 s to 6.33 s (1.86×), and disk size from 79.65 MB to 15.48 MB. File count falls from 9,001 to 10: nine Parquet tables plus the manifest. Writing is now comparable to existing PSM Parquet; the single-write measurements do not establish a statistically significant advantage. The single-run full loader has no demonstrated speedup. Peak memory for the many-run owning load increases modestly because decoded row groups are reused.

All common-content digests agree before and after conversion/roundtrip. All 11 focused C++ suites and 16 Python binding tests pass. Added regressions exercise shared row groups crossing run boundaries, differing supplementary schemas, single-run selection, per-run metadata dictionaries, filtering, and multiple inference results in shared files. No production timing instrumentation remains.

Raw observations are in `tools/benchmarks/identification_data/results/`: these record phase timings, process RSS, count/score/content digests and dataset dimensions. Absolute temporary output paths have been shortened to basenames; measured values are unchanged.

## Subsequent format simplification

Protein group members are now stored directly in each group row as a typed ordered
list. The separate group_members table and redundant member counts/ordinals are
removed. This initially reduced the benchmark layout above to nine files.

Exact PSM input membership has also been removed from the owning model and format.
Inference inputs retain the contributing run, optional score definition and selection
description. Inferred match-to-parent assignments have subsequently also been
removed from the base model and format. Protein/group results and original search
evidence remain supported; a future PSM-to-peptide-to-protein inference graph can
use a separate format. The current layout has **seven files** (six Parquet tables
plus the manifest). The timing and byte-size figures above predate these
simplifications and have not been relabelled as new measurements.


The current format also requires one ordered score schema across configured runs,
with one shared primary column and nullable supplementary values. It writes plain
filenames such as matches.parquet and rejects incompatible layouts instead of
splitting them into multiple files. The standard layout remains seven files; the
historical benchmark measurements above have not been relabelled.
