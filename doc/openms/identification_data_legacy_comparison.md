# IdentificationData: Release comparison against legacy

Measured 2026-10-05 on implementation commit `d1fcc143eaf64431da9b25bda9ba078dfd79118e` with the accompanying benchmark additions.

**The new native full loader is slower than existing PSM Parquet.** The initial layout incurred substantial many-run overhead. The shared-file optimization at the end of this report removes most of that penalty; full loading remains slower than existing PSM Parquet.

## Baseline results before consolidation

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

## Interpretation

- For one run, native loading takes 44% longer than existing PSM Parquet, while using 21% less peak RAM. It is 1.92× faster than idXML and 1.50× faster than OMS.
- Across 1,000 runs, native loading takes 2.87× as long as existing PSM Parquet. Native write time rises from 5.88 s to 16.86 s.
- The native legacy import produces three run tables and six inference tables per run, plus a manifest: 9,001 files for 1,000 runs. The importer preserves one protein result per input analysis run. Existing PSM Parquet produces four files. This makes file/table overhead an optimization target; these measurements do not isolate its exact share of CPU or I/O time.
- Native disk size grows from 17.11 MB to 79.65 MB across 1,000 runs. Tiny tables and repeated run/protein metadata affect compression and file overhead. A producer without protein results would produce fewer tables.
- These are warm-cache API load measurements, not raw disk speed. They include decoding and construction of the respective native object models.
- Streaming score scans are a different operation and are deliberately excluded from this full-load comparison.

## Conversion overhead

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

## Native writer phase profile

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

The writer now shares physical files and row groups across compatible runs/results, writes through references instead of temporary owning records, and reads slices through a bounded file/row-group cache. Complete supplementary score layouts may differ and are grouped into compatible match files. The common primary PSM score contract is checked.

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
