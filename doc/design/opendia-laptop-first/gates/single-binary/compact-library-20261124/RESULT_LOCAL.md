# Compact library: local validation, 2026-11-24

## Status and scope

Local architecture, micro-measurement and raw01 end-to-end identity evidence. **Direct-Parquet attempt2 and the final release-binary cached-PQP gate both passed all three baseline byte comparisons.** The direct-Parquet arm used the earlier integrated binary identified below. The automatic-cache arm used the exact final committed binary. No release speedup, full-library capacity result or container validation is claimed here.

## Architecture

Source: `8c0a89f053284c9aebb620d36e351b14a642d15c` (`src/topp/OpenDIA/README.md` and associated implementation). Parquet and valid automatic caches stream row groups into a compact analysis table. Search fragments retain m/z, intensity and numeric transition ID. Precursor IDs and protein sets are interned, with folds, pair hashes, decoy partners and sibling groups computed once. A separate scoring library streams into the twelve lowest transition IDs per precursor instead of remaining as a second materialized library. Indexed coordinates serve calibration and trace collectors. Search applies calibration and RT corrections on the fly without a library copy per pass. RT refinement retains training feature rows and predicts the rest individually. The first uncached PQP read still materializes the source library to create the cache.

## Completed local tests

The retained `evidence/ctest_19.log` has **19 actual `Test Passed.` records**, not merely 19 scheduled tests. `ctest_19_symlink.log` independently also has 19 passing records. These are completed local CTest records, not full raw01 acceptance.

The tests are OpenDIAWorkflow, Search, Traces, Calibration, Stage2, PeakPicker, Interference, Scorer, GBM, FragmentGrid, MzIndex, RTRefine and LibraryParquet (`_test` suffix), plus OpenDIA_no_openswath, OpenDIA_fails_closed_on_tiny_library, OpenDIA_rejects_duplicate_run_names, TOPPWRITEINI_OpenDIA, TOPPWRITEINI_OpenDIA_SectionName and TOPPWRITECTD_OpenDIA. Coverage includes compact ownership/moves and independent materialized-library comparisons, scoring selection, indexed collector equivalence, training-only RT refinement and cache fallback/symlink protection. Exact invocations, executable paths and per-test outputs are in the retained test logs. `build_ctest.log` preserves build/configure context but its start-only summary is not the passing-test evidence.

## Heap measurements

Inputs: 603,421 precursors, 3,620,526 search fragments and 7,113,750 scoring fragments. All memory values are MiB (2^20 bytes), including old harness fields labeled MB. Values below preserve the precision printed by the harnesses.

| Allocated-heap measurement | Old | Compact |
|---|---:|---:|
| Standalone search Library | 724.379 | 724.379 |
| Standalone scoring Library | 1266.16 | 1266.16 |
| AnalysisLibrary incremental total | 1344.23 | 372.519 |
| sizeof Precursor, bytes | 424 | 328 |
| Precursor vector copy | 347.514 | 208.397 |
| Index copy | 42.2674 | 42.2172 |
| Analysis traml_ids copy | 36.8299 | removed |
| Compact fragments copy | n/a | 82.8711 |
| Fragment offsets copy | n/a | 4.60376 |
| Protein sets copy | n/a | 20.0644 |
| Calibration pairs copy | n/a | 2.29602 |
| Sibling group sizes copy | n/a | 1.0184 |

Observed analysis-table heap difference: 971.711 MiB, approximately 72.29% lower. Heap is `mallinfo2().uordblks + mallinfo2().hblkhd` minus a pre-load process baseline. Component copies measure incremental allocation, not exact owned capacity or additive accounting. Vector copies allocate size rather than original capacity, and compact precursor string views borrow original map keys. Standalone Library costs are unchanged. The compact heap harness reads/frees standalone libraries before the analysis measurement, so its retained RSS is not the clean-process comparison.

## Three sequential paired clean-process analysis probes

| Trial | Old load s | Compact load s | Old retained MiB | Compact retained MiB | Old peak MiB | Compact peak MiB |
|---|---:|---:|---:|---:|---:|---:|
| 1 | 41.7814 | 34.2904 | 1987 | 466.352 | 3189 | 943.383 |
| 2 | 40.8977 | 17.0210 | 1938 | 466.363 | 3211 | 942.914 |
| 3 | 29.4294 | 17.1778 | 1988 | 466.438 | 3190 | 943.309 |

Median retained RSS: 1987 versus 466.363 MiB, approximately 76.53% lower. Median peak: 3190 versus 943.309 MiB, approximately 70.43% lower. Old retained RSS uses integer-truncated SysInfo MiB, compact uses `/proc/self/statm` pages. Old internal peak is sampled before an additional search-Library copy, while external `/usr/bin/time -v` covers that copy too. Both forms are preserved in original logs. Compact clean-process allocated heap is 372.634 to 372.645 MiB, slightly above the heap harness's 372.519 MiB.

Median observed load is 40.8977 versus 17.1778 seconds. **These are concurrent micro timings, not isolated speedup measurements.** A simultaneous 16-thread raw01 identity run and other shared-machine activity can perturb CPU and memory bandwidth. An initial old probe measured 27.2194 seconds. Trials ran old then compact sequentially, without randomized order or an isolated cache-state protocol. `run-commands.log` records commands, timestamps (16:38:14 to 16:42:22 UTC) and competing processes. Dynamic dependencies were the then-current build libraries for both arms. The old implementation was statically embedded in its preserved executable.

## Provenance and reproduction

Final source is the commit above. The **assigned final binary** SHA256 is `83c408e5ad42e4865283a0f8586f18b223139f0f9a1b7d5c9f3ea491c9235eb7`, supplied by the coordinating root. This is not the executable used by the micro probes. Probe linkage used an earlier integrated library, with exact recorded identities:

| Artifact | SHA256 |
|---|---|
| Linked compact `libOpenDIA_impl.a` | `d651ace1b1ab16b8fa722e3b4bb000a47bd9563d58db6545c3fefc7203387bbe` |
| Old heap executable | `7a6e1b0c049782f91be1f796e694026ba86e01bb1c2d4278c80a3b3d46a3d143` |
| Old probe executable | `79cf5d9cd8267138ed0f58337a91930e9563311c37a902af182a74049f8ea114` |
| Compact heap executable | `3be1d668a148b95fed346cddf7f8b28a1633cdd1be64fa35b9a0031c38133f29` |
| Compact probe executable | `af23acb3df0df18fd17fdf893108d3e92abc7926c0f4ccab0e8a4c75c10b9660` |
| Earlier `OpenDIA.compact` in `binary.sha256` | `b0a6c5eba5a3b41478808f266485e9c061f217305a3885aab1e55d68e625510a` |
| Search input Parquet | `77adfb9bd76054856057a5109598ffd76d990cc7b8761510544b7f40fc9a040e` |
| Scoring input Parquet | `36429944a81c08e7adb304613e2a5b9c92663f3efc1e2698a38df1b9031bc202` |

These hashes come from original `RESULTS.md`, `initial-sha256.txt`, `binary.sha256` and `run-commands.log`, not an assumption that probes equal the final release binary. The old embedded implementation has an executable hash but no independently established source commit here. Test logs identify their build paths, not a hashed test executable manifest.

`evidence/` copies small stable files from `/home/singjc/.jcode/scratch/compact-20261124/`: original summaries, raw heap/probe logs, harness C++ sources, compile/link commands (`build_compact.sh`, `build.log`, `old/build_*.sh`), measurement commands, test records and planned identity/container command scripts. Scripts contain original absolute local paths and are provenance, not portable launchers or evidence that every planned command completed. No executables, objects or output TSVs are included.

Original `final-sha256.txt` is preserved verbatim as a historical manifest. It references artifacts intentionally not copied, including binaries, objects and a historical attempt2 log hash. That historical hash is not a final attempt2 result. The now-completed attempt2 log and timing file have been added separately. The original summaries' statements about scratch-only work and preserved executables describe the scratch session, not this documentation directory.

## Raw01 seed0 whole-workflow gate

`parquet_attempt1.log` and `parquet_attempt1.time.txt` retain the initial failure: missing `/home/singjc/.jcode/scratch/parquet/astral/lib_final/search.opendia.parquet`, exit status 1 after 5.05 seconds. This is an input-path failure, not a completed analysis or identity comparison.

Direct-Parquet attempt2 completed at 16:55:26 UTC with executable exit 0 and an empty `diff -rq` comparison against the retained `3b1d90ed67` baseline. All three output hashes match exactly (`baseline_outputs.sha256` and `direct_parquet_outputs.sha256`): the combined precursor TSV, run precursor TSV and calibration JSON. The raw01 SHA256 is `65d7c288949ac22290e93b6a4ad69d6eae5f944f8c0cbdc93dba75729f97a125`, also the frozen Codon raw01 input hash. The baseline executable and local PQP hashes are in `local_identity_inputs.sha256`.

| Whole-run GNU time | Baseline PQP | Earlier compact, direct Parquet | Final compact, cached PQP |
|---|---:|---:|---:|
| Wall seconds | 1136.26 | 1151.88 | 4618, includes system suspend |
| Maximum RSS, KiB | 17,963,204 | 16,307,904 | 16,324,332 |
| Actual exit | 0 | 0 | 0 |
| Exact baseline output files | reference | 3/3 | 3/3 |

These are single observed runs on the same busy desktop, not an isolated timing experiment. Baseline and direct-Parquet timing/cache conditions differ, including concurrent micro probes during the compact run. The direct-Parquet run did not become faster in this observation. Its peak RSS was lower, but this is not a full-library capacity result. Calibration anchors, RT fit and scorer diagnostics also agree with baseline. No FDR or quantification claim follows from byte identity.

The final committed binary cached-PQP gate started sequentially at 16:55:31 UTC after attempt2 ended, using the same raw, seed0, sixteen threads and retained baseline. It completed at 18:12:29 UTC with actual exit 0 and exact baseline hashes for all three files (`release_cached_outputs.sha256`). `release_cached_pqp.log` explicitly records both automatic Parquet cache hits, and before/after file size/mtime records show the PQPs and caches unchanged. `release_binary.sha256` identifies the exact release executable.

The machine suspended at 17:08:24 UTC and returned at 18:09:11 UTC (`release_suspend_provenance.log`). The original GNU wall time includes that interval and is retained without subtraction or replacement. It must not be used as a runtime comparison. Byte identity remains a valid result of the unchanged, completed attempt. Source release gates are passed. Container builds and the full-library Codon experiment are separate, subsequent gates.
