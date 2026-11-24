# Compact library heap/probe comparison

Date: 2026-11-24. Evidence is local scratch only. No repository source edits or commits.
Inputs: 603,421 precursors, 3,620,526 search fragments, 7,113,750 scoring fragments.
All sizes below are MiB (2^20 bytes), including old harness fields labeled MB.

## Heap

| Measurement | Old | Compact |
|---|---:|---:|
| Standalone search Library | 724.379 | 724.379 |
| Standalone scoring Library | 1266.16 | 1266.16 |
| AnalysisLibrary total incremental allocated heap | 1344.23 | 372.519 |
| sizeof Precursor, bytes | 424 | 328 |
| Precursor vector copy heap | 347.514 | 208.397 |
| Index copy heap | 42.2674 | 42.2172 |
| Analysis traml_ids copy heap | 36.8299 | removed |
| Compact fragments copy heap | n/a | 82.8711 |
| Compact fragment offsets copy heap | n/a | 4.60376 |
| Compact protein sets copy heap | n/a | 20.0644 |
| Compact calibration pairs copy heap | n/a | 2.29602 |
| Compact sibling group sizes copy heap | n/a | 1.0184 |

Heap total reduction: approximately 72.29%, 971.711 MiB. Components are copy-allocation measurements, not exact owned-capacity accounting. See README.txt for caveats.

## Clean process analysis probes, sequential paired trials

| Trial | Old load seconds | Compact load seconds | Old retained RSS | Compact retained RSS | Old peak RSS | Compact peak RSS |
|---|---:|---:|---:|---:|---:|---:|
| 1 | 41.7814 | 34.2904 | 1987 | 466.352 | 3189 | 943.383 |
| 2 | 40.8977 | 17.0210 | 1938 | 466.363 | 3211 | 942.914 |
| 3 | 29.4294 | 17.1778 | 1988 | 466.438 | 3190 | 943.309 |

Median observed load: old 40.8977 s, compact 17.1778 s. This is NOT an isolated speedup claim. Label: concurrent micro timings. Identity raw01 run with 16 threads and other shared-machine load may perturb timings. Initial old probe before paired trials measured 27.2194 s, illustrating timing variability.
Median retained RSS: old 1987 MiB, compact 466.363 MiB, approximately 76.53% lower. Peak approximately 3190 vs 943.309 MiB, 70.43% lower. Old internal peak sampled before the extra search-Library copy, with external time peak consistent. Compact clean-process allocated heap 372.634 to 372.645 MiB, slightly above post-read/freed-process 372.519 MiB.

## Evidence

- old/ preserves original executables, sources, objects and build scripts.
- initial-sha256.txt includes original artifacts and input Parquet hashes.
- build.log includes actual compile/link command lines against latest full-ninja-build static archive.
- run-commands.log includes trial commands, timestamps, observed competing processes and linked archive SHA256.
- old-heap.log and compact-heap.log include all component data without truncation.
- *probe-analysis*.log contain retained RSS, wall time and /usr/bin/time -v statistics.
- final-sha256.txt records final evidence artifact identities.

New linked libOpenDIA_impl.a SHA256: d651ace1b1ab16b8fa722e3b4bb000a47bd9563d58db6545c3fefc7203387bbe
Old heap executable SHA256: 7a6e1b0c049782f91be1f796e694026ba86e01bb1c2d4278c80a3b3d46a3d143
Old probe executable SHA256: 79cf5d9cd8267138ed0f58337a91930e9563311c37a902af182a74049f8ea114
Compact heap executable SHA256: 3be1d668a148b95fed346cddf7f8b28a1633cdd1be64fa35b9a0031c38133f29
Compact probe executable SHA256: af23acb3df0df18fd17fdf893108d3e92abc7926c0f4ccab0e8a4c75c10b9660
