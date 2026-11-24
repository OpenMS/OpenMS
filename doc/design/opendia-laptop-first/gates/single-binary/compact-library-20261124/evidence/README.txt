Measurement session 2026-11-24, all paths under ~/.jcode/scratch/compact-20261124.
No repository edits or commits.
Before building: cp -a ~/.jcode/scratch/perf/{heap,probe,heap.cpp,probe.cpp,heap.o,probe.o,build_heap.sh,build_probe.sh} old/
Old initial commands, sequential:
/usr/bin/time -v stdbuf -oL old/heap ~/.jcode/scratch/perf/lib/search.opendia.parquet ~/.jcode/scratch/perf/lib/scoring.opendia.parquet > old-heap.log 2>&1
/usr/bin/time -v stdbuf -oL old/probe analysis ~/.jcode/scratch/perf/lib/search.opendia.parquet ~/.jcode/scratch/perf/lib/scoring.opendia.parquet > old-probe-analysis.log 2>&1
Build: bash build_compact.sh > build.log 2>&1
Trials: bash run_measurements.sh > run-commands.log 2>&1
Memory MB labels in old harness actually mean MiB (2^20 bytes). Old probe rss_mb is integer-truncated MiB from SysInfo. New probe retained_rss_MiB uses /proc/self/statm pages.
Heap = mallinfo2().uordblks + mallinfo2().hblkhd, with process pre-load baseline subtracted. Components are incremental allocation from copies, not exact owned-capacity accounting. Copies of vectors allocate size rather than source capacity. Precursor string_view IDs in compact precursor copy borrow the original map keys, as intended.
Compact heap runs after standalone search/scoring reads and frees, so use clean compact_probe analysis for representative retained RSS/load time. Old probe reports analysis before an extra search Library copy, so /usr/bin/time peak covers the copy too. Internal old peak_mb was measured before that copy. Old heap's total reported before component copies.
All micro timings are concurrent, not isolated benchmarking: root identity raw01 run uses 16 threads and can perturb CPU and memory bandwidth. Exact observed processes recorded in run-commands.log. Shared library dependencies remain current build libraries for both old and new. Old static OpenDIA implementation is embedded in preserved binary. Input hashes in initial-sha256.txt. Linked new archive hash in run-commands.log.
