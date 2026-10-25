# OpenDIA layout experiment

Planning-branch experiment only. Do not merge these files with the WP-01 implementation.
This measures alternatives, not a production-layout change or an assertion that SoA is optimal.

## Question and source trace

Compare the current owning v0 records with practical AoS/SoA alternatives, separating
**metadata factoring**, **numeric layout**, and **accessor/copy overhead**.

The implementation under test is supplied explicitly with `OPENDIA_SOURCE_DIR`:

- `src/topp/OpenDIA/source/OpenDIARecords.cpp`: `RawTile` owns `vector<RawPoint>` and returns complete points by value, including native-ID strings.
- `src/topp/OpenDIA/test/OpenDIARecordsInterop_test.cpp`: current record consumers are boundary tests, not production extraction adapters.
- `src/openms/source/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.cpp`: existing spectrum access materializes separate m/z and intensity arrays.
- `src/openms/source/ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.cpp`: existing extraction visits spectra and sorted target m/z, using separate m/z/intensity/optional IM iterators and narrow open-interval tophat windows.

There is no complete future-engine record access trace yet. These real legacy interfaces
and explicit row-access scenarios constrain this experiment, but cannot establish a
universal winner or future end-to-end engine speedup.

## Comparisons

1. Actual current `RawTile`, compiled from its unchanged source, accessed through its public API.
2. The same full `RawPoint` AoS payload accessed by const reference, separating layout from API copies.
3. Flat SoA retaining per-point metadata, testing column separation alone.
4. Acquisition-factored AoS: metadata once per acquisition, hot numeric point records and offsets.
5. Acquisition-factored SoA: metadata once per acquisition, contiguous numeric columns and offsets.

Every variant also owns the same acquisition directory and offsets, preserving empty
spectra. The current point-only `RawTile` cannot represent empty acquisitions by itself.
Consequently the baseline memory number includes this necessary experimental directory,
not just the unaugmented `RawTile` object.

No float downcasts, native-ID truncation, coordinate rounding, mobility/FAIMS substitution,
or dropped rows are allowed. Metadata is not inherently duplicated by SoA. We compare
factored AoS as well, so savings from deduplication are not misattributed to SoA.

## Checks before timing

- Reconstruct every field of every row and compare to independently retained source data.
- Preserve optional mobility presence/value/unit, discrete FAIMS (including explicit zero),
  native IDs, frame/scan identities, RT, isolation bounds, MS level and intensity.
- Always-active self-checks cover empty/singleton scans, duplicate m/z and exact boundaries,
  missing mobility and long native IDs. NDEBUG must not disable correctness checks.
- Check tophat outputs against the actual public OpenMS extractor, not just another
  experimental loop. Experimental kernels are supporting evidence, not pipeline acceptance.
- Time the actual public tophat kernel separately, including any packing needed by each
  representation to satisfy its existing separate-vector interface.

## Measurement protocol

Use one optimized binary and identical flags for the experimental layouts, actual record
implementation and actual extractor translation unit. Loading and correctness checks are
outside timed query loops. The remaining shared SDK can be a Debug build, which must be
recorded and must not be advertised as full-engine Release measurement.

Each case runs in a fresh child process with one pinned allowed CPU and one thread.
Randomize the order of variants/cases with a recorded seed, use warm-up followed by repeated
measurements, and retain every sample, checksum, command, return code and fixture hash.
Report medians and spreads, not only the best run. Do not change system power settings.

Use actual non-IM, per-point IM and FAIMS mzML fixtures at native size and explicitly
labeled repeated-fixture scales large enough to stress caches. Whole acquisitions are
repeated without inventing a new scientific data distribution. This scaling is not a
representative full-run or giant diaPASEF-frame benchmark.

Workloads include dense numeric scans, sorted narrow/wide extraction, the actual public
kernel with packing, sparse gathers, and complete-row metadata access. Keep their results
separate. A synthetic workload-weighted overall score would hide the tradeoffs.

Memory measurements distinguish the live allocator delta attributable to retained layout
construction from total process RSS and peak RSS. Keep the same loaded source data alive
for all variants, record its baseline, and do not present total process peak as layout
bytes. Linux/glibc allocator statistics include allocator overhead and are not exact
requested payload bytes. The libc documentation defines `uordblks` as allocated chunks
and `hblkhd` as mmap allocation size. Memory snapshots are outside timed loops.
Construction time is diagnostic only: the actual `RawTile` constructor validates inputs,
whereas experimental packing need not repeat its complete production validation.

## Decision rule

Prefer the smallest representation on the speed/memory Pareto frontier for the access
paths actually needed. Favor a simple scan/tile view API over callers manually joining
metadata fields. Keep full-row access where needed and report its reconstruction cost.
Do not rewrite all record families because raw-point scans favor a particular layout.
Construction cost, temporary packing, sparse metadata reads and unsupported future
workloads remain visible in the decision. The best tested raw layout does not settle
hypothesis/evidence/family/score layouts or byte-budget admission.

## Reproduce on Linux/glibc

Use an out-of-tree build and an existing compatible OpenMS SDK. Replace these placeholders
with local paths. Output directories/files must be new, the runner refuses to overwrite evidence.

```bash
cmake -S <planning-checkout>/doc/design/opendia-laptop-first/experiments/layout \
  -B <scratch>/layout-release -DCMAKE_BUILD_TYPE=Release \
  -DOPENDIA_SOURCE_DIR=<implementation-checkout> -DOpenMS_DIR=<OpenMS-build>
cmake --build <scratch>/layout-release -j3
ctest --test-dir <scratch>/layout-release --output-on-failure
python3 <experiment>/run.py --self-test
python3 <experiment>/run.py --binary <scratch>/layout-release/layout_benchmark \
  --source <implementation-checkout> --output <scratch>/layout-runs --cpu <allowed-CPU>
python3 <experiment>/summarize.py <scratch>/layout-runs --output <scratch>/layout-summary.json
```

Default matrix: five layouts, three actual fixtures, native and at-least-1,048,576-point
scales, plus at-least-4,194,304 IM points, three randomized fresh-process trials per case
and seven measurements after a warm-up per workload: **105 child processes**.
Raw stdout/stderr, commands, source/fixture/binary hashes, ordering and CPU settings stay
in the selected evidence directory. `summarize.py` rejects incomplete/unbalanced cases,
different geometry/workloads or cross-layout/trial checksum disagreements.

Observed results will be recorded after the local runs.
