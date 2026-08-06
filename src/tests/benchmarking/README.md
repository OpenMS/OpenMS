# Benchmark Tests (Issue #8788)

This directory will contain **benchmark tests** for OpenMS: tests that run
real TOPP tool pipelines on **large, realistic datasets** (multi-GB, with
known ground truth) and compare scientific + resource metrics against
committed baselines.

## Status: PR 1 skeleton

Currently only the skeleton is present:

- CMake option `ENABLE_BENCHMARK_TESTING` (OFF by default) in
  `src/tests/CMakeLists.txt` (matching the `ENABLE_TOPP_TESTING` convention)
- This directory with a single placeholder test (`BENCH_SKELETON`) proving
  the CTest wiring

## How to enable

```bash
cmake -DENABLE_BENCHMARK_TESTING=ON ...
ctest -N -L benchmark   # lists the registered benchmark tests
ctest -L benchmark      # runs them
```

## Naming scheme

Benchmark tests use the `BENCH_` prefix (`BENCH_<tool>_<dataset>[_<num>]`),
mirroring the `TOPP_` prefix convention: it distinguishes benchmark tests
from unit/TOPP tests in nightly builds and CDash so benchmark results can be
isolated and trended over time (the "run regularly" goal of Issue #8788).

## Design notes

- **Datasets are not committed.** They are referenced in a dataset manifest
  (accession / URL / SHA256) and downloaded on demand, following the existing
  `FetchContent` + `URL_HASH` pattern already used for
  `ENABLE_OPENTIMS_TESTS` / `ENABLE_THERMO_RAW_TESTS`
  (see top-level `CMakeLists.txt`).
- **Benchmarks extend the validation stage**, not the execution stage:
  run the existing TOPP tools, extract metrics from their outputs
  (mzTab / featureXML / idXML), compare against baselines.
- **Scheduled runs** (nightly CI) come in a later PR; these tests are
  intentionally excluded from per-PR test runs.

- **TOPP tools dependency:** `ENABLE_BENCHMARK_TESTING` is independent of
  `BUILD_TOPP_TOOLS` for now (the skeleton runs no TOPP tool). Once real
  benchmark tests execute TOPP tools, they must also be covered by the
  `BUILD_TOPP_TOOLS` FATAL_ERROR guard in `src/tests/CMakeLists.txt`.

## Adding a benchmark (later PRs)

1. Add the dataset to the manifest with checksum.
2. Register a `BENCH_<tool>_<dataset>` test in `CMakeLists.txt`.
3. Commit a baseline YAML with expected metric ranges.
