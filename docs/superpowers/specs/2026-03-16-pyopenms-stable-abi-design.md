# pyOpenMS Stable ABI Migration

**Issue:** [#8185](https://github.com/OpenMS/OpenMS/issues/8185)
**Date:** 2026-03-16
**Status:** Draft

## Goal

Build pyOpenMS against Python's Stable ABI (abi3) so that one wheel per platform works across Python minor versions, eliminating the need to build for each Python version separately.

## Constraints

- Minimum Python version rises from 3.9/3.10 (CMake/pyproject.toml) to 3.12 (nanobind's `STABLE_ABI` requires 3.12+)
- CMake minimum version rises to 3.26 (`Development.SABIModule` component requires it)
- Must benchmark before merging — especially zero-copy Arrow paths
- Performance regression > 15% on zero-copy/Arrow or > 10% across the board is a no-go

## Approach: Benchmark-Gated Migration

Write benchmarks first, establish baseline, make the (small) code changes, re-benchmark, merge only if performance is acceptable.

## Phase 1: Benchmark Baseline

### Existing coverage

The benchmark suite at `src/pyOpenMS/tests/benchmark_pyopenms.py` already covers:
- Memory Views / Zero-Copy (get_data vs get_data_view, capsule ownership)
- Arrow Export (to_arrow with format/filter variants)
- Type Conversions (AASequence, EmpiricalFormula, etc.)

### New benchmarks to add

**Type Caster Round-Trip category:**
- Bulk `OpenMS::String` <-> Python `str` conversion (1000+ strings)
- `DataValue` variant round-trips (int, float, string, string-list variants)
- `DPosition<2>` tuple packing/unpacking
- `std::vector<String>` and `std::map<String, V>` conversions

**Parquet I/O sub-category:**
- `XIMParquetFile` read and query operations (exercises Arrow C Data Interface bridge)

### Baseline capture

Run full benchmark suite on current non-abi3 build:
```bash
PYTHONPATH=OpenMS-build/pyOpenMS python3 src/pyOpenMS/tests/benchmark_pyopenms.py --json baseline.json
```

## Phase 2: Code Changes

### CMakeLists.txt (`src/pyOpenMS/CMakeLists.txt`)

Add `STABLE_ABI` flag to all 15 `nanobind_add_module()` calls (13 domain + 1 main + 1 arrow, the last conditional on `WITH_PARQUET`):

```cmake
# Domain modules (in loop over PYOPENMS_DOMAINS):
nanobind_add_module(${MODULE_NAME}
  NB_DOMAIN pyopenms
  STABLE_ABI
  NOMINSIZE
  ${PYOPENMS_BINDINGS_DIR}/bind_${domain}.cpp
)

# Main module:
nanobind_add_module(_pyopenms
  NB_DOMAIN pyopenms
  STABLE_ABI
  NOMINSIZE
  ${PYOPENMS_BINDINGS_DIR}/main_module.cpp
)

# Arrow module (conditional on WITH_PARQUET):
nanobind_add_module(_arrow_zerocopy
  NB_DOMAIN pyopenms
  STABLE_ABI
  NOMINSIZE
  ${CMAKE_CURRENT_SOURCE_DIR}/bindings/arrow_zerocopy.cpp
)
```

Update `find_package` to use `Development.SABIModule` (required by nanobind for stable ABI) and bump minimum Python to 3.12:
```cmake
find_package(Python 3.12 REQUIRED COMPONENTS Interpreter Development.SABIModule NumPy)
```

Bump `cmake_minimum_required` to 3.26 (needed for `Development.SABIModule` support).

### Type caster fixes (3 files in `bindings/type_casters/`)

Replace `PyList_SET_ITEM()` with `PyList_SetItem()` in 8 call sites:
- `openms_stl_caster.h`
- `openms_datavalue_caster.h`
- `openms_dposition_caster.h`

Both functions steal a reference, so no additional ref-counting changes needed. `PyList_SetItem` is the stable-ABI-compatible equivalent that adds bounds checking.

No other CPython C API changes required. All other functions used in the type casters (`PyList_New`, `PyDict_New`, `PySequence_GetItem`, `PyUnicode_AsUTF8AndSize`, `Py_DECREF`, `PyXxx_Check`, etc.) are part of the stable ABI as of Python 3.12.

### pyproject.toml (`src/pyOpenMS/pyproject.toml`)

- `requires-python = ">=3.12"`
- cibuildwheel build targets: change to `cp312-*` only (abi3 wheel built once, works on all 3.12+ versions)
- Remove `Programming Language :: Python :: 3.10` and `3.11` classifiers
- Wheel repair tools (auditwheel/delocate/delvewheel) handle abi3 wheels natively

### CI workflow (`.github/workflows/pyopenms-wheels-cibuildwheel.yml`)

The workflow currently derives both build selectors and test matrix from `python_versions.json`. These must be split:

- **Build matrix:** `cp312-*` only (one abi3 wheel per platform). Either hardcode `CIBW_BUILD=cp312-*` or introduce a separate `build_versions.json` with `["3.12"]`.
- **Test matrix:** `["3.12", "3.13", "3.14"]` — tests the single abi3 wheel on each Python version.
- Wheel filename changes from `*-cp312-cp312-*.whl` to `*-cp312-abi3-*.whl`
- Update test-wheels install command: `pip install wheelhouse/*cp${PYVER}*.whl` → `pip install wheelhouse/*.whl` (only one wheel per platform)

## Phase 3: Benchmark Comparison & Go/No-Go

1. Build pyOpenMS with `STABLE_ABI` enabled
2. Run full benchmark suite: `--json stable_abi.json`
3. Compare using `compare_benchmarks.py` utility

### Go/no-go criteria

| Threshold | Action |
|-----------|--------|
| Median regression < 5% across all categories | Green light: merge |
| 5-15% regression in specific categories | Yellow: investigate, document, maintainer decision |
| > 15% regression in zero-copy/Arrow OR > 10% across the board | Red: don't merge, file upstream issue if warranted |

### Comparison utility

New script `src/pyOpenMS/tests/compare_benchmarks.py`:
- Loads two JSON benchmark result files
- Prints category-level median deltas
- Color-coded pass/warn/fail output
- Reusable for future regression detection

## Phase 4: Testing

- Full pytest suite (`src/pyOpenMS/tests/`) must pass on the abi3 build
- Verify abi3 wheel installs and imports on Python 3.12, 3.13, and 3.14
- CI tests the single wheel on all supported versions
- **Explicit extension module import checks:** The current wheel test only exercises modules auto-imported by `pyopenms/__init__.py` (names matching `_pyopenms_*`). Add explicit import verification for `_pyopenms` (main module) and `_arrow_zerocopy` (conditional, when `WITH_PARQUET` is on) across all target Python versions to catch abi3 linking issues early.

## Risk Assessment

| Risk | Likelihood | Mitigation |
|------|-----------|------------|
| Performance regression in type casters | Low | Benchmark-gated; most hot paths use same C API functions |
| nanobind `NB_DOMAIN` + `STABLE_ABI` interaction | Low | nanobind handles cross-module sharing internally; tested upstream |
| Vectorcall perf difference under stable ABI | Low | Vectorcall IS available in stable ABI on 3.12+; nanobind #1286 tracks edge cases with `tp_call` fallback on < 3.14 |
| Dropping 3.10/3.11 users | Medium | Python 3.12 is 3+ years old; Bioconda/conda-forge can pin |
| `py-build-cmake` standalone builds | Low | Verify that `pip wheel . --no-build-isolation` produces correct abi3 wheel tags |
| NumPy ndarray stable ABI interaction | Low | nanobind's ndarray support abstracts NumPy C API; verify no direct struct access |

## Files Changed

| File | Change |
|------|--------|
| `src/pyOpenMS/CMakeLists.txt` | Add `STABLE_ABI` to module definitions, bump Python minimum |
| `src/pyOpenMS/bindings/type_casters/openms_stl_caster.h` | `PyList_SET_ITEM` -> `PyList_SetItem` |
| `src/pyOpenMS/bindings/type_casters/openms_datavalue_caster.h` | `PyList_SET_ITEM` -> `PyList_SetItem` |
| `src/pyOpenMS/bindings/type_casters/openms_dposition_caster.h` | `PyList_SET_ITEM` -> `PyList_SetItem` |
| `src/pyOpenMS/pyproject.toml` | `requires-python >= 3.12`, update cibuildwheel targets |
| `.github/workflows/pyopenms-wheels-cibuildwheel.yml` | Update Python version matrix |
| `src/pyOpenMS/tests/benchmark_pyopenms.py` | Add type caster + parquet benchmark categories |
| `src/pyOpenMS/tests/compare_benchmarks.py` | New: benchmark comparison utility |
