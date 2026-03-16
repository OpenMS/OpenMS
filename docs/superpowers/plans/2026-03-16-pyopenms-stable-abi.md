# pyOpenMS Stable ABI Migration Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build pyOpenMS against Python's Stable ABI (abi3) so one wheel per platform works across Python 3.12+, eliminating per-version wheel builds.

**Architecture:** Enable nanobind's `STABLE_ABI` flag on all 15 extension modules, replace 8 `PyList_SET_ITEM` calls with stable-ABI-compatible `PyList_SetItem`, and update Python version requirements from 3.9/3.10 to 3.12. A benchmark suite captures before/after performance to gate the merge.

**Tech Stack:** nanobind 2.10.0, CMake 3.26+, Python 3.12+ Stable ABI, cibuildwheel, pytest

**Spec:** `docs/superpowers/specs/2026-03-16-pyopenms-stable-abi-design.md`

---

## Chunk 1: Benchmarks

### Task 1: Add type caster round-trip benchmarks

**Files:**
- Modify: `src/pyOpenMS/tests/benchmark_pyopenms.py`

These benchmarks exercise the custom type caster hot paths most likely affected by stable ABI overhead.

- [ ] **Step 1: Write type caster benchmark function**

Add a new function `bench_type_casters(suite)` following the existing pattern (see `bench_type_conversions` at line 942 for reference). Insert before `def main():` (around line 1079). Uses `suite.bench(name, category, fn, iterations=N)` — note argument order is `(name, category)`.

```python
def bench_type_casters(suite: BenchmarkSuite):
    """Benchmark custom type caster round-trip performance (stable ABI sensitive)."""
    import pyopenms
    print("\n[TYPE CASTERS]")

    iters = 500 if suite.quick else 2000

    # --- OpenMS::String <-> Python str ---
    ef = pyopenms.EmpiricalFormula("C6H12O6")
    suite.bench("String round-trip (1000x toString)", "Type Casters",
                lambda: [ef.toString() for _ in range(1000)], iterations=iters)

    # --- DataValue variant round-trips ---
    spec = pyopenms.MSSpectrum()
    spec.setMetaValue("test_int", 42)
    suite.bench("DataValue int round-trip (1000x)", "Type Casters",
                lambda: [spec.getMetaValue("test_int") for _ in range(1000)], iterations=iters)

    spec.setMetaValue("test_float", 3.14159)
    suite.bench("DataValue float round-trip (1000x)", "Type Casters",
                lambda: [spec.getMetaValue("test_float") for _ in range(1000)], iterations=iters)

    spec.setMetaValue("test_str", "hello world")
    suite.bench("DataValue string round-trip (1000x)", "Type Casters",
                lambda: [spec.getMetaValue("test_str") for _ in range(1000)], iterations=iters)

    # --- DataValue list types (exercise PyList_SetItem paths) ---
    spec.setMetaValue("test_sl", ["alpha", "beta", "gamma", "delta", "epsilon"])
    suite.bench("DataValue StringList round-trip (500x)", "Type Casters",
                lambda: [spec.getMetaValue("test_sl") for _ in range(500)], iterations=iters)

    spec.setMetaValue("test_il", [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    suite.bench("DataValue IntList round-trip (500x)", "Type Casters",
                lambda: [spec.getMetaValue("test_il") for _ in range(500)], iterations=iters)

    spec.setMetaValue("test_dl", [1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8, 9.9, 10.0])
    suite.bench("DataValue DoubleList round-trip (500x)", "Type Casters",
                lambda: [spec.getMetaValue("test_dl") for _ in range(500)], iterations=iters)

    # --- vector<String> conversion (openms_stl_caster.h) ---
    spec.setMetaValue("keys_test", "value")
    suite.bench("getKeys() vector<String> conversion", "Type Casters",
                lambda: spec.getKeys(), iterations=iters)

    # --- map<String, DataValue> conversion (openms_stl_caster.h) ---
    for i in range(20):
        spec.setMetaValue(f"key_{i}", float(i))

    def get_all_meta():
        keys = spec.getKeys()
        return {k: spec.getMetaValue(k) for k in keys}

    suite.bench("map<String,DataValue> 20-entry round-trip", "Type Casters",
                get_all_meta, iterations=iters)
```

- [ ] **Step 2: Wire up the new category in main()**

In `main()` (around line 1199), add before the summary section:

```python
    if should_run("type caster"):
        bench_type_casters(suite)
```

- [ ] **Step 3: Run benchmarks to verify they work**

Run:
```bash
PYTHONPATH=OpenMS-build/pyOpenMS python3 src/pyOpenMS/tests/benchmark_pyopenms.py --quick --categories type_casters
```
Expected: All type_casters benchmarks complete without errors, timing output printed.

- [ ] **Step 4: Commit**

```bash
git add src/pyOpenMS/tests/benchmark_pyopenms.py
git commit -m "bench: add type caster round-trip benchmarks for stable ABI comparison"
```

---

### Task 2: Add Parquet I/O benchmarks

**Files:**
- Modify: `src/pyOpenMS/tests/benchmark_pyopenms.py`

**Depends on:** Task 1 (same file, extends benchmark suite)

- [ ] **Step 1: Write Parquet I/O benchmark function**

Add after `bench_type_casters`. These exercise the Arrow C Data Interface bridge. Guard with early return if `XIMParquetFile` is unavailable.

```python
def bench_parquet_io(suite: BenchmarkSuite, exp):
    """Benchmark Parquet/Arrow I/O paths (requires WITH_PARQUET build)."""
    import pyopenms
    print("\n[PARQUET I/O]")

    if not hasattr(pyopenms, 'XIMParquetFile'):
        print("  (skipped - WITH_PARQUET not enabled)")
        return

    import os
    iters = 100 if suite.quick else 500

    # Try to find the XIM test data file
    script_dir = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.join(script_dir, "..", "..", "tests", "class_tests", "openms", "data",
                     "XIMParquetFile_23_input.xim"),
    ]
    xim_path = None
    for c in candidates:
        if os.path.exists(c):
            xim_path = os.path.abspath(c)
            break

    if xim_path:
        xim = pyopenms.XIMParquetFile()
        xim.load(xim_path)

        suite.bench("XIMParquetFile.to_df()", "Parquet I/O",
                    lambda: xim.to_df(), iterations=iters)

        suite.bench("XIMParquetFile.to_arrow()", "Parquet I/O",
                    lambda: xim.to_arrow(), iterations=iters)

        suite.bench("XIMParquetFile.query_mobilograms().to_df()", "Parquet I/O",
                    lambda: xim.query_mobilograms().to_df(), iterations=iters)
    else:
        print("  (skipped XIM benchmarks - test data not found)")

    # Arrow zero-copy export from MSExperiment (uses existing exp fixture)
    if hasattr(exp, 'to_arrow'):
        suite.bench("MSExperiment.to_arrow(long)", "Parquet I/O",
                    lambda: exp.to_arrow(format="long"), iterations=iters)

        suite.bench("MSExperiment.to_arrow(semi_wide)", "Parquet I/O",
                    lambda: exp.to_arrow(format="semi_wide"), iterations=iters)
```

- [ ] **Step 2: Wire up in main()**

In `main()`, add after the type casters block:

```python
    if should_run("parquet"):
        bench_parquet_io(suite, exp)
```

- [ ] **Step 3: Run to verify**

Run:
```bash
PYTHONPATH=OpenMS-build/pyOpenMS python3 src/pyOpenMS/tests/benchmark_pyopenms.py --quick --categories parquet_io
```
Expected: parquet_io benchmarks complete (or print skip message if no parquet support).

- [ ] **Step 4: Commit**

```bash
git add src/pyOpenMS/tests/benchmark_pyopenms.py
git commit -m "bench: add Parquet/Arrow I/O benchmarks for stable ABI comparison"
```

---

### Task 3: Add benchmark comparison utility

**Files:**
- Create: `src/pyOpenMS/tests/compare_benchmarks.py`

- [ ] **Step 1: Write the comparison script**

```python
#!/usr/bin/env python3
"""Compare two benchmark JSON result files and report regressions.

Usage:
    python compare_benchmarks.py baseline.json stable_abi.json
    python compare_benchmarks.py baseline.json stable_abi.json --threshold 5.0
"""
import argparse
import json
import sys
from collections import defaultdict


def load_results(path):
    """Load benchmark JSON and return {category: {name: median_s}}."""
    with open(path) as f:
        data = json.load(f)
    by_category = defaultdict(dict)
    for result in data.get("results", []):
        by_category[result["category"]][result["name"]] = result["median_s"]
    return dict(by_category)


def compare(baseline, current, warn_threshold=5.0, fail_threshold=15.0):
    """Compare two result dicts. Returns (rows, overall_delta, passed)."""
    rows = []
    all_deltas = []

    all_categories = sorted(set(baseline) | set(current))
    for cat in all_categories:
        base_benchmarks = baseline.get(cat, {})
        curr_benchmarks = current.get(cat, {})
        all_names = sorted(set(base_benchmarks) | set(curr_benchmarks))

        for name in all_names:
            base_val = base_benchmarks.get(name)
            curr_val = curr_benchmarks.get(name)

            if base_val is None:
                rows.append((cat, name, None, curr_val, None, "NEW"))
                continue
            if curr_val is None:
                rows.append((cat, name, base_val, None, None, "REMOVED"))
                continue

            if base_val > 0:
                delta_pct = ((curr_val - base_val) / base_val) * 100
            else:
                delta_pct = 0.0

            all_deltas.append(delta_pct)

            if delta_pct > fail_threshold:
                status = "FAIL"
            elif delta_pct > warn_threshold:
                status = "WARN"
            else:
                status = "PASS"

            rows.append((cat, name, base_val, curr_val, delta_pct, status))

    overall_delta = sum(all_deltas) / len(all_deltas) if all_deltas else 0.0
    passed = all(r[5] in ("PASS", "NEW") for r in rows)
    return rows, overall_delta, passed


def print_report(rows, overall_delta, passed):
    """Print color-coded comparison table."""
    colors = {
        "PASS": "\033[32m",  # green
        "WARN": "\033[33m",  # yellow
        "FAIL": "\033[31m",  # red
        "NEW": "\033[36m",   # cyan
        "REMOVED": "\033[90m",  # gray
    }
    reset = "\033[0m"

    print(f"\n{'Category':<25} {'Benchmark':<50} {'Base (ms)':>10} {'Curr (ms)':>10} {'Delta':>8} {'Status':>8}")
    print("-" * 115)

    for cat, name, base_val, curr_val, delta_pct, status in rows:
        color = colors.get(status, "")
        base_str = f"{base_val * 1000:.3f}" if base_val is not None else "N/A"
        curr_str = f"{curr_val * 1000:.3f}" if curr_val is not None else "N/A"
        delta_str = f"{delta_pct:+.1f}%" if delta_pct is not None else "---"
        print(f"{cat:<25} {name:<50} {base_str:>10} {curr_str:>10} {color}{delta_str:>8} {status:>8}{reset}")

    print("-" * 115)
    overall_color = colors["PASS"] if passed else colors["FAIL"]
    print(f"Overall mean delta: {overall_color}{overall_delta:+.1f}%{reset}")
    print(f"Result: {overall_color}{'PASSED' if passed else 'FAILED'}{reset}\n")


def main():
    parser = argparse.ArgumentParser(description="Compare benchmark results")
    parser.add_argument("baseline", help="Path to baseline JSON results")
    parser.add_argument("current", help="Path to current JSON results")
    parser.add_argument("--warn-threshold", type=float, default=5.0,
                        help="Percent regression to trigger warning (default: 5.0)")
    parser.add_argument("--fail-threshold", type=float, default=15.0,
                        help="Percent regression to trigger failure (default: 15.0)")
    args = parser.parse_args()

    baseline = load_results(args.baseline)
    current = load_results(args.current)
    rows, overall_delta, passed = compare(baseline, current, args.warn_threshold, args.fail_threshold)
    print_report(rows, overall_delta, passed)

    sys.exit(0 if passed else 1)


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Test with a synthetic pair of JSON files**

Run:
```bash
python3 -c "
import json
base = {'results': [
    {'category': 'test', 'name': 'bench1', 'median_s': 0.010, 'mean_s': 0.011, 'stdev_s': 0.001, 'min_s': 0.009},
    {'category': 'test', 'name': 'bench2', 'median_s': 0.020, 'mean_s': 0.021, 'stdev_s': 0.002, 'min_s': 0.019},
]}
curr = {'results': [
    {'category': 'test', 'name': 'bench1', 'median_s': 0.0105, 'mean_s': 0.011, 'stdev_s': 0.001, 'min_s': 0.009},
    {'category': 'test', 'name': 'bench2', 'median_s': 0.025, 'mean_s': 0.026, 'stdev_s': 0.002, 'min_s': 0.024},
]}
json.dump(base, open('/tmp/base.json', 'w'))
json.dump(curr, open('/tmp/curr.json', 'w'))
"
python3 src/pyOpenMS/tests/compare_benchmarks.py /tmp/base.json /tmp/curr.json
```
Expected: Table showing bench1 as PASS (~5%), bench2 as WARN (~25%), overall FAILED.

- [ ] **Step 3: Commit**

```bash
git add src/pyOpenMS/tests/compare_benchmarks.py
git commit -m "bench: add benchmark comparison utility for stable ABI regression detection"
```

---

### Task 4: Capture baseline benchmarks

**Files:** None (manual step, produces JSON artifact)

- [ ] **Step 1: Run full benchmark suite and save baseline**

Run:
```bash
PYTHONPATH=OpenMS-build/pyOpenMS python3 src/pyOpenMS/tests/benchmark_pyopenms.py --json /tmp/stable_abi_baseline.json
```
Expected: All categories complete, JSON file written.

- [ ] **Step 2: Verify JSON contains all categories including new ones**

Run:
```bash
python3 -c "import json; data=json.load(open('/tmp/stable_abi_baseline.json')); cats=set(r['category'] for r in data['results']); print(sorted(cats))"
```
Expected: Output includes `type_casters` and `parquet_io` among the categories.

---

## Chunk 2: Code Changes

### Task 5: Replace PyList_SET_ITEM with PyList_SetItem in type casters

**Files:**
- Modify: `src/pyOpenMS/bindings/type_casters/openms_stl_caster.h:119`
- Modify: `src/pyOpenMS/bindings/type_casters/openms_datavalue_caster.h:255,270,285,500,515,530`
- Modify: `src/pyOpenMS/bindings/type_casters/openms_dposition_caster.h:243`

Both `PyList_SET_ITEM` and `PyList_SetItem` steal a reference. `PyList_SetItem` additionally `Py_XDECREF`s the old item at that index, but since all call sites operate on freshly created lists from `PyList_New`, every slot is `NULL` so the extra `Py_XDECREF` is a no-op. This is a safe drop-in replacement.

- [ ] **Step 1: Replace in openms_stl_caster.h**

At line 119, change:
```cpp
PyList_SET_ITEM(list, i, item);
```
to:
```cpp
PyList_SetItem(list, i, item);
```

- [ ] **Step 2: Replace in openms_datavalue_caster.h**

Replace all 6 occurrences. At lines 255, 270, 285, 500, 515, 530, change:
```cpp
PyList_SET_ITEM(list, i, item);
```
to:
```cpp
PyList_SetItem(list, i, item);
```

- [ ] **Step 3: Replace in openms_dposition_caster.h**

At line 243, change:
```cpp
PyList_SET_ITEM(list, i, tuple);  // steals reference
```
to:
```cpp
PyList_SetItem(list, i, tuple);  // steals reference
```

- [ ] **Step 4: Build and run tests to verify no regressions**

Run:
```bash
cmake --build OpenMS-build --target pyopenms -j$(nproc)
PYTHONPATH=OpenMS-build/pyOpenMS python3 -m pytest src/pyOpenMS/tests/ -x -q --import-mode=importlib 2>&1 | tail -5
```
Expected: All tests pass (2904 passed, 0 failed).

- [ ] **Step 5: Commit**

```bash
git add src/pyOpenMS/bindings/type_casters/openms_stl_caster.h \
        src/pyOpenMS/bindings/type_casters/openms_datavalue_caster.h \
        src/pyOpenMS/bindings/type_casters/openms_dposition_caster.h
git commit -m "fix: replace PyList_SET_ITEM with PyList_SetItem for stable ABI compat

PyList_SET_ITEM is not in Python's Limited API. PyList_SetItem is the
stable-ABI-compatible equivalent (both steal a reference). All call
sites operate on freshly created lists so the behavioral difference
(bounds checking + old-item decref) is a no-op.

Ref: #8185"
```

---

### Task 6: Enable STABLE_ABI in CMakeLists.txt

**Files:**
- Modify: `src/pyOpenMS/CMakeLists.txt:14,127,267,276,314`

- [ ] **Step 1: Bump cmake_minimum_required to 3.26**

At line 14, change:
```cmake
cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
```
to:
```cmake
cmake_minimum_required(VERSION 3.26 FATAL_ERROR)
```

- [ ] **Step 2: Update find_package to use Development.SABIModule and Python 3.12**

At line 127, change:
```cmake
find_package(Python 3.9 REQUIRED COMPONENTS Interpreter Development.Module NumPy)
```
to:
```cmake
find_package(Python 3.12 REQUIRED COMPONENTS Interpreter Development.SABIModule NumPy)
```

- [ ] **Step 3: Add STABLE_ABI to domain module loop**

At lines 267-271, change:
```cmake
  nanobind_add_module(${MODULE_NAME}
    NB_DOMAIN pyopenms  # Share types across modules
    NOMINSIZE            # Use CMake's default -O3 instead of nanobind's -Os
    ${PYOPENMS_BINDINGS_DIR}/bind_${domain}.cpp
  )
```
to:
```cmake
  nanobind_add_module(${MODULE_NAME}
    NB_DOMAIN pyopenms   # Share types across modules
    STABLE_ABI           # Python Limited/Stable ABI (abi3)
    NOMINSIZE            # Use CMake's default -O3 instead of nanobind's -Os
    ${PYOPENMS_BINDINGS_DIR}/bind_${domain}.cpp
  )
```

- [ ] **Step 4: Add STABLE_ABI to main module**

At lines 276-280, change:
```cmake
nanobind_add_module(_pyopenms
  NB_DOMAIN pyopenms
  NOMINSIZE
  ${PYOPENMS_BINDINGS_DIR}/main_module.cpp
)
```
to:
```cmake
nanobind_add_module(_pyopenms
  NB_DOMAIN pyopenms
  STABLE_ABI
  NOMINSIZE
  ${PYOPENMS_BINDINGS_DIR}/main_module.cpp
)
```

- [ ] **Step 5: Add STABLE_ABI to arrow module**

At lines 314-318, change:
```cmake
nanobind_add_module(_arrow_zerocopy
  NB_DOMAIN pyopenms
  NOMINSIZE
  ${CMAKE_CURRENT_SOURCE_DIR}/bindings/arrow_zerocopy.cpp
)
```
to:
```cmake
nanobind_add_module(_arrow_zerocopy
  NB_DOMAIN pyopenms
  STABLE_ABI
  NOMINSIZE
  ${CMAKE_CURRENT_SOURCE_DIR}/bindings/arrow_zerocopy.cpp
)
```

- [ ] **Step 6: Rebuild and verify abi3 extension suffixes**

Run:
```bash
cmake --build OpenMS-build --target pyopenms -j$(nproc)
ls OpenMS-build/pyOpenMS/pyopenms/_pyopenms*.so 2>/dev/null || ls OpenMS-build/pyOpenMS/pyopenms/_pyopenms*.pyd 2>/dev/null
```
Expected: Extensions have `.abi3.so` suffix (e.g., `_pyopenms.abi3.so`) instead of `.cpython-3XX-linux-gnu.so`.

- [ ] **Step 7: Run full test suite**

Run:
```bash
PYTHONPATH=OpenMS-build/pyOpenMS python3 -m pytest src/pyOpenMS/tests/ -x -q --import-mode=importlib 2>&1 | tail -5
```
Expected: All tests pass.

- [ ] **Step 8: Commit**

```bash
git add src/pyOpenMS/CMakeLists.txt
git commit -m "feat: enable Python Stable ABI (abi3) for pyOpenMS nanobind modules

Add STABLE_ABI flag to all 15 nanobind_add_module() calls. Bump
cmake_minimum_required to 3.26 (needed for Development.SABIModule)
and find_package Python minimum to 3.12 (nanobind STABLE_ABI requirement).

This produces .abi3.so extensions that work across Python 3.12+,
eliminating the need to build separate wheels per Python version.

Ref: #8185"
```

---

### Task 7: Update Python version requirements in pyproject.toml

**Files:**
- Modify: `src/pyOpenMS/pyproject.toml:33,47-60,151-155`

- [ ] **Step 1: Update requires-python**

At line 33, change:
```toml
requires-python = ">=3.10"
```
to:
```toml
requires-python = ">=3.12"
```

- [ ] **Step 2: Remove 3.10 and 3.11 classifiers**

At lines 47-60, remove the two lines:
```toml
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
```

- [ ] **Step 3: Update cibuildwheel build targets**

At line 153, change:
```toml
build = "cp310-* cp311-* cp312-* cp313-* cp314-*"
```
to:
```toml
build = "cp312-*"
```

Note: cibuildwheel builds only for cp312, producing an abi3 wheel. The test matrix (handled separately in CI) still verifies the wheel on 3.12, 3.13, 3.14.

- [ ] **Step 4: Commit**

```bash
git add src/pyOpenMS/pyproject.toml
git commit -m "feat: update Python version requirements to 3.12+ for stable ABI

Drop Python 3.10/3.11 support. Build only cp312 wheels (abi3 wheels
work on all 3.12+ versions). Update classifiers accordingly.

Ref: #8185"
```

---

## Chunk 3: CI & Verification

### Task 8: Update CI workflow for abi3 wheel build/test split

**Files:**
- Modify: `.github/workflows/python_versions.json`
- Modify: `.github/workflows/pyopenms-wheels-cibuildwheel.yml:56-66,475-479`

- [ ] **Step 1: Update python_versions.json**

Change `.github/workflows/python_versions.json` from:
```json
[
  "3.10", "3.11", "3.12", "3.13", "3.14"
]
```
to:
```json
{
  "build": ["3.12"],
  "test": ["3.12", "3.13", "3.14"]
}
```

- [ ] **Step 2: Update setup job to use split build/test versions**

In the workflow file setup job (around lines 51-66), make these specific changes:

**Line 51** — currently outputs the full JSON for the test matrix:
```bash
# Before:
echo "python-versions=$(jq -c . .github/workflows/python_versions.json)" >> $GITHUB_OUTPUT
# After:
echo "python-versions=$(jq -c '.test' .github/workflows/python_versions.json)" >> $GITHUB_OUTPUT
```

**Line 53** — display line, update to show build versions:
```bash
# Before:
echo "Python versions: $(jq -r '.[]' .github/workflows/python_versions.json | tr '\n' ' ')"
# After:
echo "Build versions: $(jq -r '.build[]' .github/workflows/python_versions.json | tr '\n' ' ')"
echo "Test versions: $(jq -r '.test[]' .github/workflows/python_versions.json | tr '\n' ' ')"
```

**Lines 56, 59, 62, 65** — four platform-specific CIBW_BUILD lines, change each `jq -r '.[]'` to `jq -r '.build[]'`:
```bash
# Before (repeated for each platform):
CIBW_LINUX_X64=$(jq -r '.[]' .github/workflows/python_versions.json | sed ...)
# After:
CIBW_LINUX_X64=$(jq -r '.build[]' .github/workflows/python_versions.json | sed ...)
```

- [ ] **Step 3: Update test-wheels job matrix to use python-versions output**

The test-wheels job matrix (around line 448) already uses `${{ fromJson(needs.setup.outputs.python-versions) }}`. Since we changed the `python-versions` output in Step 2 to use `.test`, this will now correctly produce `["3.12", "3.13", "3.14"]`. No change needed here — verify it references the setup job output.

- [ ] **Step 4: Update test-wheels pip install command**

At line 475, change:
```bash
pip install wheelhouse/*cp${PYVER}*.whl
```
to:
```bash
pip install wheelhouse/*.whl
```

Since there's only one abi3 wheel per platform, no version-specific glob is needed.

- [ ] **Step 5: Add explicit extension module import checks**

After the existing import checks (lines 476-478), add:
```bash
python -c "from pyopenms._pyopenms import __doc__; print('_pyopenms module OK')"
python -c "
try:
    from pyopenms._arrow_zerocopy import featuremap_features_to_arrow
    print('_arrow_zerocopy module OK')
except ImportError:
    print('_arrow_zerocopy not available (WITH_PARQUET off)')
"
```

- [ ] **Step 6: Commit**

```bash
git add .github/workflows/python_versions.json \
        .github/workflows/pyopenms-wheels-cibuildwheel.yml
git commit -m "ci: split build/test matrix for abi3 wheels

Build only cp312 (abi3 wheel works across versions). Test on 3.12,
3.13, 3.14 to verify cross-version compatibility. Add explicit
import checks for _pyopenms and _arrow_zerocopy extension modules.

Ref: #8185"
```

---

### Task 9: Capture stable ABI benchmarks and compare

**Files:** None (manual step, uses artifacts from Tasks 1-4)

- [ ] **Step 1: Run full benchmark suite on stable ABI build**

Run:
```bash
PYTHONPATH=OpenMS-build/pyOpenMS python3 src/pyOpenMS/tests/benchmark_pyopenms.py --json /tmp/stable_abi_current.json
```

- [ ] **Step 2: Compare against baseline**

Run:
```bash
python3 src/pyOpenMS/tests/compare_benchmarks.py /tmp/stable_abi_baseline.json /tmp/stable_abi_current.json
```
Expected: Table showing per-benchmark deltas. Go/no-go:
- All PASS (< 5% regression): proceed to merge
- Any WARN (5-15%): investigate specific benchmarks, document findings
- Any FAIL (> 15%) or overall > 10%: do not merge, investigate root cause

- [ ] **Step 3: Document results**

If results pass, add a brief summary to the PR description with the key numbers (overall delta, any notable category-level results).
