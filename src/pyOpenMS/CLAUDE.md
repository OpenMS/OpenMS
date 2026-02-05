# pyOpenMS Development

## Architecture

- **Bindings** (`bindings/`) — hand-maintained nanobind C++ binding files (10 domain modules + main)
- **Type casters** (`bindings/type_casters/`) — custom nanobind type casters (C++ ↔ Python conversion)
- **Addon system** (`pyopenms/addons/`) — pure Python methods injected into wrapped classes at import time
- **DataFrame wrappers** in `pyopenms/_dataframes.py` add pandas methods (keeps pandas optional)

## Build

```bash
# Build (as part of OpenMS)
cmake --build OpenMS-build --target pyopenms -j$(nproc)

# Test
PYTHONPATH=OpenMS-build/pyOpenMS python3 -m pytest src/pyOpenMS/tests/ -v

# Standalone wheel build
cd src/pyOpenMS && pip wheel . --no-build-isolation
```

## Wrapping New C++ Classes

1. Pick the right `bindings/bind_<domain>.cpp` based on the C++ header path:
   - `KERNEL/` → `bind_kernel.cpp`
   - `FORMAT/` → `bind_format.cpp`
   - `ANALYSIS/` → `bind_analysis.cpp`
   - `CHEMISTRY/` → `bind_chemistry.cpp`
   - `METADATA/` → `bind_metadata.cpp`
   - `PROCESSING/` → `bind_processing.cpp`
   - `FEATUREFINDER/` → `bind_featurefinder.cpp`
   - `DATASTRUCTURES/`, `MATH/`, `CONCEPT/` → `bind_datastructures.cpp`
   - `ML/` → `bind_ml.cpp`
   - Everything else → `bind_misc.cpp`
2. Add the `#include` for the C++ header
3. Add `nb::class_<...>(m, "ClassName", "docstring")` with `.def()` chains
4. Build and test

Each class has a section comment (`// --- ClassName ---`) for navigation.

## Addon Rules

- Pure Python methods in `pyopenms/addons/`
- Use the `@addon("ClassName")` decorator from `pyopenms.addons`
- Keep addons minimal — only add non-auto-generatable methods
- Performance-critical methods should be C++ lambdas in the bindings

## Module Split

10 domain-based nanobind modules for parallel compilation:
`_pyopenms_kernel`, `_pyopenms_metadata`, `_pyopenms_chemistry`, `_pyopenms_analysis`,
`_pyopenms_featurefinder`, `_pyopenms_format`, `_pyopenms_processing`,
`_pyopenms_datastructures`, `_pyopenms_ml`, `_pyopenms_misc`

Plus `_pyopenms` (main entry) and `_arrow_zerocopy` (when WITH_PARQUET=ON).

All modules use `NB_DOMAIN "pyopenms"` to share type information.

## Type Handling

- nanobind type casters auto-convert `OpenMS::String` ↔ Python `str`
- `PeptideIdentificationList` required (not Python list) for `setPeptideIdentifications()`
- `AASequence.fromString()`: valid amino acids only (A-Z except B, J, O, U, X, Z)
- `DPosition<1>` accepts `float`, `DPosition<2>` accepts `tuple`

## Gotchas

- `getDriftTimeUnitAsString()` returns `'<NONE>'` not empty string
- Mobilogram has no MetaInfoInterface (no `getKeys()`, `getMetaValue()`)
- Test from `/tmp` to avoid path shadowing with source `pyopenms/`

## Testing

```bash
# All tests (from /tmp to avoid import shadowing)
cd /tmp && PYTHONPATH=.../OpenMS-build/pyOpenMS python3 -m pytest .../src/pyOpenMS/tests/ -v

# Legacy tests only
python3 -m pytest src/pyOpenMS/tests/unittests/ -v

# Nanobind-specific tests
python3 -m pytest src/pyOpenMS/tests/nanobind/ -v
```
