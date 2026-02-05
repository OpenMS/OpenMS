# pyOpenMS Development

## Architecture

- **Generator** (`generator/`) parses `.pxd` files from `pxds/` and emits nanobind C++ code using libclang
- **Type casters** (`bindings/type_casters/`) handle C++ ↔ Python type conversion
- **Generated bindings** (`bindings/generated/`) contain pre-committed nanobind C++ wrapper code (11 files)
- **Addon system** (`pyopenms/addons/`) injects pure Python methods into wrapped classes at import time
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

Bindings compile from pre-committed generated sources in `bindings/generated/*.cpp` — no code generation at build time.

## Regenerating Bindings (maintainer-only)

Only needed when `.pxd` files change or new classes are added:

```bash
cd src/pyOpenMS
python -m generator --pxd-dir pxds \
  --output-dir bindings/generated \
  --openms-include-dir ../../src/openms/include ../../OpenMS-build/src/openms/include
```

Requires `clang` Python bindings (libclang). Clear cache if needed: `rm -rf OpenMS-build/libclang_cache`

## Addon Rules

- Pure Python methods in `pyopenms/addons/` (not Cython `.pyx`)
- Use the `@addon("ClassName")` decorator from `pyopenms.addons`
- Keep addons minimal - only add non-auto-generatable methods
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

## Wrapping New C++ Classes

1. Add a `.pxd` declaration file in `pxds/`
2. Regenerate bindings (see above)
3. Commit the updated `bindings/generated/*.cpp` files

## wrap- Directive Reference

Directives parsed from .pxd files:

| Directive | Purpose |
|-----------|---------|
| `wrap-doc:` | Documentation string |
| `wrap-as:` | Rename method in Python |
| `wrap-ignore` | Skip wrapping this method |
| `wrap-inherits:` | Specify base classes |
| `wrap-instances:` | Template specializations |
| `wrap-upper-limit:` | Bounds checking for `[]` |
| `wrap-hash:` | Enable `__hash__` |
| `wrap-manual-memory` | Custom memory management (singletons) |
| `wrap-attach:ClassName` | Attach namespace function as static method |

## Gotchas

- `getDriftTimeUnitAsString()` returns `'<NONE>'` not empty string
- Mobilogram has no MetaInfoInterface (no `getKeys()`, `getMetaValue()`)
- `# wrap-doc:` must be properly indented (strict whitespace parsing)
- Test from `/tmp` to avoid path shadowing with source `pyopenms/`
- Clear libclang cache after header changes: `rm -rf OpenMS-build/libclang_cache`

## Testing

```bash
# All tests (from /tmp to avoid import shadowing)
cd /tmp && PYTHONPATH=.../OpenMS-build/pyOpenMS python3 -m pytest .../src/pyOpenMS/tests/ -v

# Legacy tests only
python3 -m pytest src/pyOpenMS/tests/unittests/ -v

# Nanobind-specific tests
python3 -m pytest src/pyOpenMS/tests/nanobind/ -v
```
