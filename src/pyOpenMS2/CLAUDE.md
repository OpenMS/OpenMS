# pyOpenMS2 Development

pyOpenMS2 is the nanobind-based replacement for pyOpenMS (which uses autowrap/Cython).

## Architecture

- **Generator** (`generator/`) parses existing `.pxd` files from `../pyOpenMS/pxds/` and emits nanobind C++ code
- **Type casters** (`bindings/type_casters/`) handle C++ ↔ Python type conversion
- **Addon system** (`pyopenms/addons/`) injects Python methods into wrapped classes at import time
- **Generated bindings** (`bindings/`) contain the nanobind C++ wrapper code

## Quick Commands

```bash
# Build
cmake --build OpenMS-build --target pyopenms2 -j$(nproc)

# Test
PYTHONPATH=OpenMS-build/pyOpenMS2 python3 -m pytest src/pyOpenMS2/tests/ -v

# Run generator only (dry-run)
cd src/pyOpenMS2 && python -m generator --pxd-dir ../pyOpenMS/pxds --addons-dir ../pyOpenMS/addons --output-dir generated --dry-run --openms-include-dir ../../src/openms/include ../../OpenMS-build/src/openms/include

# Generate bindings
cd src/pyOpenMS2 && python -m generator --pxd-dir ../pyOpenMS/pxds --addons-dir ../pyOpenMS/addons --output-dir bindings/generated --num-modules 8 --all-classes --openms-include-dir ../../src/openms/include ../../OpenMS-build/src/openms/include
```

## Build Gotchas

- **Generator needs both include dirs**: `--openms-include-dir` must include both the source tree (`src/openms/include`) and the build tree (`OpenMS-build/src/openms/include`). The build tree contains generated headers like `OpenMS/config.h`. Without it, libclang fails to parse most headers (42 classes instead of 772+).
- **CMake handles this automatically** via `INTERFACE_INCLUDE_DIRECTORIES` on the OpenMS target. The issue only arises when running the generator manually.
- **Feature parity testing**: Run old pyOpenMS tests against the pyOpenMS2 build with `PYTHONPATH=OpenMS-build/pyOpenMS2-build python3 -m pytest src/pyOpenMS/tests/unittests/ -v`. This requires a shim that redirects `import pyopenms` to pyOpenMS2.

## Generator Components

| Module | Purpose |
|--------|---------|
| `pxd_parser.py` | Parses .pxd files, extracts class/method declarations and wrap- directives |
| `cpp_parser.py` | Uses libclang to parse C++ headers for accurate type information |
| `type_registry.py` | Maps C++ types to Python types, tracks custom casters needed |
| `nanobind_emitter_v2.py` | Generates nanobind C++ binding code using libclang-merged info |
| `addon_processor.py` | Classifies addons as pure Python vs C++ required |

## Type Casters

Custom type casters in `bindings/type_casters/`:

| Header | Types Handled |
|--------|---------------|
| `string.h` | `OpenMS::String` ↔ Python `str` |
| `dposition.h` | `DPosition<1>` ↔ `float`, `DPosition<2>` ↔ `tuple` |
| `datavalue.h` | `DataValue`, `ParamValue` ↔ Python variants |
| `stl_containers.h` | `vector<String>`, `set<String>`, `map<String, V>` |

## Addon System

Addons are pure Python methods injected into C++ classes at import time:

```python
# pyopenms/addons/msspectrum.py
from . import addon

@addon("MSSpectrum")
def __repr__(self):
    return f"MSSpectrum(n_peaks={len(self)}, rt={self.getRT():.2f}s)"
```

Performance-critical methods (`get_peaks`, `set_peaks`) are implemented in C++.

## Key Differences from pyOpenMS (autowrap)

| Aspect | pyOpenMS (autowrap) | pyOpenMS2 (nanobind) |
|--------|---------------------|----------------------|
| Build time | Slow (Cython compilation) | Fast (direct C++) |
| Binary size | Large | Small (5-10x smaller) |
| Addons | Cython .pyx files | Pure Python + C++ lambdas |
| Type conversion | autowrap converters | nanobind type casters |
| Multi-module | 8 .so files | 8 .so files (same) |

## Module Split

Bindings are split across 8 modules for parallel compilation and memory management:
- `_pyopenms2_1.so` through `_pyopenms2_8.so`
- `_pyopenms2.so` (main entry that imports all)

Classes are distributed using hash-based assignment for even split.

## NB_DOMAIN

All modules use `NB_DOMAIN "pyopenms"` to share type information across the split modules. This is critical for cross-module type visibility.

## Testing

```bash
# Run all tests
pytest src/pyOpenMS2/tests/ -v

# Run specific test file
pytest src/pyOpenMS2/tests/test_msspectrum.py -v

# Run with coverage
pytest src/pyOpenMS2/tests/ --cov=pyopenms --cov-report=html
```

## Migration Notes

When porting features from pyOpenMS:
1. Check if the method uses `cdef` or direct `.inst.get()` access
2. If yes → implement in C++ (lambda in binding or type caster)
3. If no → convert to pure Python addon

## Common Patterns

```cpp
// Simple method binding
.def("getRT", &OpenMS::MSSpectrum::getRT, "Returns retention time")

// Method with overloads
.def("findNearest", nb::overload_cast<double>(&OpenMS::MSSpectrum::findNearest, nb::const_))

// Lambda for custom behavior
.def("get_peaks", [](const OpenMS::MSSpectrum& self) {
    // ... numpy array creation
})

// Iterator support
.def("__iter__", [](OpenMS::MSSpectrum& self) {
    return nb::make_iterator(self.begin(), self.end());
}, nb::keep_alive<0, 1>())
```

## Wrap- Directive Reference

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
