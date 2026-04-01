# pyOpenMS Development

## Architecture

- **Bindings** (`bindings/`) — hand-maintained nanobind C++ binding files (13 domain modules + main)
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

### 1. Choose the binding file

Pick `bindings/bind_<domain>.cpp` based on the C++ header path:

| Header path | Binding file |
|---|---|
| `KERNEL/` | `bind_kernel.cpp` |
| `KERNEL/MSSpectrum.h` | `bind_spectrum.cpp` |
| `KERNEL/MSChromatogram.h` | `bind_chromatogram.cpp` |
| `KERNEL/MSExperiment.h` | `bind_experiment.cpp` |
| `FORMAT/` | `bind_format.cpp` |
| `ANALYSIS/` | `bind_analysis.cpp` |
| `CHEMISTRY/` | `bind_chemistry.cpp` |
| `METADATA/` | `bind_metadata.cpp` |
| `PROCESSING/` | `bind_processing.cpp` |
| `FEATUREFINDER/` | `bind_featurefinder.cpp` |
| `DATASTRUCTURES/`, `MATH/`, `CONCEPT/` | `bind_datastructures.cpp` |
| `ML/` | `bind_ml.cpp` |
| `IONMOBILITY/` enums (`DriftTimeUnit`, `IMFormat`, `IMPeakType`) | `bind_kernel.cpp` |
| Everything else (incl. `IONMOBILITY/IMTypes`) | `bind_misc.cpp` |

### 2. Add the `#include` and class binding

```cpp
#include <OpenMS/MODULE/MyClass.h>

// Inside NB_MODULE:
auto myclass = nb::class_<OpenMS::MyClass>(m, "MyClass", "docstring")
    .def(nb::init<>())
    .def(nb::init<const OpenMS::MyClass &>())
    .def("__copy__", [](const OpenMS::MyClass& self) { return OpenMS::MyClass(self); })
    .def("__deepcopy__", [](const OpenMS::MyClass& self, nb::dict) { return OpenMS::MyClass(self); }, "memo"_a)
    .def("someMethod", [](const OpenMS::MyClass& self) { return self.someMethod(); })
    ;
```

### 3. Inheritance rules

**Simple base (no vtable issues):** If the base class has NO virtual destructor, declare it as nanobind base:
```cpp
// OK: MetaInfoInterface has no virtual destructor, Acquisition has no virtual destructor
nb::class_<OpenMS::Acquisition, OpenMS::MetaInfoInterface>(m, "Acquisition", ...)
```

**Virtual destructor mismatch — CRITICAL:** If the derived class has a virtual destructor but the base class (`MetaInfoInterface`) does NOT, you **must NOT** declare the base in nanobind. Instead, use the helper template from `binding_utils.h`:
```cpp
// WRONG — will segfault! PeptideHit has virtual ~PeptideHit() but MetaInfoInterface has no vtable
nb::class_<OpenMS::PeptideHit, OpenMS::MetaInfoInterface>(m, "PeptideHit", ...)

// CORRECT — use helper template
auto cls = nb::class_<OpenMS::PeptideHit>(m, "PeptideHit", ...)
    .def(...)
    ;
def_MetaInfoInterface<OpenMS::PeptideHit>(cls);
```

**Why:** When a derived class introduces a vtable pointer but the base has none, nanobind computes incorrect pointer offsets during base class casts, causing segfaults in `setMetaValue`, `isMetaEmpty`, etc.

**How to check:** Look for `virtual ~ClassName()` or `~ClassName() override` in the C++ header. If present and the base (MetaInfoInterface) is non-virtual, use the helper template pattern.

**Classes currently using the helper pattern (no nanobind base):**
- CVTermList, ChromatogramSettings, PeptideHit, PeptideIdentification, ProteinIdentification
- MSSpectrum, MSChromatogram, ExperimentalSettings, and others with multiple inheritance

**Available helper templates** (in `binding_utils.h`):
- `def_MetaInfoInterface<T>(cls)` — binds getMetaValue, setMetaValue, isMetaEmpty, etc.
- `def_UniqueIdInterface<T>(cls)` — binds getUniqueId, setUniqueId, etc.
- `def_CVTermList<T>(cls)` — binds addCVTerm, getCVTerms, hasCVTerm, etc.
- `def_DefaultParamHandler<T>(cls)` — binds setParameters, getParameters, getDefaults, etc.
- `def_ProgressLogger<T>(cls)` — binds setLogType, getLogType, startProgress, etc.
- `def_DocumentIdentifier<T>(cls)` — binds setIdentifier, getLoadedFilePath, etc.

### 4. Enums

```cpp
// Basic enum
nb::enum_<OpenMS::MyEnum>(m, "MyEnum")
    .value("VALUE1", OpenMS::MyEnum::VALUE1)
    .value("VALUE2", OpenMS::MyEnum::VALUE2)
    ;

// Enum supporting int() conversion — add nb::is_arithmetic()
nb::enum_<OpenMS::DriftTimeUnit>(m, "DriftTimeUnit", nb::is_arithmetic())
    .value("NONE", OpenMS::DriftTimeUnit::NONE)
    ...

// Nested enum (scoped under a class)
nb::enum_<OpenMS::MyClass::InnerEnum>(myclass, "InnerEnum")
    .value("A", OpenMS::MyClass::InnerEnum::A)
    ;
```

### 5. Methods returning numpy arrays

For C++ methods that fill output vectors by reference, wrap them to return numpy arrays:
```cpp
.def("get2DPeakData", [](const OpenMS::MSExperiment& self, ...) {
    std::vector<float> rt, mz, intensity;
    self.get2DPeakData(..., rt, mz, intensity);
    const size_t n = rt.size();
    // Allocate, copy, create capsule owners, return ndarray tuple
    ...
}, ...)
```

### 6. Build and test

```bash
cmake --build OpenMS-build --target pyopenms -j$(nproc)
cd /tmp && PYTHONPATH=.../OpenMS-build/pyOpenMS python3 -m pytest .../src/pyOpenMS/tests/ -v
```

Each class has a section comment (`// --- ClassName ---`) for navigation.

## Addon Rules

- Pure Python methods in `pyopenms/addons/`
- Use the `@addon("ClassName")` decorator from `pyopenms.addons`
- Keep addons minimal — only add non-auto-generatable methods
- Performance-critical methods should be C++ lambdas in the bindings

## Module Split

13 domain-based nanobind modules for parallel compilation:
`_pyopenms_kernel`, `_pyopenms_spectrum`, `_pyopenms_chromatogram`, `_pyopenms_experiment`,
`_pyopenms_metadata`, `_pyopenms_chemistry`, `_pyopenms_analysis`,
`_pyopenms_featurefinder`, `_pyopenms_format`, `_pyopenms_processing`,
`_pyopenms_datastructures`, `_pyopenms_ml`, `_pyopenms_misc`

Plus `_pyopenms` (main entry) and `_arrow_zerocopy`.

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

# Specific test file
python3 -m pytest src/pyOpenMS/tests/unittests/test_type_casters.py -v
```
