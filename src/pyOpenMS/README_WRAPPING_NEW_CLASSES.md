# How to Wrap a New C++ Class in pyOpenMS

This guide explains how to expose an OpenMS C++ class to Python via nanobind.

## Overview

pyOpenMS uses hand-maintained nanobind C++ binding files in `bindings/bind_<domain>.cpp`.
Each file is compiled into a separate shared library module (`_pyopenms_<domain>.so`).
All modules share types via `NB_DOMAIN "pyopenms"`.

## Step 1: Choose the Binding File

Pick `bindings/bind_<domain>.cpp` based on the C++ header path:

| Header path | Binding file |
|---|---|
| `KERNEL/` (general) | `bind_kernel.cpp` |
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
| Everything else | `bind_misc.cpp` |

## Step 2: Add Include and Class Binding

Add the C++ header include at the top of the file, then add the class binding
inside the `NB_MODULE(...)` block:

```cpp
#include <OpenMS/MODULE/MyClass.h>

// Inside NB_MODULE:

// -----------------------------------------------------------------------
// MyClass
// -----------------------------------------------------------------------
auto myclass = nb::class_<OpenMS::MyClass>(m, "MyClass", "Short description")
    .def(nb::init<>())
    .def(nb::init<const OpenMS::MyClass &>())
    .def("__copy__", [](const OpenMS::MyClass& self) { return OpenMS::MyClass(self); })
    .def("__deepcopy__", [](const OpenMS::MyClass& self, nb::dict) {
        return OpenMS::MyClass(self);
    }, "memo"_a)
    .def("getSomething", [](const OpenMS::MyClass& self) {
        return self.getSomething();
    }, "Returns something")
    .def("setSomething", [](OpenMS::MyClass& self, int value) {
        self.setSomething(value);
    }, "value"_a, "Sets something")
    ;
```

### Method binding patterns

```cpp
// Simple getter (const reference return → copy by default)
.def("getName", [](const OpenMS::MyClass& self) { return self.getName(); })

// Getter returning a reference (keep alive with parent)
.def("getSettings", [](const OpenMS::MyClass& self) -> const OpenMS::Settings& {
    return self.getSettings();
}, nb::rv_policy::reference_internal)

// Setter with named argument
.def("setName", [](OpenMS::MyClass& self, const OpenMS::String& name) {
    self.setName(name);
}, "name"_a, "Sets the name")

// Method with default argument
.def("sort", [](OpenMS::MyClass& self, bool reverse) {
    self.sort(reverse);
}, "reverse"_a = false)

// Static method
.def_static("create", []() { return OpenMS::MyClass::create(); })

// Read-write property (direct member access)
.def_rw("filename", &OpenMS::MyClass::filename)

// Read-only property
.def_ro("size", &OpenMS::MyClass::size)

// Operators
.def(nb::self == nb::self)
.def(nb::self != nb::self)

// Hash support
.def("__hash__", [](const OpenMS::MyClass& self) {
    return std::hash<OpenMS::MyClass>{}(self);
})
```

## Step 3: Handle Inheritance

### Simple inheritance (safe)

If both the base and derived class have **matching vtable status** (both have virtual
destructors, or neither does), declare the base in nanobind directly:

```cpp
// OK: Both Acquisition and MetaInfoInterface have non-virtual destructors
nb::class_<OpenMS::Acquisition, OpenMS::MetaInfoInterface>(m, "Acquisition", ...)

// OK: Both Precursor and CVTermList have virtual destructors
nb::class_<OpenMS::Precursor, OpenMS::CVTermList>(m, "Precursor", ...)
```

### Virtual destructor mismatch — CRITICAL

If the **derived class has a virtual destructor** but the base class (`MetaInfoInterface`)
does **NOT**, you **must NOT** declare the base in nanobind. Instead, use the helper
template from `binding_utils.h`:

```cpp
// WRONG — will segfault!
// PeptideHit has virtual ~PeptideHit() but MetaInfoInterface has no vtable
nb::class_<OpenMS::PeptideHit, OpenMS::MetaInfoInterface>(m, "PeptideHit", ...)

// CORRECT — use helper template instead of nanobind base
auto cls = nb::class_<OpenMS::PeptideHit>(m, "PeptideHit", ...)
    .def(...)
    ;
def_MetaInfoInterface<OpenMS::PeptideHit>(cls);
```

**Why this happens:** When a derived class introduces a vtable pointer (via virtual
destructor) but the base has none, the vtable pointer shifts the base subobject by 8
bytes. nanobind's internal base class pointer cast uses the wrong offset, causing
segfaults in `setMetaValue`, `isMetaEmpty`, `clearMetaInfo`, etc.

**How to check:** Look for `virtual ~ClassName()` or `~ClassName() override` in the
C++ header. If present and the base (`MetaInfoInterface`) is non-virtual, use the
helper template pattern.

### Classes currently using the helper pattern

These classes have virtual destructors and use `def_MetaInfoInterface` instead of
declaring MetaInfoInterface as a nanobind base:

- `CVTermList`, `ChromatogramSettings`
- `PeptideHit`, `PeptideIdentification`, `ProteinIdentification`
- `MSSpectrum`, `MSChromatogram`, `ExperimentalSettings`
- Various other classes with multiple inheritance

### Available helper templates

All defined in `binding_utils.h`:

| Helper | Methods it adds |
|---|---|
| `def_MetaInfoInterface<T>(cls)` | `getMetaValue`, `setMetaValue`, `removeMetaValue`, `metaRegistry`, `getKeys`, `isMetaEmpty`, `clearMetaInfo`, `metaValueExists` |
| `def_UniqueIdInterface<T>(cls)` | `getUniqueId`, `setUniqueId` (3 overloads), `clearUniqueId`, `hasValidUniqueId`, `hasInvalidUniqueId`, `ensureUniqueId`, `isValid` |
| `def_CVTermList<T>(cls)` | `addCVTerm`, `getCVTerms`, `setCVTerms`, `hasCVTerm`, `replaceCVTerm`, `replaceCVTerms`, `consumeCVTerms`, `empty` |
| `def_DefaultParamHandler<T>(cls)` | `setParameters`, `getParameters`, `getDefaults`, `getName`, `setName`, `getSubsections` |
| `def_ProgressLogger<T>(cls)` | `setLogType`, `getLogType`, `startProgress`, `setProgress`, `endProgress`, `nextProgress` |
| `def_DocumentIdentifier<T>(cls)` | `setIdentifier`, `getIdentifier`, `setLoadedFilePath`, `getLoadedFilePath`, `setLoadedFileType`, `getLoadedFileType` |

## Step 4: Bind Enums

```cpp
// Basic enum
nb::enum_<OpenMS::MyEnum>(m, "MyEnum")
    .value("VALUE_A", OpenMS::MyEnum::VALUE_A)
    .value("VALUE_B", OpenMS::MyEnum::VALUE_B)
    ;

// Enum that supports int() conversion — add nb::is_arithmetic()
nb::enum_<OpenMS::DriftTimeUnit>(m, "DriftTimeUnit", nb::is_arithmetic())
    .value("NONE", OpenMS::DriftTimeUnit::NONE)
    .value("MILLISECOND", OpenMS::DriftTimeUnit::MILLISECOND)
    ;

// Nested enum (scoped under a class — pass class variable instead of m)
nb::enum_<OpenMS::MyClass::InnerEnum>(myclass, "InnerEnum")
    .value("A", OpenMS::MyClass::InnerEnum::A)
    ;
```

Use `nb::is_arithmetic()` when the enum should be convertible to/from `int` in Python.

## Step 5: Return NumPy Arrays

For C++ methods that fill output vectors by reference, wrap them to return numpy arrays:

```cpp
.def("get2DPeakData", [](const OpenMS::MSExperiment& self,
                          double min_rt, double max_rt,
                          double min_mz, double max_mz,
                          size_t ms_level) {
    std::vector<float> rt, mz, intensity;
    self.get2DPeakData(min_rt, max_rt, min_mz, max_mz, ms_level, rt, mz, intensity);
    const size_t n = rt.size();

    // Allocate owning arrays
    std::unique_ptr<float[]> rt_uptr(new float[n]);
    std::copy(rt.begin(), rt.end(), rt_uptr.get());

    // Create capsule to manage lifetime
    float* rt_data = rt_uptr.release();
    nb::capsule rt_owner(rt_data, [](void* p) noexcept { delete[] static_cast<float*>(p); });

    // Return as numpy array
    return nb::ndarray<nb::numpy, float, nb::ndim<1>>(rt_data, {n}, rt_owner);
}, "min_rt"_a, "max_rt"_a, "min_mz"_a, "max_mz"_a, "ms_level"_a)
```

**Key pattern:** allocate with `unique_ptr`, release into raw pointer, wrap in `nb::capsule`
for Python-managed lifetime, then construct `nb::ndarray` with the capsule as owner.

## Step 6: Container Protocols

For classes that behave like containers:

```cpp
// Iteration
.def("__iter__", [](OpenMS::MyContainer& self) {
    return nb::make_iterator<nb::rv_policy::reference_internal>(
        nb::handle(), "MyContainer_iter", self.begin(), self.end());
})

// Length
.def("__len__", [](const OpenMS::MyContainer& self) { return self.size(); })

// Indexing
.def("__getitem__", [](OpenMS::MyContainer& self, size_t i) -> OpenMS::Element& {
    if (i >= self.size()) throw nb::index_error();
    return self[i];
}, nb::rv_policy::reference_internal)

.def("__setitem__", [](OpenMS::MyContainer& self, size_t i, const OpenMS::Element& val) {
    if (i >= self.size()) throw nb::index_error();
    self[i] = val;
}, "i"_a, "val"_a)
```

## Step 7: Build and Test

```bash
# Rebuild
cmake --build OpenMS-build --target pyopenms -j$(nproc)

# Test (from /tmp to avoid import shadowing)
cd /tmp && PYTHONPATH=.../OpenMS-build/pyOpenMS python3 -m pytest .../src/pyOpenMS/tests/ -v

# Quick smoke test
cd /tmp && PYTHONPATH=.../OpenMS-build/pyOpenMS python3 -c "
import pyopenms as oms
obj = oms.MyClass()
print(dir(obj))
"
```

## Adding a New Domain Module

If you need to create an entirely new binding file (rare):

1. Create `bindings/bind_<domain>.cpp` with:
   ```cpp
   #include "all_casters.h"
   #include <nanobind/nanobind.h>
   #include "binding_utils.h"

   namespace nb = nanobind;
   using namespace nb::literals;

   NB_MODULE(_pyopenms_<domain>, m) {
       m.doc() = "pyOpenMS <domain> bindings";
       // ... class bindings ...
   }
   ```
2. Add `<domain>` to `PYOPENMS_DOMAINS` in `src/pyOpenMS/CMakeLists.txt`
3. If load order matters, add `_pyopenms_<domain>` to `_priority_modules` in `pyopenms/__init__.py`
4. Re-run cmake: `cd OpenMS-build && cmake ..`
5. Build: `cmake --build OpenMS-build --target pyopenms -j$(nproc)`

## Type Casters

Custom type casters in `bindings/type_casters/` handle automatic C++ ↔ Python conversion:

| C++ type | Python type |
|---|---|
| `OpenMS::String` | `str` |
| `OpenMS::DataValue` | `int` / `float` / `str` / `list` |
| `OpenMS::ParamValue` | `int` / `float` / `str` / `list` |
| `DPosition<1>` | `float` |
| `DPosition<2>` | `tuple` |
| `std::string` | `str` (also accepts `bytes`) |

## Adding Python-Only Methods (Addons)

For methods better implemented in Python (e.g., DataFrame export), add an addon file
in `pyopenms/addons/`:

```python
from pyopenms.addons import addon

@addon("MyClass")
def my_python_method(self, arg):
    """This method is injected into MyClass at import time."""
    return self.someMethod(arg) + 1
```

## Cross-Module Type Sharing

All domain modules use `NB_DOMAIN "pyopenms"` so types defined in one module
(e.g., `MSSpectrum` in `_pyopenms_spectrum`) are available in another module
(e.g., `MzMLFile` in `_pyopenms_format`). No special imports needed.

## Common Pitfalls

1. **Segfault on MetaInfoInterface methods**: Virtual destructor mismatch. See Step 3.
2. **Build doesn't pick up new files**: Re-run cmake after adding to `PYOPENMS_DOMAINS`.
3. **Import shadowing**: Test from `/tmp` with `PYTHONPATH`, not from the source tree.
4. **Enum not convertible to int**: Add `nb::is_arithmetic()` to the enum definition.
5. **`setUniqueId(int)` missing**: Use `def_UniqueIdInterface<T>()` helper which includes
   all overloads, not just the no-arg version.
6. **Object lifetime issues**: Use `nb::rv_policy::reference_internal` for methods
   returning references to internal data.
7. **"incompatible function arguments"**: Check parameter types match exactly. Use lambdas
   for explicit conversion.
8. **"nanobind: type not found"**: The class might be in a different domain module. Check
   `NB_DOMAIN` is set.
9. **Linker errors about undefined symbols**: Make sure the `#include` is added and the
   class is in the OpenMS library.
