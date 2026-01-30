# pyOpenMS2 — Nanobind-Based Python Bindings for OpenMS

pyOpenMS2 is the next-generation Python binding layer for the [OpenMS](https://www.openms.de/) C++ library for mass spectrometry. It replaces the original pyOpenMS (based on autowrap/Cython) with [nanobind](https://github.com/wentingj/nanobind)-based bindings that are faster to build, produce smaller binaries (5-10x), and are easier to maintain.

## Table of Contents

- [Goals and Scope](#goals-and-scope)
- [How Wrapping Works](#how-wrapping-works)
  - [Overview](#overview)
  - [Stage 1: Parsing .pxd Files](#stage-1-parsing-pxd-files)
  - [Stage 2: C++ Header Parsing with libclang](#stage-2-c-header-parsing-with-libclang)
  - [Stage 3: Type Resolution](#stage-3-type-resolution)
  - [Stage 4: Code Emission](#stage-4-code-emission)
  - [Stage 5: Module Splitting](#stage-5-module-splitting)
- [Supported Features](#supported-features)
  - [Classes and Methods](#classes-and-methods)
  - [Operators and Special Methods](#operators-and-special-methods)
  - [Enums](#enums)
  - [Templates](#templates)
  - [Type Conversions](#type-conversions)
  - [Memory Management](#memory-management)
  - [Data Export](#data-export)
- [Addon System](#addon-system)
- [Type Casters](#type-casters)
- [Wrap-Directive Reference](#wrap-directive-reference)
- [Multi-Module Architecture](#multi-module-architecture)
- [Stub Generation (PEP 561)](#stub-generation-pep-561)
- [Feature Parity with pyOpenMS](#feature-parity-with-pyopenms)
- [Future Direction](#future-direction)
- [Building](#building)
- [Testing](#testing)

---

## Goals and Scope

The primary goal of the current code generation pipeline is to achieve **complete feature parity** with the existing pyOpenMS Cython bindings. The emitters and generator are designed so that all former pyOpenMS tests pass against pyOpenMS2 without modification. This means:

- Every class, method, enum, and operator exposed by pyOpenMS must also be exposed by pyOpenMS2.
- The Python-level API (class names, method signatures, return types) must be identical.
- Addon methods (`get_peaks`, `set_peaks`, `to_df`, etc.) must behave the same way.
- File I/O helper signatures and behavior must match.

The `.pxd` files from the original pyOpenMS serve as the **allowlist and contract**: only classes and methods declared there are wrapped, and the directives within them control renaming, skipping, template instantiation, and documentation.

## How Wrapping Works

### Overview

pyOpenMS2 uses a multi-stage code generation pipeline that reads the existing pyOpenMS `.pxd` declarations, parses the actual C++ headers with libclang for accurate type information, and emits nanobind C++ binding code. The pipeline runs at build time and produces a set of `.cpp` files that are compiled into shared libraries.

```
┌─────────────────┐     ┌──────────────────┐     ┌──────────────────┐
│  .pxd files     │     │  C++ headers     │     │  .pyx addons     │
│  (allowlist)    │     │  (source of      │     │  (legacy Cython  │
│                 │     │   truth)         │     │   addons)        │
└────────┬────────┘     └────────┬─────────┘     └────────┬─────────┘
         │                       │                         │
         ▼                       ▼                         ▼
   ┌───────────┐          ┌────────────┐          ┌──────────────────┐
   │ PxdParser │          │ CppHeader  │          │ AddonProcessor   │
   │           │          │ Parser     │          │ (classify as     │
   │           │          │ (libclang) │          │  pure-Python or  │
   │           │          │            │          │  C++-required)   │
   └─────┬─────┘          └─────┬──────┘          └────────┬─────────┘
         │                       │                         │
         └───────────┬───────────┘                         │
                     ▼                                     │
              ┌─────────────┐                              │
              │ merge_with_ │                              │
              │ pxd()       │◄─────────────────────────────┘
              └──────┬──────┘
                     │
                     ▼
              ┌─────────────┐
              │ TypeRegistry │
              │ (resolve     │
              │  all types)  │
              └──────┬──────┘
                     │
                     ▼
              ┌──────────────┐
              │ Nanobind     │
              │ Emitter      │──► module_1.cpp ... module_8.cpp
              │ (5500 lines) │──► main_module.cpp
              └──────────────┘
```

### Stage 1: Parsing .pxd Files

The `PxdParser` (`generator/pxd_parser.py`) reads Cython `.pxd` files from `src/pyOpenMS/pxds/`. These files are **not** used to generate Cython code — they serve as a structured allowlist that declares:

- Which C++ classes to wrap and their C++ namespaces
- Which methods and constructors to expose
- Which methods to skip (`wrap-ignore`)
- How to rename methods (`wrap-as:`)
- Template instantiations (`wrap-instances:`)
- Inheritance relationships (`wrap-inherits:`)
- Documentation strings (`wrap-doc:`)

The parser extracts class declarations, method signatures, enum definitions, and all `wrap-` directives into an intermediate data structure.

### Stage 2: C++ Header Parsing with libclang

The `CppHeaderParser` (`generator/cpp_parser.py`) uses **libclang** (the Clang compiler's C API) to parse the actual C++ headers in `src/openms/include/`. This provides:

- Accurate parameter types, including const qualifiers and reference/pointer distinctions
- Default argument values
- Accurate return types (the .pxd files sometimes simplify these)
- Detection of pure virtual methods (abstract classes are auto-skipped)
- Detection of incomplete/forward-declared types (methods using them are auto-blocked)
- Scoped vs. unscoped enum detection
- Template parameter information

The libclang-parsed information is **merged** with the .pxd-parsed information via `merge_with_pxd()`. The .pxd acts as the filter (only declared items are wrapped), while libclang provides the ground truth for types and signatures.

### Stage 3: Type Resolution

The `TypeRegistry` (`generator/type_registry.py`) maintains a complete mapping from C++ types to their Python equivalents and nanobind binding strategies. Types are categorized as:

| Category | Examples | Strategy |
|----------|----------|----------|
| `PRIMITIVE` | `int`, `double`, `bool`, `size_t` | Direct nanobind conversion |
| `STRING` | `OpenMS::String` | Custom type caster |
| `CONTAINER` | `vector<T>`, `map<K,V>`, `set<T>` | nanobind STL casters + custom for String containers |
| `CUSTOM_CASTER` | `DPosition<2>`, `DataValue` | Custom type caster headers |
| `CLASS` | `MSSpectrum`, `Feature` | Bound as `nb::class_<T>` |
| `ENUM` | `MSLevel`, `PeakType` | Bound as `nb::enum_<T>` |
| `OPAQUE` | `shared_ptr<T>`, forward-declared types | Passed by reference, not converted |

The registry also tracks which custom caster headers must be included for each generated module.

### Stage 4: Code Emission

The `NanobindEmitter` (`generator/nanobind_emitter.py`, ~5500 lines) is the core of the generator. It takes the merged class/method data and type registry and emits valid nanobind C++ code. For each class it generates:

- `nb::class_<ClassName>` declaration with proper base classes
- Constructor bindings via `nb::init<>()`
- Method bindings with proper overload resolution (`nb::overload_cast`)
- Const/non-const method disambiguation
- Operator overloads (`__eq__`, `__lt__`, `__getitem__`, `__len__`, etc.)
- Enum definitions attached to classes or at module level
- Public data member access via `.def_rw()` / `.def_ro()`
- Iterator support via `nb::make_iterator()` with `nb::keep_alive`

For methods that are too complex for automatic generation (e.g., methods requiring numpy array creation, custom memory layouts, or type-punning), the emitter contains a `SPECIAL_METHODS` dictionary with ~500+ handwritten C++ lambda implementations. These cover performance-critical methods like `get_peaks()`, `set_peaks()`, file I/O helpers, and methods with complex parameter conversions.

### Stage 5: Module Splitting

The generated bindings are split across **8 submodules** (`module_1.cpp` through `module_8.cpp`) plus a `main_module.cpp` entry point. Classes are distributed using hash-based assignment for even load balancing. This splitting is necessary because:

- A single translation unit would exceed compiler memory limits
- Parallel compilation of 8 modules is significantly faster
- Incremental rebuilds only recompile affected modules

All modules share type information via `NB_DOMAIN "pyopenms"`, which is nanobind's mechanism for cross-module type visibility. Without this, a `MSSpectrum` created in module 1 would not be recognized by functions in module 5.

## Supported Features

### Classes and Methods

- **770+ C++ classes** wrapped (dynamically detected, filtered by .pxd allowlist)
- Instance methods (const and non-const)
- Static methods
- Multiple constructors
- Method overloading with automatic disambiguation
- Default argument values preserved from C++ headers
- Public data members (read-write and read-only)

### Operators and Special Methods

| Python Method | C++ Source |
|---------------|------------|
| `__eq__`, `__ne__`, `__lt__`, `__le__`, `__gt__`, `__ge__` | Comparison operators |
| `__add__`, `__sub__`, `__mul__`, `__truediv__` | Arithmetic operators |
| `__iadd__`, `__isub__`, `__imul__` | In-place arithmetic |
| `__getitem__`, `__setitem__` | `operator[]` with bounds checking |
| `__len__` | `.size()` |
| `__iter__`, `__next__` | `begin()`/`end()` iterators |
| `__contains__` | Container search |
| `__hash__` | Enabled via `wrap-hash:` directive |
| `__repr__`, `__str__` | Addon-provided or auto-generated |
| `__bool__` | Boolean conversion |
| `__neg__`, `__invert__` | Unary operators |

### Enums

Both traditional C enums and C++11 scoped enums (`enum class`) are supported. The emitter detects scoped vs. unscoped enums from C++ headers and generates the appropriate nanobind enum binding. Enums can be:

- Defined at module level
- Attached to a class via `wrap-attach:` directive
- Accessed via both the class and module level for backward compatibility

### Templates

Template classes are instantiated via the `wrap-instances:` directive in .pxd files. For example:

```cython
cdef cppclass MRMTransitionGroup[SpectrumT, TransitionT]:
    # wrap-instances:
    #   MRMTransitionGroupCP := MRMTransitionGroup[MSChromatogram, ReactionMonitoringTransition]
```

The generator maps Cython type names to C++ types and generates a separate Python class for each instantiation. Template mappings are maintained in `_TEMPLATE_PXD_TO_CPP` within the emitter.

### Type Conversions

Custom type casters handle automatic bidirectional conversion between C++ and Python types:

| C++ Type | Python Type | Notes |
|----------|-------------|-------|
| `OpenMS::String` | `str` | Full Unicode support, bytes fallback |
| `DPosition<1>` | `float` | Scalar position |
| `DPosition<2>` | `tuple[float, float]` | 2D position |
| `DataValue` | `None`, `int`, `float`, `str`, `list` | Variant-like union |
| `ParamValue` | Same as DataValue | Same caster |
| `std::string` | `str` or `bytes` | Accepts both |
| `vector<T>` | `list[T]` | Standard nanobind STL caster |
| `map<K, V>` | `dict[K, V]` | Standard nanobind STL caster |
| `set<T>` | `set[T]` | Standard nanobind STL caster |
| `vector<String>` | `list[str]` | Custom caster for String elements |
| `shared_ptr<T>` | Managed reference | Automatic reference counting |

### Memory Management

- **Automatic reference counting** via nanobind for all wrapped classes
- **`shared_ptr<T>`** and **`unique_ptr<T>`** holders detected and handled automatically
- **`keep_alive` policies** on iterators to prevent dangling references
- **Singleton pattern** support: classes with `wrap-manual-memory` are wrapped as callables that return the singleton instance (e.g., `ElementDB()` calls `getInstance()`)
- **Abstract class detection**: classes with pure virtual methods are automatically excluded from binding

### Data Export

pyOpenMS2 provides multiple data export pathways:

- **NumPy arrays**: `get_peaks()` returns `(mz_array, intensity_array)` as numpy arrays (zero-copy where possible)
- **Pandas DataFrames**: `to_df()` and `get_df()` methods on spectra, chromatograms, features, and maps
- **Dictionary export**: `get_data_dict()` for flexible column selection
- **Apache Arrow tables** (optional, requires `WITH_PARQUET=ON`): `spectra_to_arrow()` and `chromatograms_to_arrow()` via the Arrow C Data Interface for true zero-copy export. Supports "long" (one row per peak) and "semi_wide" (grouped) formats with filtering by RT, m/z, and MS level.

## Addon System

The addon system allows pure Python methods to be injected into C++ classes at import time. This provides a clean separation between the C++ binding layer and Python convenience methods.

```python
# pyopenms/addons/msspectrum.py
from . import addon

@addon("MSSpectrum")
def __repr__(self):
    return f"MSSpectrum(n_peaks={len(self)}, rt={self.getRT():.2f}s)"

@addon("MSSpectrum")
def get_tic(self):
    """Return the total ion current (sum of intensities)."""
    import numpy as np
    _, intensities = self.get_peaks()
    return float(np.sum(intensities))
```

Addons are loaded via `apply_addons(globals())` in `pyopenms/__init__.py`. The `AddonProcessor` classifies legacy Cython `.pyx` addon files into:

- **`PURE_PYTHON`** — Can be converted to a Python addon directly
- **`CPP_PERFORMANCE`** — Should be implemented in C++ for performance (e.g., `get_peaks`)
- **`CPP_REQUIRED`** — Must be in C++ (uses `cdef`, memoryviews, `.inst.get()`, etc.)

Performance-critical methods like `get_peaks`, `set_peaks`, `get_data_dict`, and `to_df` are implemented as C++ lambdas in the emitter's `SPECIAL_METHODS` dictionary rather than as Python addons, to avoid per-element Python overhead on large arrays.

## Type Casters

Custom nanobind type casters live in `bindings/type_casters/` and are automatically included in generated modules:

| Header | Types | Description |
|--------|-------|-------------|
| `openms_string_caster.h` | `OpenMS::String` | Bidirectional str conversion with bytes fallback |
| `openms_dposition_caster.h` | `DPosition<1>`, `DPosition<2>` | Scalar float or (float, float) tuple |
| `openms_datavalue_caster.h` | `DataValue`, `ParamValue` | Variant: None/int/float/str/list |
| `openms_stl_caster.h` | `vector<String>`, `map<String, V>`, `set<String>` | STL containers with String conversion |
| `std_string_bytes_caster.h` | `std::string` | Accept both str and bytes input |
| `all_casters.h` | (master include) | Includes all of the above plus nanobind STL support |

The type registry scans these casters during generation to determine which types should be auto-converted rather than bound as classes.

## Wrap-Directive Reference

Directives are embedded as comments in `.pxd` files and control the code generator:

| Directive | Example | Description |
|-----------|---------|-------------|
| `wrap-doc:` | `# wrap-doc: Returns the retention time` | Documentation string for the method or class |
| `wrap-as:` | `# wrap-as: getName` | Rename the method in the Python API |
| `wrap-ignore` | `# wrap-ignore` | Skip wrapping this method entirely |
| `wrap-instances:` | `# wrap-instances: Peak2D := Peak1D[DPosition[2]]` | Instantiate a template with specific type arguments |
| `wrap-inherits:` | `# wrap-inherits: MetaInfoInterface` | Declare base classes for inheritance |
| `wrap-upper-limit:` | `# wrap-upper-limit: size()` | Add bounds checking to `__getitem__` |
| `wrap-hash:` | `# wrap-hash` | Enable `__hash__` for this class |
| `wrap-manual-memory` | `# wrap-manual-memory` | Singleton pattern: bind as callable returning `getInstance()` |
| `wrap-attach:` | `# wrap-attach: MSSpectrum` | Attach an enum to a class rather than module level |
| `wrap-iter:` | `# wrap-iter: begin(), end(), Peak1D` | Generate `__iter__` using the specified begin/end methods |
| `wrap-const:` | `# wrap-const` | Mark method as const (for overload resolution) |
| `wrap-static:` | `# wrap-static` | Mark method as static |

## Multi-Module Architecture

The binding code is compiled into 9 shared libraries:

```
_pyopenms2.so        ← Main entry point, imports all submodules
_pyopenms2_1.so      ← Classes A-D (hash-distributed)
_pyopenms2_2.so      ← Classes E-H
...
_pyopenms2_8.so      ← Classes W-Z
```

All modules declare `NB_DOMAIN "pyopenms"` so that types defined in one module are visible in all others. The main module imports all submodules at initialization, ensuring the full type registry is populated before user code runs.

The number of modules is configurable via `PY_NUM_MODULES` in CMake (default: 8).

## Stub Generation (PEP 561)

pyOpenMS2 generates `.pyi` type stub files for IDE autocompletion and static type checking:

1. **nanobind.stubgen** generates initial stubs from the compiled `.so` modules
2. **`stubs.patch`** applies pattern fixes (e.g., singleton DB class type annotations)
3. **`fix_stubs.py`** post-processes stubs to fix:
   - Python keywords used as parameter names (appends `_`)
   - Malformed `NDArray` annotations → `numpy.ndarray`
   - Singleton type hints
4. A `py.typed` marker file is created for PEP 561 compliance

This enables full type checking with mypy, pyright, and IDE autocompletion in VS Code, PyCharm, etc.

## Feature Parity with pyOpenMS

The emitters are designed to produce bindings that are **feature-parity complete** with the original Cython-based pyOpenMS. The goal is that all existing pyOpenMS unit tests pass against pyOpenMS2 without modification. This is validated by running:

```bash
PYTHONPATH=OpenMS-build/pyOpenMS2-build python3 -m pytest src/pyOpenMS/tests/unittests/ -v
```

A compatibility shim redirects `import pyopenms` to the pyOpenMS2 package so that existing test code works unchanged.

**What feature parity covers:**

- All ~770+ wrapped classes with identical Python names
- All methods, constructors, and overloads
- All enums (scoped and unscoped) with the same values
- All addon methods (get_peaks, set_peaks, to_df, etc.)
- All operator overloads (__eq__, __getitem__, __iter__, etc.)
- Identical exception behavior and error messages where possible
- Template instantiations producing the same set of concrete classes
- File I/O helper method signatures

**What is intentionally different:**

- Internal implementation (nanobind vs. Cython/autowrap)
- Binary size (5-10x smaller)
- Build time (significantly faster)
- Better error messages from nanobind's type checking
- PEP 561 type stubs (not available in pyOpenMS)

## Future Direction

The current code generation pipeline exists to bootstrap pyOpenMS2 to full feature parity using the existing `.pxd` files as a contract. Once parity is achieved and validated, the development model may shift:

- **Direct editing of generated files**: Instead of regenerating from `.pxd` files, the generated `.cpp` binding files may be checked in and edited directly. This would allow fine-grained control over individual bindings without the indirection of the generator pipeline.
- **Deprecation of the generator**: Once the generated code is stable and manually maintained, the generator and `.pxd` parser would become optional development tools rather than required build steps.
- **Independent evolution**: pyOpenMS2 could expose new C++ API surface not covered by the original `.pxd` files, add Python-specific APIs, or adopt different design patterns without being constrained by the Cython legacy.
- **Removal of SPECIAL_METHODS**: As the generator matures or manual editing takes over, hardcoded C++ blocks can be refined or replaced with cleaner abstractions.

The `.pxd` files will remain useful as documentation of the original API contract even after the generator is no longer the primary development tool.

## Building

pyOpenMS2 is a **separate CMake project**, not a target in the main OpenMS build.

### Prerequisites

- OpenMS built (provides `OpenMSConfig.cmake`)
- Python 3.9+
- NumPy
- libclang (for the generator)
- nanobind v2.4.0 (fetched automatically by CMake)

### Configuration and Build

```bash
# Configure (from repo root, assuming OpenMS built in OpenMS-build/)
cmake -S src/pyOpenMS2 -B OpenMS-build/pyOpenMS2-build \
  -DOpenMS_DIR=OpenMS-build \
  -DCMAKE_BUILD_TYPE=Release

# Build
cmake --build OpenMS-build/pyOpenMS2-build -j$(nproc)
```

### Running the Generator Manually

```bash
cd src/pyOpenMS2
python -m generator \
  --pxd-dir ../pyOpenMS/pxds \
  --addons-dir ../pyOpenMS/addons \
  --output-dir bindings/generated \
  --num-modules 8 \
  --all-classes \
  --openms-include-dir ../../src/openms/include ../../OpenMS-build/src/openms/include
```

Both include directories (source tree and build tree) are required. The build tree contains generated headers like `OpenMS/config.h` — without it, libclang can only parse ~42 classes instead of 770+.

## Testing

```bash
# pyOpenMS2 tests
PYTHONPATH=OpenMS-build/pyOpenMS2-build python3 -m pytest src/pyOpenMS2/tests/ -v

# Feature parity: run original pyOpenMS tests against pyOpenMS2
PYTHONPATH=OpenMS-build/pyOpenMS2-build python3 -m pytest src/pyOpenMS/tests/unittests/ -v

# Specific test file
PYTHONPATH=OpenMS-build/pyOpenMS2-build python3 -m pytest src/pyOpenMS2/tests/test_msspectrum.py -v
```
