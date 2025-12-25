# OpenMS Agent Notes

This file provides context and instructions for AI coding agents working on OpenMS. It follows the [AGENTS.md](https://agents.md) standard and is intentionally concise; use the linked docs for full details.

## Critical Constraints

**NEVER do these things:**
- Build the project unless explicitly asked (extremely resource-intensive)
- Modify files in `src/openms/extern/` (third-party vendored code)
- Commit secrets, credentials, or `.env` files
- Use `std::cout`/`std::cerr` directly (use OpenMS logging macros)
- Use `std::endl` (use `\n` for performance)
- Add `using namespace` or `using std::...` in header files
- Modify the contrib tree or third-party dependencies
- Skip tests when making code changes

## Quick Commands

```bash
# Configure (from OpenMS-build/ directory, adjust paths as needed)
cmake -DCMAKE_BUILD_TYPE=Debug ../OpenMS

# Build everything (includes tests)
cmake --build . -j$(nproc)

# Run all tests
ctest -j$(nproc)

# Run specific test by name pattern
ctest -R FeatureMap -j4

# Run tests with verbose output
ctest -R MyTest -V

# Run style checks
cmake --build . --target test_style

# Regenerate pyOpenMS after changes
rm pyOpenMS/.cpp_extension_generated
cmake --build . --target pyopenms -j4

# Run pyOpenMS tests
ctest -R pyopenms

# Check code formatting
clang-format --dry-run -Werror src/openms/source/MYFILE.cpp
```

## Key Docs in This Repo

- `README.md`, `CONTRIBUTING.md`, `ARCHITECTURE.MD`, `CODE_OF_CONDUCT.md`, `PULL_REQUEST_TEMPLATE.md`.
- `src/pyOpenMS/README.md`, `src/pyOpenMS/README_WRAPPING_NEW_CLASSES`.
- `share/OpenMS/examples/external_code/README.md`, `src/tests/external/README.md`.
- `dockerfiles/README.md`, `cmake/MacOSX/README.md`, `tools/jenkins/README.MD`.
- Doxygen (if built) in `OpenMS-build/doc/html/` including `index.html`, `developer_coding_conventions.html`, `developer_cpp_guide.html`, `developer_how_to_write_tests.html`, `howto_commit_messages.html`, `developer_faq.html`, `developer_tutorial.html`, `install_linux.html`, `install_mac.html`, `install_win.html`, `pyOpenMS.html`.

## Repo Layout

- Default build directory: `OpenMS-build/` (out-of-tree).
- Core C++: `src/openms/`, `src/openms_gui/`, `src/openswathalgo/`, `src/topp/`.
- Tests: `src/tests/class_tests/openms/`, `src/tests/class_tests/openms_gui/`, `src/tests/topp/`.
- pyOpenMS: `src/pyOpenMS/` with `pxds/`, `addons/`, `pyopenms/`, `tests/`.

```
OpenMS/
├── src/
│   ├── openms/              # Core C++ library
│   │   ├── include/OpenMS/  # Headers (.h)
│   │   └── source/          # Implementation (.cpp)
│   ├── openms_gui/          # Qt-based GUI components
│   ├── openswathalgo/       # OpenSWATH algorithms
│   ├── topp/                # Command-line tools (TOPP)
│   ├── pyOpenMS/            # Python bindings
│   │   ├── pxds/            # .pxd declarations for autowrap
│   │   ├── addons/          # Python-only method additions
│   │   └── tests/           # Python tests
│   └── tests/
│       ├── class_tests/openms/source/  # C++ unit tests
│       └── topp/            # TOPP integration tests
├── cmake/                   # CMake modules
├── doc/                     # Documentation source
└── share/OpenMS/            # Runtime data files
```

## Build and Install

- **CMake minimum**: 3.21; **C++ standard**: C++20
- Out-of-tree build expected in `OpenMS-build/`; build in place for development (install prefixes are for system installs).
- Use `CMAKE_BUILD_TYPE=Debug` for development to keep assertions/pre/post-conditions.
- Dependencies via distro packages or the contrib tree; set `OPENMS_CONTRIB_LIBS` and `CMAKE_PREFIX_PATH` as needed (Qt, contrib).
- pyOpenMS build deps: `src/pyOpenMS/requirements_bld.txt`; enable with `-DPYOPENMS=ON` and optional `-DPY_NUM_THREADS`/`-DPY_NUM_MODULES`.
- Style checks: `ENABLE_STYLE_TESTING=ON` runs cpplint at `src/tests/coding/cpplint.py`.

**Required dependencies:**
- XercesC, Boost (date_time, regex, iostreams), Eigen3 (3.4.0+), libSVM (2.91+), COIN-OR or GLPK, ZLIB, BZip2, Qt6 (6.1.0+)

**Optional:** HDF5 (`-DWITH_HDF5=ON`), Apache Arrow/Parquet (`-DWITH_PARQUET=ON`)

## Platform-Specific Build Gotchas

### Windows
- **MSYS/MinGW NOT supported** — must use Visual Studio environment
- **MSVC 2019+ required** (version 1920+); AddressSanitizer needs this minimum
- **64-bit only**; use Visual Studio generator (not Ninja/Make)
- **Keep build paths short** to avoid path length issues
- **Never mix Release/Debug libraries** — causes stack corruption and segfaults
- Compiler must match between contrib and OpenMS builds
- HDF5 forced to static linking on MSVC
- OpenMP requires `/openmp:experimental` flag (set automatically) for SIMD support
- Nested OpenMP (`MT_ENABLE_NESTED_OPENMP`) defaults to OFF on MSVC

### macOS
- **Apple Clang (Xcode) required**; Homebrew for dependencies
- **AppleClang >= 15.0.0**: Requires `-ld_classic` linker flag (set automatically)
- Remove older Qt versions if they interfere with Qt6
- Qt6 requires `PrintSupport` component for platform plugin
- `QT_QPA_PLATFORM=minimal` helps for headless/remote GUI runs
- Code signing and notarization required for distribution (see `cmake/MacOSX/README.md`)
- `fix_dependencies.rb` script fixes RPATH for relocatable binaries

### Linux
- Package manager preferred for dependencies; contrib is fallback
- `-fPIC` flag applied automatically for shared library compatibility
- `QT_QPA_PLATFORM=minimal` for headless GUI test runs
- STL debug mode (`_GLIBCXX_DEBUG`) only supported with GCC in Debug builds
- System libraries (libc, libstdc++, libpthread, etc.) excluded from packaging

### Qt6 Issues
- **Minimum version**: 6.1.0
- If Qt6 not found: `-DCMAKE_PREFIX_PATH='<path_to_Qt6_lib_parent>'`
- WebEngineWidgets optional; if missing, JavaScript views disabled in TOPPView (warning only)
- Required components: Core, Network; GUI components need Widgets, Svg, OpenGLWidgets

### Boost from Homebrew Warning
- Statically linked Boost from system installs (brew) NOT fully supported
- Issue: Boost CMake doesn't expose transitive dependencies as targets
- Workaround: Use `-DBOOST_USE_STATIC=OFF` for shared libraries OR build Boost with contrib

### Common CMake Issues
- **CMAKE_SIZEOF_VOID_P bug**: Variable vanishes on CMake version updates → delete `CMakeFiles/` and `CMakeCache.txt`, rerun cmake
- **Eigen3 version detection**: Build system handles CMake's version checking quirks with Eigen3 4.0+ automatically

## Testing

- Unit/class tests: `src/tests/class_tests/<lib>/source/`, add to `executables.cmake`; data in `src/tests/class_tests/libs/data/` (prefix files with class name).
- TOPP tests: add to `src/tests/topp/CMakeLists.txt`, data in `src/tests/topp/`.
- GUI tests: `src/tests/class_tests/openms_gui/source/` (Qt TestLib).
- Build `all`/`ALL_BUILD` to include tests and `FuzzyDiff` (TOPP tests depend on it).
- Use `NEW_TMP_FILE` for each output file in tests; avoid side effects in comparison macros.
- Run with `ctest`, use `-R` for subset, `-V/-VV` for verbosity, `-C` for multi-config generators.
- Use `FuzzyDiff` for numeric comparisons; keep test data small; use whitelist for unstable lines.
- Test templates: `tools/create_test.php` (requires `make xml`).
- `START_SECTION` macro pitfalls: wrap template methods with 2+ arguments in parentheses.
- pyOpenMS tests: `ctest -R pyopenms` or `pytest` with `PYTHONPATH=/path/to/OpenMS-build/pyOpenMS` (run outside the source tree to avoid shadowing).

**Unit test example:**
```cpp
// src/tests/class_tests/openms/source/MyClass_test.cpp
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/PATH/TO/MyClass.h>

START_TEST(MyClass, "$Id$")

MyClass* ptr = nullptr;

START_SECTION(MyClass())
  ptr = new MyClass();
  TEST_NOT_EQUAL(ptr, nullptr)
END_SECTION

START_SECTION(void process(const MSSpectrum&))
  MSSpectrum spec;
  spec.push_back(Peak1D(100.0, 1000.0));
  ptr->process(spec);
  TEST_EQUAL(spec.size(), 1)
END_SECTION

delete ptr;
END_TEST
```

## Coding Conventions

- Indentation: 2 spaces, no tabs; Unix line endings.
- Spacing: after keywords (`if`, `for`) and around binary operators.
- Braces: opening/closing braces align; use braces even for single-line blocks (trivial one-liners may stay single-line).
- File names: class name matches file name; one class per file; always pair `.h` with `.cpp`.
- Templates: use `_impl.h` only when needed; `.h` must not include `_impl.h`.
- Names: classes/types/namespaces in PascalCase; methods lowerCamel; variables snake_case; private/protected members end with `_`.
- Enums and macros uppercase with underscores; avoid the preprocessor; prefer `enum class`.
- Parameters: lower_case with underscores; document ranges/units.
- File extensions: lowercase, except `ML`/`XML` and `mzData`.
- Use OpenMS primitive types from `OpenMS/CONCEPT/Types.h`.
- No `using namespace` or `using std::...` in headers; allowed in `.cpp`.
- Follow Rule-of-0 or Rule-of-6.
- Accessors: get/set pairs for protected/private members; no reference getters for primitive types.
- Exceptions: derive from `Exception::Base`; throw with file/line/`OPENMS_PRETTY_FUNCTION`; catch by reference; document possible exceptions.
- Doxygen: `@brief` + blank line + details; use `@defgroup/@ingroup`; use `.doxygen` files for free-standing docs; `@todo` includes assignee name.
- Comments: at least ~5% of code, use `//` style, plain English describing the next few lines.
- Each file preamble contains the `$Maintainer:$` marker.
- Formatting: use `./.clang-format` in supporting IDEs.

## Doxygen Documentation Style

OpenMS uses `/** */` block comments with `@` tags (not `\` backslash). `@brief` is **required** (not auto-generated from first line).

**File header (required in every .h file):**
```cpp
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Your Name $
// $Authors: Original Author, Your Name $
// --------------------------------------------------------------------------
```

**Class documentation:**
```cpp
/**
  @brief An algorithm to decharge features (i.e. as found by FeatureFinder).

  Detailed description goes here after a blank line.
  Can span multiple lines.

  @htmlinclude OpenMS_FeatureDeconvolution.parameters

  @ingroup Analysis
*/
class OPENMS_DLLAPI FeatureDeconvolution : public DefaultParamHandler
```

**Method documentation with parameters:**
```cpp
/**
  @brief Compute a zero-charge feature map from charged features.

  Find putative ChargePairs, then score them and hand over to ILP.

  @param[in] fm_in      Input feature-map
  @param[out] fm_out    Output feature-map (sorted by position)
  @param[in,out] cons   Consensus map modified in place

  @return The number of charge groups found

  @throws Exception::MissingInformation if RT/MZ data missing
  @throws Exception::InvalidParameter if threshold < 0

  @note The original sequence is saved as MetaValue.
  @warning This method modifies fm_out in place.
*/
Size compute(const FeatureMap& fm_in, FeatureMap& fm_out, ConsensusMap& cons);
```

**Parameter direction tags:** Always use `[in]`, `[out]`, or `[in,out]` for all parameters.

**Grouping constructors/destructors:**
```cpp
/** @name Constructors and Destructors
*/
//@{
/// Default constructor
FeatureDeconvolution();

/// Copy constructor
FeatureDeconvolution(const FeatureDeconvolution& source);

/// Destructor
~FeatureDeconvolution() override;
//@}
```

**Simple inline documentation:** Use `///` for brief single-line docs:
```cpp
/// Fragment mass tolerance for spectrum comparisons
double fragment_mass_tolerance_;

/// Is fragment mass tolerance given in ppm (or Da)?
bool fragment_tolerance_ppm_;
```

**Common Doxygen tags:**
| Tag | Usage |
|-----|-------|
| `@brief` | Required first line summary |
| `@param[in/out]` | Parameter with direction |
| `@return` | Return value description |
| `@throws` / `@exception` | Exceptions that may be thrown |
| `@note` | Important notes |
| `@warning` | Warnings about usage |
| `@ingroup` | Category grouping (e.g., `Analysis_ID`) |
| `@see` | Cross-references |
| `@todo` | Include assignee name: `@todo JohnDoe fix this` |

**Naming examples:**
```cpp
// Classes/Types/Namespaces: PascalCase
class FeatureMap;
namespace OpenMS { }

// Methods: lowerCamelCase
void processSpectrum();

// Variables: snake_case
int peak_count = 0;

// Private/protected members: trailing underscore
double intensity_;

// Enums/macros: UPPER_SNAKE_CASE
enum class Status { RUNNING, COMPLETE };
#define OPENMS_DLLAPI
```

**File structure example:**
```cpp
// MyClass.h - Header file
#pragma once
#include <OpenMS/KERNEL/MSSpectrum.h>

namespace OpenMS
{
  class OPENMS_DLLAPI MyClass  // Export macro required
  {
  public:
    MyClass();
    void process(const MSSpectrum& spectrum);

  private:
    double threshold_;  // Trailing underscore
  };
}

// MyClass.cpp - Implementation file
#include <OpenMS/PATH/TO/MyClass.h>
using namespace OpenMS;  // OK in .cpp files

MyClass::MyClass() : threshold_(0.0) {}

void MyClass::process(const MSSpectrum& spectrum)
{
  // 2-space indentation, braces on own lines
  if (spectrum.empty())
  {
    OPENMS_LOG_WARN << "Empty spectrum\n";  // Use logging macros
    return;
  }
}
```

## C++ Guide (OpenMS-specific)

- `OPENMS_DLLAPI` on all non-template exported classes/structs/functions/vars; not on templates; include in friend operator declarations.
- Use OpenMS logging macros and `OpenMS::LogStream`; avoid `std::cout/err` directly.
- Use `ProgressLogger` in tools for progress reporting.
- Avoid `std::endl` for performance; prefer `\n`.
- Prefer `OpenMS::String` for numeric formatting and parsing (precision and speed).
- Use `Size`/`SignedSize` for STL `.size()` values.
- Avoid pointers; prefer references.
- Prefer forward declarations in headers; include only base class headers, non-pointer members, and templates.

## TOPP Tool Development

- Add new tool source (e.g., `src/topp/<Tool>.cpp`) and register in `src/topp/executables.cmake`.
- Register tool in `src/openms/source/APPLICATIONS/ToolHandler.cpp` to generate Doxygen help output.
- Define parameters in `registerOptionsAndFlags_()`; read with `getStringOption_` and related helpers.
- Document the tool and add to `doc/doxygen/public/TOPP.doxygen` where applicable.
- Add TOPP tests in `src/tests/topp/CMakeLists.txt`.

## pyOpenMS Wrapping

- Autowrap reads `.pxd` in `src/pyOpenMS/pxds/` and generates `pyopenms/pyopenms.pyx` -> `pyopenms.cpp` -> module.
- Addons in `src/pyOpenMS/addons/` inject Python-only methods (indent only; no `cdef class`).
- Keep `.pxd` signatures in sync with C++ APIs; update or remove `wrap-ignore` when wrapping changes.
- Always declare default and copy constructors in `.pxd`; use `cimport`, not Python `import`.
- For non-inheriting classes use `cdef cppclass ClassName:` with no base.
- Autowrap hints: `wrap-ignore`, `wrap-as`, `wrap-iter-begin/end`, `wrap-instances`, `wrap-attach`, `wrap-upper-limit`, `wrap-inherits`.
- Avoid custom `__init__` unless required; it overrides autowrap dispatchers.
- Use snake_case for Python-facing names and DataFrame columns.
- Do not add Python-only methods to `.pxd`; use addons or `_dataframes.py` wrappers.
- DataFrame pattern: `get_data_dict()` in addon returns numpy arrays; `get_df()` in `src/pyOpenMS/pyopenms/_dataframes.py` wraps with pandas.
- Type converters: implement in `src/pyOpenMS/converters/special_autowrap_conversionproviders.py`, register in `src/pyOpenMS/converters/__init__.py`.
- Gotchas: autowrap returns Python strings; do not `.decode()`. Avoid `cdef` for autowrap string returns. Avoid `cdef` typed variables for autowrap return values inside `def` methods; use Python type checks. Keep addons minimal; avoid redundant aliases. `# wrap-doc:` indentation is strict.
- Regenerate after addon changes:
  ```bash
  rm OpenMS-build/pyOpenMS/.cpp_extension_generated
  cmake --build OpenMS-build --target pyopenms -j4
  ```

## Change-Impact Checklist

- New C++ class: add `.h`/`.cpp`, Doxygen docs, class test, `OPENMS_DLLAPI`, register in CMake lists.
- C++ API change: update `.pxd`/addons, pyOpenMS tests, and relevant docs; tag commits with `API` as needed.
- New/changed TOPP tool: register in `src/topp/executables.cmake` and `ToolHandler.cpp`, add docs, add TOPP tests and data.
- Parameter or I/O change: update tool docs/CTD, tests, and `CHANGELOG`; use `PARAM`/`IO` commit tags.
- File format change: update `FileHandler::NamesOfTypes[]`, schemas/validators, and tests.

## Contribution Workflow and Commit Messages

- Development follows Gitflow; use forks and open PRs against `develop`.
- Commit format: `[TAG1,TAG2] short summary` (<=120 chars, <=80 preferred), blank line, longer description, and `Fixes #N`/`Closes #N` when applicable.
- Commit tags: NOP, DOC, COMMENT, API, INTERNAL, FEATURE, FIX, TEST, FORMAT, PARAM, IO, LOG, GUI, RESOURCE, BUILD.
- PR checklist: update `AUTHORS` and `CHANGELOG`, run/extend tests, update pyOpenMS bindings when needed.
- Minimize pushes on open PRs (CI is heavy).
- Run `tools/checker.php` and/or `ENABLE_STYLE_TESTING` for local checks.

**Commit message example:**
```
[TAG1,TAG2] Short summary (<=80 chars preferred)

Longer description explaining why, not what.

Fixes #123
```

## Debugging and Profiling

- Linux: use `ldd` to inspect shared libs; `nm -C` for symbols; `perf`/`hotspot` for profiling.
- Windows: Dependency Walker or `dumpbin /DEPENDENTS` and `dumpbin /EXPORTS`.
- Memory checks: AddressSanitizer or valgrind with `tools/valgrind/openms_external.supp`.

```bash
# Linux: inspect shared libraries
ldd /path/to/binary
nm -C /path/to/library.so | grep MySymbol

# Memory checking
valgrind --suppressions=tools/valgrind/openms_external.supp ./MyTest

# Profile with perf
perf record -g ./MyTool input.mzML
perf report
```

## External Projects and Examples

- Example external CMake project: `share/OpenMS/examples/external_code/`.
- External test project: `src/tests/external/`.
- Use the same compiler/generator as OpenMS; set `OPENMS_CONTRIB_LIBS` and `OpenMS_DIR` when configuring.

## CI, Packaging, and Containers

- CI runs in GitHub Actions; CDash collects nightly results.
- Jenkins packaging uses `tools/jenkins/os_compiler_matrix.tsv` (edit only if needed).
- PR commands/labels: `/reformat`, label `NoJenkins`, comment `rebuild jenkins`.
- Container images: see `dockerfiles/README.md` and GHCR packages.
- macOS code signing/notarization: see `cmake/MacOSX/README.md`.

## Documentation Links (External)

### OpenMS Docs
- http://www.openms.org/
- http://www.OpenMS.de
- https://openms.readthedocs.io/en/latest
- https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/index.html
- https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/index.html
- http://www.openms.de/current_doxygen/html/
- https://pyopenms.readthedocs.io/en/latest/index.html
- https://pyopenms.readthedocs.io/en/latest/apidocs/index.html
- https://abibuilder.cs.uni-tuebingen.de/archive/openms/OpenMSInstaller/
- https://abibuilder.cs.uni-tuebingen.de/archive/openms/OpenMSInstaller/nightly/
- http://www.psidev.info/

### Doxygen Developer Pages (release/latest)
- https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/developer_tutorial.html
- https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/developer_coding_conventions.html
- https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/developer_cpp_guide.html
- https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/developer_how_to_write_tests.html
- https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/howto_commit_messages.html
- https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/developer_faq.html

### Developer Workflow and Contribution
- https://github.com/OpenMS/OpenMS
- https://github.com/OpenMS/OpenMS/issues
- https://github.com/OpenMS/OpenMS/wiki#-for-developers
- https://github.com/OpenMS/OpenMS/wiki/Coding-conventions
- https://github.com/OpenMS/OpenMS/wiki/Write-tests
- https://github.com/OpenMS/OpenMS/wiki/pyOpenMS#wrap
- https://pyopenms.readthedocs.io/en/latest/wrap_classes.html
- https://openms.readthedocs.io/en/latest/contribute-to-openms/pull-request-checklist.html
- https://github.com/OpenMS/OpenMS/wiki/Pull-Request-Checklist
- https://github.com/OpenMS/OpenMS/wiki/Preparation-of-a-new-OpenMS-release#release_developer
- http://nvie.com/posts/a-successful-git-branching-model/
- https://help.github.com/articles/fork-a-repo
- https://help.github.com/articles/syncing-a-fork
- https://help.github.com/articles/using-pull-requests
- http://cdash.seqan.de/index.php?project=OpenMS
- https://github.com/OpenMS/OpenMS/tags

### Build/Install Guides
- https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/install_linux.html
- https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/install_mac.html
- https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/install_win.html
- https://github.com/OpenMS/THIRDPARTY
- https://pkgs.org/search/?q=openms
- http://manpages.ubuntu.com/manpages/hardy/man1/ctest.1.html
- http://www.cmake.org
- http://cmake.org/
- https://visualstudio.microsoft.com/de/downloads/?q=build+tools
- http://www.7-zip.org/
- https://www.qt.io/download
- https://wiki.qt.io/Building_Qt_6_from_Git
- https://developer.apple.com/xcode/
- https://brew.sh/
- http://www.OpenMS.de/download/

### Coding and Tooling
- https://clang.llvm.org/docs/ClangFormat.html
- https://devblogs.microsoft.com/cppblog/clangformat-support-in-visual-studio-2017-15-7-preview-1/
- https://git-scm.com/
- http://www.doxygen.org
- http://www.doxygen.org/index.html
- https://llvm.org/builds/
- https://docs.microsoft.com/en-us/cpp/error-messages/compiler-errors-1/compiler-error-c2471?view=msvc-170
- https://github.com/OpenMS/autowrap/blob/master/docs/README.md

### Testing and Profiling Tools
- https://openms.readthedocs.io/en/latest/docs/topp/adding-new-tool-to-topp.html#how-do-I-add-a-new-TOPP-test
- https://perf.wiki.kernel.org/index.php/Main_Page
- https://github.com/KDAB/hotspot
- http://sandsoftwaresound.net/perf/perf-tutorial-hot-spots/
- http://valgrind.org/docs/manual/
- https://github.com/cbielow/wintime
- http://www.dependencywalker.com/

### Packaging and Containers
- https://github.com/orgs/OpenMS/packages
- https://github.com/OpenMS/NSIS
- http://miktex.org/
- http://www.ghostscript.com/
- http://www.graphviz.org
