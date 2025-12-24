# OpenMS Agent Notes

This file provides context and instructions for AI coding agents working on OpenMS. It follows the [AGENTS.md](https://agents.md) standard.

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

## Project Stack

- **Language**: C++20, Python 3.9+
- **Build**: CMake 3.24+, out-of-tree builds in `OpenMS-build/`
- **Testing**: CTest, GoogleTest-style macros, pytest for Python
- **Style**: `.clang-format` in repo root, cpplint via `ENABLE_STYLE_TESTING=ON`
- **Platforms**: Linux, macOS (Apple Clang), Windows (MSVC 2019+)

## Repository Layout

```
OpenMS/
├── src/
│   ├── openms/           # Core C++ library
│   │   ├── include/OpenMS/  # Headers (.h)
│   │   └── source/          # Implementation (.cpp)
│   ├── openms_gui/       # Qt-based GUI components
│   ├── openswathalgo/    # OpenSWATH algorithms
│   ├── topp/             # Command-line tools (TOPP)
│   ├── pyOpenMS/         # Python bindings
│   │   ├── pxds/         # .pxd declarations for autowrap
│   │   ├── addons/       # Python-only method additions
│   │   └── tests/        # Python tests
│   └── tests/
│       ├── class_tests/openms/source/  # C++ unit tests
│       └── topp/         # TOPP integration tests
├── cmake/                # CMake modules
├── doc/                  # Documentation source
└── share/OpenMS/         # Runtime data files
```

## Code Style (with Examples)

**Naming conventions:**
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

**File structure:**
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

**Formatting rules:**
- 2 spaces indentation, no tabs
- Unix line endings (LF)
- Braces on their own lines, aligned
- Space after keywords (`if`, `for`, `while`)
- Always use braces, even for single-line blocks

## Testing Patterns

**Unit test structure:**
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

**Adding tests:**
1. Create `src/tests/class_tests/openms/source/ClassName_test.cpp`
2. Add to `src/tests/class_tests/openms/executables.cmake`
3. Use `NEW_TMP_FILE(filename)` for temp output files
4. Test data goes in `src/tests/class_tests/libs/data/` (prefix with class name)

## Git Workflow

**Commit message format:**
```
[TAG1,TAG2] Short summary (<=80 chars preferred)

Longer description explaining why, not what.

Fixes #123
```

**Valid tags:** `NOP`, `DOC`, `COMMENT`, `API`, `INTERNAL`, `FEATURE`, `FIX`, `TEST`, `FORMAT`, `PARAM`, `IO`, `LOG`, `GUI`, `RESOURCE`, `BUILD`

**Branch workflow:**
- Fork the repo, branch from `develop`
- Open PRs against `develop` (Gitflow)
- Minimize pushes on open PRs (CI is resource-heavy)

## Change Impact Checklist

When you change | Also update
----------------|------------
C++ class (new) | Add `.h`/`.cpp`, Doxygen docs, class test, `OPENMS_DLLAPI`, CMake registration
C++ API | `.pxd` files, pyOpenMS addons, tests, docs
TOPP tool (new) | `src/topp/executables.cmake`, `ToolHandler.cpp`, docs, TOPP tests
Parameters | Tool docs, CTD, tests, `CHANGELOG`
File format | `FileHandler::NamesOfTypes[]`, schemas, tests

## pyOpenMS Wrapping

**Key files:**
- `.pxd` declarations: `src/pyOpenMS/pxds/`
- Python addons: `src/pyOpenMS/addons/`
- Type converters: `src/pyOpenMS/converters/`

**Common patterns:**
```python
# In addons/MyClass.pyx - inject Python-only methods
def get_df(self):
    """Return pandas DataFrame."""
    import pandas as pd
    return pd.DataFrame(self.get_data_dict())
```

**Gotchas:**
- Always declare default and copy constructors in `.pxd`
- Use `cimport`, not Python `import` for Cython imports
- Autowrap returns Python strings; do NOT call `.decode()`
- Use snake_case for Python-facing names

## Verification Commands

After making changes, verify with:
```bash
# Check formatting
clang-format --dry-run -Werror <changed-files>

# Run relevant tests
ctest -R <ClassName> -V

# For pyOpenMS changes
cd OpenMS-build && ctest -R pyopenms -V

# Style check
cmake --build OpenMS-build --target test_style
```

## Key Documentation

**In-repo docs:**
- `README.md` - Project overview
- `CONTRIBUTING.md` - Contribution guidelines
- `src/pyOpenMS/README.md` - pyOpenMS development
- `src/pyOpenMS/README_WRAPPING_NEW_CLASSES` - Wrapping guide

**Online resources:**
- [OpenMS Documentation](https://openms.readthedocs.io/en/latest)
- [pyOpenMS API Reference](https://pyopenms.readthedocs.io/en/latest/apidocs/index.html)
- [Developer Coding Conventions](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/developer_coding_conventions.html)
- [How to Write Tests](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/developer_how_to_write_tests.html)
- [GitHub Wiki](https://github.com/OpenMS/OpenMS/wiki)

## Common Gotchas

1. **Template methods with 2+ args in tests**: Wrap in parentheses for `START_SECTION`
2. **GUI tests need display**: Set `QT_QPA_PLATFORM=minimal` for headless runs
3. **pyOpenMS tests shadow imports**: Run from outside source tree with `PYTHONPATH` set
4. **Windows paths**: Keep build paths short; use 64-bit only
5. **FuzzyDiff for numeric tests**: Build `all`/`ALL_BUILD` to include it

## Debugging Tips

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
