# pyopenms

pyOpenMS is a Python library for the analysis of mass spectrometry data. It uses nanobind for C++ bindings to the OpenMS C++ library. To see which classes and functions are currently wrapped, check the `bindings/bind_*.cpp` files or consult our [API documentation](https://pyopenms.readthedocs.io/en/latest/apidocs/index.html).

Additionally, it provides convenience functions for plotting and converting from/to DataFrames or NumPy/pyarrow arrays.

## Build Instructions

There are two main ways to build pyOpenMS:

### 1. Standalone Build with `uv` (For Users and Distribution)

**Prerequisites:**

- A **built OpenMS library** (the C++ library must be compiled or installed first)
- OpenMS build directory or installation

**Required CMake Variables:**

If OpenMS is not installed system-wide, you may need to specify:

- `-DOpenMS_DIR="/path/to/OpenMS-build"` - Path to OpenMS build/install directory with CMake config files
- `-DCMAKE_PREFIX_PATH="/path/to/qt6;/path/to/other/deps"` - Semicolon-separated paths to dependencies (Qt6, etc.)

Depending on how OpenMS was built, you may also need:

- `-DOpenMP_ROOT="/path/to/openmp"` - Path to OpenMP installation (e.g., libomp on macOS)

#### Usage

The recommended way to build pyOpenMS as a standalone Python package is using `uv`:

```bash
# Install uv if you don't have it
pip install uv

# Build wheel from the OpenMS/src/pyOpenMS directory
cd src/pyOpenMS
uv build  # add py-build-cmake options here to discover OpenMS+deps and tweak compilation options

# Install the wheel
pip install dist/pyopenms-*.whl
```

This method uses the PEP 517 build system defined in `pyproject.toml` and will:

- Automatically discover or download an OpenMS installation
- Handle all build dependencies
- Create a distributable wheel package

### 2. In-Tree Development Build with CMake (For OpenMS Developers)

When developing OpenMS and pyOpenMS together, use the integrated CMake build from the same build
directory as OpenMS.

**Prerequisites:**

We recommend using uv which will be looked up at CMake configure time and creates a nice venv for you
automatically.

```bash
# Clone OpenMS to e.g. OpenMS-source
# gh repo clone OpenMS/OpenMS OpenMS-source

# Familiarize with OpenMS C++ build settings

# Install uv
pipx install uv
```

**Build Commands:**

```bash
# Optional: Create a .python-version file in /path/to/OpenMS-source/src/pyOpenMS to steer the python version used by UV
echo "3.12" > /path/to/OpenMS-source/src/pyOpenMS/.python-version

# CMake Configure from your OpenMS build directory (e.g., OpenMS-build/)
# If you have not yet set up an OpenMS build directory, you need to create one first and probably also want to configure
# some OpenMS C++ related build options.
cd /path/to/OpenMS-build
cmake -DPYOPENMS=ON /path/to/OpenMS-source

# UV will create a .venv under /path/to/OpenMS-build/pyOpenMS/.venv with an appropriate python and
# all python build dependencies and sets the CMake Cache variable Python_EXECUTABLE to the corresponding python exe

# Build pyOpenMS
cmake --build . --target pyopenms # or pyopenms_wheel if you want to build wheels, too
```

**Optional CMake Configuration:**

In addition to the usual CMake options you can set for the OpenMS C++ toolkit, enabling pyopenms offers the following:

- `-DPython_EXECUTABLE="/path/to/python"` - Specify Python interpreter
- `-DNO_DEPENDENCIES=ON`- When not distributing a wheel, you can use this to avoid copying dependencies into the pyopenms build folder. Make sure that the pyopenms shared modules find their dependencies at the original places with correct RPATH/INSTALL_NAME_DIR CMake settings.
- `-DWITH_UV=OFF` - Do not use uv to create a new venv. If disabled, make sure the found (or specified, see Python_EXECUTABLE) Python executable has access to all required dependencies.
- `-DPYOPENMS_UV_PYTHON_VERSION=3.12` - Specify the python version that uv should use to create the venv. This will decide with which python version the extension module and the pyopenms wheel will be compatible with. Note: If such a python version is not available on the system, uv will download it for you.

**Available CMake Targets when `PYOPENMS=ON`:**

| Target | Description |
|--------|-------------|
| `pyopenms` | Main target: builds complete pyOpenMS module (runs all sub-targets) |
| `pyopenms_compile` | Compile all C++ nanobind extension modules (10 domain modules + main + arrow) |
| `_pyopenms_kernel`, `_pyopenms_chemistry`, etc. | Individual domain module targets |
| `pyopenms_copy_deps` | Copy OpenMS libraries to pyOpenMS build directory (when `NO_DEPENDENCIES=OFF`) |
| `pyopenms_fix_deps` | Fix library dependencies on macOS (when `NO_DEPENDENCIES=OFF`) |
| `pyopenms_wheel` | Package pyOpenMS as a wheel file (depends on `pyopenms`, creates wheel in `$BUILDDIR/pyopenms_wheels/`) |

**Run tests:**

   ```bash
   # With ctest (pyopenms specific tests only)
   cd /path/to/OpenMS-build
   ctest -R pyopenms # add -V for verbose output

   # Or run pytest directly for faster iteration:
   cd /path/to/OpenMS-build/pyOpenMS
   python -m pytest tests/unittests
   python -m pytest tests/integration_tests
   # Note that if you want to make changes to non-compile .py files, you should
   # make these changes in the OpenMS-source folder instead and reconfigure the
   # with CMake from your build directory such that it syncs.
   ```

**Install for development:**

   ```bash
   # in OpenMS-build/pyOpenMS
   pip install -e pyopenms --no-cache-dir --no-binary=pyopenms
   ```

## How pyOpenMS is Built Under the Hood

The build process compiles hand-maintained nanobind C++ binding files:

### Step 1: C++ Compilation and Linking

- **Input:** Hand-maintained `bindings/bind_*.cpp` files (11 files across 10 domains + main)
- **Tools:** C++ compiler (gcc/clang/MSVC), linker, nanobind
- **Output:** Domain-based shared modules: `_pyopenms.*.so` (Linux), `.dylib` (macOS), `.pyd` (Windows)
  - 10 domain modules: `_pyopenms_kernel`, `_pyopenms_chemistry`, `_pyopenms_analysis`, etc.
  - 1 main module: `_pyopenms`
  - 1 optional Arrow module: `_arrow_zerocopy` (when `WITH_PARQUET=ON`)
- **What happens:** nanobind C++ code is compiled and linked against OpenMS, OpenSwathAlgo, and Python. All modules share types via `NB_DOMAIN "pyopenms"`.

### Step 2: Addon Injection (at import time)

- **Input:** Pure Python addon files in `pyopenms/addons/`
- **What happens:** When `import pyopenms` runs, the addon system injects Python convenience methods (like `to_df()`, `__repr__()`) into the C++ wrapper classes using the `@addon("ClassName")` decorator.

### Step 3: Dependency Bundling (optional, when `NO_DEPENDENCIES=OFF`)

- **Targets:** `pyopenms_copy_deps`, `pyopenms_fix_deps` (macOS only)
- **What happens:**
  - OpenMS shared libraries are copied to the pyOpenMS module directory
  - On macOS, `install_name_tool` adjusts library paths for relocatability
  - This creates a standalone Python package with all dependencies included

The result is a native Python extension module that can be imported with `import pyopenms`.

## Wrapping New Classes

Bindings are hand-maintained in `bindings/bind_<domain>.cpp` files. See [README_WRAPPING_NEW_CLASSES](./README_WRAPPING_NEW_CLASSES) for detailed instructions.

**Quick overview:**

1. Pick the right `bindings/bind_<domain>.cpp` based on the C++ header path
2. Add the `#include` for the C++ header
3. Add `nb::class_<...>(m, "ClassName", "docstring")` with `.def()` chains
4. Rebuild: `cmake --build OpenMS-build --target pyopenms`

## Development Patterns

### Naming Conventions

Use **lowercase snake_case** for all Python-facing names to follow PEP 8:

| Type | Convention | Examples |
|------|------------|----------|
| DataFrame columns | `snake_case` | `precursor_mz`, `native_id`, `ion_mobility` |
| Method names | `snake_case` | `get_peaks()`, `get_data_dict()`, `to_df()` |
| Variables | `snake_case` | `peak_count`, `meta_values` |

**Note**: C++ OpenMS uses camelCase (e.g., `getPrecursorMZ()`), but Python convenience methods
and DataFrame columns should use snake_case for Pythonic consistency.

### DataFrame Export Pattern (get_data_dict + get_df)

To add DataFrame export to a class, implement methods in a pure Python addon file:

1. **`df_columns()`**: Returns list of available column names (for discovery)
2. **`get_data_dict(columns=None)`**: Returns dict of numpy arrays (works without pandas)
3. **`get_df(columns=None)`**: Returns pandas DataFrame (imports pandas lazily)

This pattern ensures:

- Users without pandas can still access data via `get_data_dict()`
- Column selection happens at the data extraction level for efficiency
- All DataFrame logic is in one place (the addon)

**Example addon** (`pyopenms/addons/myclass.py`):

```python
from pyopenms.addons import addon

@addon("MyClass")
def df_columns(self, columns='default'):
    """Returns list of column names that to_df() would produce."""
    cols = ['mz', 'intensity']
    if columns == 'all':
        cols.append('extra_data')
    return cols

@addon("MyClass")
def get_data_dict(self, columns=None):
    """Returns dict of numpy arrays for DataFrame conversion."""
    import numpy as np
    requested = set(columns) if columns is not None else None

    def want(col):
        return requested is None or col in requested

    data = {}
    if want('mz'):
        data['mz'] = self.get_mz_array()
    if want('intensity'):
        data['intensity'] = self.get_intensity_array()
    return data

@addon("MyClass")
def get_df(self, columns=None):
    """Returns pandas DataFrame."""
    import pandas as pd
    return pd.DataFrame(self.get_data_dict(columns=columns))
```

**Standalone utility functions** are kept in `pyopenms/_dataframes.py`:

- `peptide_identifications_to_df()`: Convert PeptideIdentification list to DataFrame
- `update_scores_from_df()`: Update scores in PeptideIdentifications from DataFrame

### Adding Pythonic Methods

Common methods to add for container-like classes:

- `__len__()`: Return `self.size()` or similar
- `__repr__()`: Return `f"ClassName(key_prop={value}, ...)"` with important properties
- `__str__()`: Delegate to `__repr__()` or return simpler output
- `get_data()`: Return safe copy of data (for DataArray classes)
- `get_data_mv()`: Return memory view (fast but unsafe, document lifetime)

### Rebuilding After Addon Changes

Addons are pure Python files in `pyopenms/addons/` — no recompilation needed. Just rebuild (which copies updated files):

```bash
cmake --build OpenMS-build --target pyopenms -j4
```

### Writing Tests

Tests are located in `src/pyOpenMS/tests/`:

- `unittests/` - Unit tests for individual components
- `integration_tests/` - Integration tests requiring test data files
- `memoryleaktests/` - Memory leak detection tests

**Accessing test data files:**

Tests automatically find test data through the `openms_test_data_dir` pytest fixture defined in `conftest.py`:

```python
import pytest
import os

class TestMyFeature(unittest.TestCase):

    @pytest.fixture(autouse=True)
    def setup_test_data(self, openms_test_data_dir):
        """Setup test with test data directory."""
        self.test_file = os.path.join(openms_test_data_dir, "my_test_file.mzML")

    def test_something(self):
        # Use self.test_file here
        pass
```

The fixture automatically locates test data using multiple strategies:

1. `OPENMS_TEST_DATA_PATH` environment variable (for custom locations)
2. Relative path from source tree (works in development)
3. CMake-configured env.py (backwards compatibility)

To override test data location:

```bash
export OPENMS_TEST_DATA_PATH=/path/to/OpenMS/src/tests/topp
python -m pytest tests/
```
