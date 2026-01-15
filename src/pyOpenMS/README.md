# pyopenms

pyOpenMS is a Python library for the analysis of mass spectrometry data. It is mainly based on Cython wrappers around the OpenMS C++ library. To see which classes and functions are currently wrapped, please check the `.pxd` files under `./pxds` or consult our [API documentation](https://pyopenms.readthedocs.io/en/latest/apidocs/index.html).

Additionally, it provides convenience functions for plotting and converting from/to DataFrames or NumPy/pyarrow arrays.

## Build Instructions

There are two main ways to build pyOpenMS:

### 1. Standalone Build with `uv` (For Users and Distribution)

The recommended way to build pyOpenMS as a standalone Python package is using `uv`:

```bash
# Install uv if you don't have it
pip install uv

# Build wheel from the OpenMS/src/pyOpenMS directory
cd src/pyOpenMS
uv build

# Install the wheel
pip install dist/pyopenms-*.whl
```

This method uses the PEP 517 build system defined in `pyproject.toml` and will:

- Automatically discover or download an OpenMS installation
- Handle all build dependencies
- Create a distributable wheel package

### 2. In-Tree Development Build with CMake (For OpenMS Developers)

When developing OpenMS and pyOpenMS together, use the integrated CMake build.

**Prerequisites:**

- A **built OpenMS library** (the C++ library must be compiled or installed first)
- OpenMS build directory or installation

**Required CMake Variables:**

If OpenMS is not installed system-wide, you may need to specify:

- `-DOpenMS_DIR="/path/to/OpenMS-build"` - Path to OpenMS build/install directory with CMake config files
- `-DCMAKE_PREFIX_PATH="/path/to/qt6;/path/to/other/deps"` - Semicolon-separated paths to dependencies (Qt6, etc.)

Depending on how OpenMS was built, you may also need:

- `-DOpenMP_ROOT="/path/to/openmp"` - Path to OpenMP installation (e.g., libomp on macOS)

**Build Commands:**

```bash
# From your OpenMS build directory (e.g., OpenMS-build/)
cmake -DPYOPENMS=ON /path/to/OpenMS/source

# Build pyOpenMS
cmake --build . --target pyopenms
```

**Optional CMake Configuration:**

- `-DPython_EXECUTABLE="/path/to/python"` - Specify Python interpreter
- `-DPY_NUM_THREADS=2` - Parallel Cython compilation (requires 16GB+ RAM)
- `-DPY_NUM_MODULES=8` - Number of module splits (default: 8)

**Available CMake Targets when `PYOPENMS=ON`:**

| Target | Description |
|--------|-------------|
| `pyopenms` | Main target: builds complete pyOpenMS module (runs all sub-targets) |
| `pyopenms_compile` | Compile all C++ extension module(s) - aggregates `pyopenms_1` through `pyopenms_N` |
| `pyopenms_1` ... `pyopenms_N` | Individual extension module targets (when `PY_NUM_MODULES > 1`) |
| `compile_pxds` | Generate `.pyx` file(s) from `.pxd` declarations using autowrap |
| `pyopenms_copy_deps` | Copy OpenMS libraries to pyOpenMS build directory (when `NO_DEPENDENCIES=OFF`) |
| `pyopenms_fix_deps` | Fix library dependencies on macOS (when `NO_DEPENDENCIES=OFF`) |
| `pyopenms_wheel` | Package pyOpenMS as a wheel file (depends on `pyopenms`, creates wheel in `pyopenms_wheels/`) |

**Available CMake Targets when `PYOPENMS=OFF` (default):**

Only the OpenMS C++ library targets are available. No pyOpenMS-related targets exist.

**Setting up for development:**

1. **(Optional) Create a virtual environment:**

   ```bash
   python -m venv /path/to/myenv
   source /path/to/myenv/bin/activate  # Linux/macOS
   # or: c:\path\to\myenv\Scripts\activate.bat  # Windows
   ```

2. **Install build dependencies:**

   ```bash
   cd src/pyOpenMS
   # Option 1: Using uv (recommended)
   pip install uv
   uv sync --only-group build
   
   # Option 2: Using pip with development dependencies
   pip install -e .[dev]
   ```

3. **Configure and build:**

   ```bash
   cd ../../OpenMS-build  # Or your build directory
   cmake -DPYOPENMS=ON ../OpenMS
   cmake --build . --target pyopenms
   ```

4. **Run tests:**

   ```bash
   ctest -R pyopenms
   
   # Or run pytest directly for faster iteration:
   cd pyOpenMS
   python -m pytest tests/unittests
   python -m pytest tests/integration_tests
   ```

5. **Install for development:**

   ```bash
   pip install -e pyopenms --no-cache-dir --no-binary=pyopenms
   ```

## How pyOpenMS is Built Under the Hood

The build process involves several steps that transform C++ code into a Python extension module:

### Step 1: Wrapper Generation (`compile_pxds`)

- **Input:** `.pxd` declaration files in `src/pyOpenMS/pxds/`
- **Tool:** `autowrap` (automated Cython wrapper generator)
- **Output:** `pyopenms/pyopenms.pyx` (Cython source file) or multiple split `.pyx` files
- **What happens:** autowrap reads the `.pxd` files containing class and method declarations and generates Cython wrapper code with proper type conversions between C++ and Python

### Step 2: Addon Injection

- **Input:** Manual additions in `src/pyOpenMS/addons/` (`.pyx` files)
- **Output:** Addons are merged into `pyopenms.pyx` or split `.pyx` files
- **What happens:** Python-specific convenience methods (like `get_df()`, `__repr__()`) are injected into the generated wrapper code

### Step 3: Cython Compilation (part of `pyopenms_compile`)

- **Input:** `pyopenms.pyx` (or `_pyopenms_1.pyx` through `_pyopenms_N.pyx` for split modules)
- **Tool:** Cython compiler
- **Output:** C++ source files (`.cpp`)
- **What happens:** Cython translates the `.pyx` file(s) into C++ code that bridges Python and C++ types. Type stub files (`.pyi`) are also generated for IDE support.

### Step 4: C++ Compilation and Linking (`pyopenms_1` ... `pyopenms_N`)

- **Input:** Generated `.cpp` files
- **Tools:** C++ compiler (gcc/clang/MSVC), linker
- **Output:** `_pyopenms.*.so` (Linux), `_pyopenms.*.dylib` (macOS), or `_pyopenms.*.pyd` (Windows)
  - For split modules: `_pyopenms_1`, `_pyopenms_2`, ... `_pyopenms_N`
- **What happens:** The C++ code is compiled and linked against:
  - OpenMS C++ library
  - OpenSwathAlgo library
  - Python C API
  - NumPy C API
  - Cython runtime

### Step 5: Dependency Bundling (optional, when `NO_DEPENDENCIES=OFF`)

- **Targets:** `pyopenms_copy_deps`, `pyopenms_fix_deps` (macOS only)
- **What happens:** 
  - OpenMS shared libraries are copied to the pyOpenMS module directory
  - On macOS, `install_name_tool` adjusts library paths for relocatability
  - This creates a standalone Python package with all dependencies included

The result is a native Python extension module that can be imported with `import pyopenms`.

## Wrapping New Classes

To add new OpenMS classes to pyOpenMS, you need to create `.pxd` declaration files that tell autowrap which classes and methods to wrap.

**Quick overview:**

1. Create or edit a `.pxd` file in `src/pyOpenMS/pxds/`
2. Declare the class and methods using Cython syntax
3. Add wrapping hints as comments (e.g., `# wrap-ignore`, `# wrap-doc:`)
4. Rebuild: `cmake --build . --target pyopenms`

**Detailed instructions:** See [README_WRAPPING_NEW_CLASSES](./README_WRAPPING_NEW_CLASSES)

**autowrap documentation:** https://github.com/OpenMS/autowrap/blob/master/docs/README.md

**Important wrapping hints:**

| Hint | Purpose | Example |
|------|---------|---------|
| `# wrap-ignore` | Skip this method | `void internal() # wrap-ignore` |
| `# wrap-doc:` | Add Python docstring | See autowrap docs for format |
| `# wrap-as:NewName` | Rename method | `void getValue() # wrap-as:get_value` |
| `# wrap-iter-begin/end` | Enable iteration | For container classes |
| `# wrap-instances:` | Template instantiation | `# wrap-instances:T:int,double` |

**Common patterns:**

- Always declare default and copy constructors
- Use `cimport` for Cython imports, not Python `import`
- Match parameter types exactly with C++ signatures
- Use `except + nogil` for methods that may throw exceptions

**After making changes:**

Force regeneration by removing the generation marker:

```bash
rm OpenMS-build/pyOpenMS/.cpp_extension_generated
cmake --build OpenMS-build --target pyopenms
```

## Development Patterns

### Naming Conventions

Use **lowercase snake_case** for all Python-facing names to follow PEP 8:

| Type | Convention | Examples |
|------|------------|----------|
| DataFrame columns | `snake_case` | `precursor_mz`, `native_id`, `ion_mobility` |
| Method names | `snake_case` | `get_peaks()`, `get_data_dict()`, `get_df()` |
| Variables | `snake_case` | `peak_count`, `meta_values` |

**Note**: C++ OpenMS uses camelCase (e.g., `getPrecursorMZ()`), but Python convenience methods
and DataFrame columns should use snake_case for Pythonic consistency.

### DataFrame Export Pattern (get_data_dict + get_df)

To add DataFrame export to a class, implement both methods directly in the Cython addon file:

1. **`get_df_columns()`**: Returns list of available column names (for discovery)
2. **`get_data_dict(columns=None)`**: Returns dict of numpy arrays (works without pandas)
3. **`get_df(columns=None)`**: Returns pandas DataFrame (imports pandas lazily)

This pattern ensures:

- Users without pandas can still access data via `get_data_dict()`
- Column selection happens at the data extraction level for efficiency
- All DataFrame logic is in one place (the Cython addon)

**Example addon** (`addons/MyClass.pyx`):

```cython
cimport numpy as np
import numpy as np
import pandas as pd

    def get_df_columns(self, columns='default'):
        """Returns list of column names that get_df() would produce."""
        cols = ['mz', 'intensity']
        if columns == 'all':
            cols.append('extra_data')
        return cols

    def get_data_dict(self, columns=None):
        """Returns dict of numpy arrays for DataFrame conversion."""
        if columns is not None:
            requested = set(columns)
        else:
            requested = None

        def want(col):
            return requested is None or col in requested

        data = {}
        if want('mz'):
            data['mz'] = self.get_mz_array()
        if want('intensity'):
            data['intensity'] = self.get_intensity_array()
        return data

    def get_df(self, columns=None):
        """Returns pandas DataFrame."""
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

After modifying addon `.pyx` files, force regeneration:

```bash
rm OpenMS-build/pyOpenMS/.cpp_extension_generated
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
