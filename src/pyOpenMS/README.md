# pyopenms

pyOpenMS is a Python library for the analysis of mass spectrometry data. It is mainly based on Cython wrappers around the OpenMS C++ library. To see which classes and functions are currently wrapped, please check the `.pxd` files under `./pxds` or consult our [API documentation](https://pyopenms.readthedocs.io/en/latest/apidocs/index.html).

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
- `-DPY_NUM_THREADS=2` - Parallel Cython compilation (only affects the generation of cythonized .cpp files, NOT compilation, steer this with CMake/build program parallelism instead)
- `-DPY_NUM_MODULES=8` - Number of module splits (default: 8)
- `-DNO_DEPENDENCIES=ON`- When not distributing a wheel, you can use this to avoid copying dependencies into the pyopenms build folder. Make sure that the pyopenms shared modules find their dependencies at the original places with correct RPATH/INSTALL_NAME_DIR CMake settings.
- `-DWITH_UV=OFF` - Do not use uv to create a new venv. If disabled, make sure the found (or specified, see Python_EXECUTABLE) Python executable has access to all required dependencies.
- `-DPYOPENMS_UV_PYTHON_VERSION=3.12` - Specify the python version that uv should use to create the venv. This will decide with which python version the extension module and the pyopenms wheel will be compatible with. Note: If such a python version is not available on the system, uv will download it for you.

**Available CMake Targets when `PYOPENMS=ON`:**

| Target | Description |
|--------|-------------|
| `pyopenms` | Main target: builds complete pyOpenMS module (runs all sub-targets) |
| `pyopenms_compile` | Compile all C++ extension module(s) - aggregates `pyopenms_1` through `pyopenms_N` |
| `pyopenms_1` ... `pyopenms_N` | Individual extension module targets (when `PY_NUM_MODULES > 1`) |
| `compile_pxds` | Generate `.pyx` file(s) from `.pxd` declarations using autowrap |
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

The build process involves several steps that transform C++ code into a Python extension module:

### Step 1: Wrapper Generation (`compile_pxds`)

- **Input:** `.pxd` declaration files in `src/pyOpenMS/pxds/`
- **Tool:** `autowrap` (automated Cython wrapper generator)
- **Output:** Internally held corresponding pyx source code for every pxd file. Designed for a 1:1 correspondence of classes and pxds.
- **What happens:** autowrap reads the `.pxd` files containing class and method declarations and generates Cython wrapper code with proper type conversions between C++ and Python

### Step 2: Addon Injection

- **Input:** Manual additions in `src/pyOpenMS/addons/` (`.pyx` files)
- **Output:** Addons are merged into `pyopenms.pyx` or split `.pyx` files
- **What happens:** Python-specific convenience methods (like `to_df()`, `__repr__()`) are injected into the generated wrapper code. Pyx files whose name corresponds to a pxd file are added to the end of the internally held, generated pyx source code for every pxd file/class. Pyx source codes are merged into N "split modules" to reduce peak memory usage during cpp compilation later. Special pyx files are added to the split modules as follows:
  - ADD_TO_FIRST.pyx: Used to create stable references for utility methods. Only gets added to the first split module (i.e. _pyopenms_1)
  - ADD_TO_ALL_OTHER.pyx: Usually used in conjunction with ADD_TO_FIRST to import the utility methods into all but the first split module.
  - ADD_TO_ALL.pyx: Usually additional non-OpenMS imports that are needed by all modules.

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
| Method names | `snake_case` | `get_peaks()`, `get_data_dict()`, `to_df()` |
| Variables | `snake_case` | `peak_count`, `meta_values` |

**Note**: C++ OpenMS uses camelCase (e.g., `getPrecursorMZ()`), but Python convenience methods
and DataFrame columns should use snake_case for Pythonic consistency.

### DataFrame Export Pattern (get_data_dict + get_df)

To add DataFrame export to a class, implement both methods directly in the Cython addon file:

1. **`df_columns()`**: Returns list of available column names (for discovery)
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

    def df_columns(self, columns='default'):
        """Returns list of column names that to_df() would produce."""
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
