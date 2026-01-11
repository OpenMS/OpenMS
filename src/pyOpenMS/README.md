General
-------

pyopenms is a Python library for the analysis of mass spectrometry data.
It is mainly based on Cython wrappers around the OpenMS C++ library. To see which classes and functions are
currently wrapped, please check the pxd files under "./pxds" or consult our
[API documentation](https://pyopenms.readthedocs.io/en/latest/apidocs/index.html).
Additionally, it provides some convenience functions for plotting or converting from/to dataframes or numpy arrays.

Wrapping new classes
--------------------

See [README_WRAPPING_NEW_CLASSES](./README_WRAPPING_NEW_CLASSES)

Build instructions
------------------

0. (optional) Create a virtual python environment:
    
    ```bash
    python -m venv /path/to/myenv
    
    # ... and activate it, e.g.
    # Linux:
    source <venv>/bin/activate
    # Windows:
    c:\path\to\myenv\Scripts\activate.bat
    ```
    
1. Get Python 3.7+ and the following Python libraries:
   
   pip install -r requirements_bld.txt

2. If running from an OpenMS build tree (recommended), just reconfigure with

   ```bash
   cmake -DPYOPENMS=ON .
   ```
   
   If it does not find the python that you just installed or complains about libraries not located although
   you just installed them, help it find the correct python executable by adding `-DPython_EXECUTABLE="/path/to/python(.exe)"`
   If your computer has a lot of RAM (16GB+) you can add `-DPY_NUM_THREADS=2`
   (or up to the number of split modules, which are by default PY_NUM_MODULES=8).

   Building with an existing, installed OpenMS library is possible but not well-tested. All you need is the current pyOpenMS
   directory (where this README is located). Configure CMake in a build dir of your choice with general CMake
   and pyopenms-related options only. It will try to find OpenMS based on its CMake config files that should be installed
   with newer versions of OpenMS (around 2.8+) and hopefully parse all necessary options and find OpenMS' transitive dependencies.
   The rest is the same.

4. Build CMake target "pyopenms" build-system agnostic with

   ```bash
   cmake --build . --target pyopenms
   ```

5. Run tests with

   ```bash
   ctest -R pyopenms
   ```

   "-R" to restrict to pyopenms* tests. If running out of the OpenMS build tree, this should not be necessary.
   
   **Alternative: Running tests directly with pytest**
   
   For development and debugging, you can run tests directly with pytest:
   
   ```bash
   cd <build-dir>/pyOpenMS
   python -m pytest tests/unittests
   python -m pytest tests/integration_tests
   ```
   
   Tests automatically locate test data files using the `conftest.py` configuration. If test data
   is in a non-standard location, set the `OPENMS_TEST_DATA_PATH` environment variable:
   
   ```bash
   export OPENMS_TEST_DATA_PATH=/path/to/OpenMS/src/tests/topp
   python -m pytest tests/
   ```

6. Install locally (and in-place for live edits [option -e]) into current Python with

   ```bash
   pip install -e pyopenms --no-cache-dir --no-binary=pyopenms
   ```
   
   `--no-binary` is used because the binaries are/were built with CMake.

Development Patterns
--------------------

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

