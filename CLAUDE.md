- the build folder is OpenMS-build/

## Known Build Issues

### Boost static libs on macOS
Boost's CMake config has incomplete `find_dependency()` calls for transitive dependencies. Use `-DBOOST_USE_STATIC_LIBS=OFF` to avoid linker errors at the final linking stage. This is a 5+ year old upstream Boost issue.

### CMAKE_PREFIX_PATH separators
Per [CMake documentation](https://cmake.org/cmake/help/latest/variable/CMAKE_PREFIX_PATH.html):
- **-D option**: Use semicolons (`;`) as list separators. Example: `-DCMAKE_PREFIX_PATH="/path/one;/path/two"`
- **Environment variable**: Use OS-native path separators (`:` on Unix, `;` on Windows). Example: `export CMAKE_PREFIX_PATH=/path/one:/path/two`

Using the wrong separator causes CMake to treat multiple paths as a single incorrect path.

## pyOpenMS Python Wrapping

### Architecture
- pyOpenMS uses **autowrap** to generate Cython bindings from `.pxd` files
- `.pxd` files are in `src/pyOpenMS/pxds/` and define the C++ class interface
- **Addon files** (`.pyx`) in `src/pyOpenMS/addons/` extend generated classes with custom Python methods

### Adding Custom Python Methods (e.g., `__repr__`, `__len__`)
1. Create/edit an addon file: `src/pyOpenMS/addons/<ClassName>.pyx`
2. Methods are injected into the generated class - indent with 4 spaces, no class declaration needed
3. Example addon file structure:
   ```cython


       def __repr__(self):
           """Return string representation."""
           return f"ClassName(property={self.getProperty()})"
   ```

### DataFrame Export Pattern
- **`_dataframes.py`** (`src/pyOpenMS/pyopenms/_dataframes.py`) wraps Cython classes to add pandas-dependent methods
- Pattern: `get_data_dict()` in Cython addon returns numpy arrays, `get_df()` in `_dataframes.py` wraps it with pandas
- This keeps pandas as an **optional dependency** - users without pandas can still use `get_data_dict()`
- Example: `_MSSpectrumDF` extends `_MSSpectrum` and adds `get_df()` that calls `get_data_dict()`

### Type Converters for Custom C++ Types
When a C++ type (like `DPosition<1>`) needs to accept Python primitives (like `float`):

1. **Create a converter class** in `src/pyOpenMS/converters/special_autowrap_conversionproviders.py`:
   ```python
   class OpenMSDPosition1(TypeConverterBase):
       def get_base_types(self):
           return "DPosition1",  # C++ type name in pxd
       def matches(self, cpp_type):
           return not cpp_type.is_ptr
       def matching_python_type(self, cpp_type):
           return ""  # Empty = no Python type in signature
       def type_check_expression(self, cpp_type, argument_var):
           return "isinstance(%s, (int, float))" % argument_var
       def input_conversion(self, cpp_type, argument_var, arg_num):
           # Convert Python value to C++ type
           code = Code().add("""
               |cdef _DPosition1 $dp
               |$dp[0] = <double>$argument_var
           """, locals())
           return code, dp, cleanup
   ```

2. **Register the converter** in `src/pyOpenMS/converters/__init__.py`:
   ```python
   special_converters.append(OpenMSDPosition1())
   ```

3. **Remove `wrap-ignore`** from methods using that type - autowrap will now handle them

### Autowrap Auto-generates `__init__` Dispatchers
- When a class has multiple constructors, autowrap generates `__init__` that dispatches to `_init_0`, `_init_1`, etc.
- **Don't add custom `__init__` in addons** unless you need special behavior (e.g., `String.pyx` uses `convString()`)
- Addon `__init__` **overwrites** autowrap's generated one - only do this intentionally

### Key Gotchas
- **autowrap methods return Python strings** - do NOT use `.decode('utf-8')` on return values from wrapped methods like `toString()`, `getSequence()`, etc.
- **`cdef` typed variables**: Use `cdef double`, `cdef int`, `cdef float` for C++ types, but NOT for strings returned by autowrap
- **Don't declare `cdef class ClassName:`** in addon files - methods are just indented code that gets injected into the generated class
- **Don't add Python-only methods to `.pxd` files** - declarations like `DataFrame get_df() except +` try to wrap non-existent C++ methods
- **Don't use `cdef` typed variables for autowrap return values** inside `def` methods - use Python type introspection instead (e.g., `isinstance(v, int)` rather than `cdef DataValue dv`)
- **Autowrap annotation format**: `# wrap-doc:` must be properly indented inside the class block (strict whitespace parsing)
- **Rebuilding after addon changes**: Delete `.cpp_extension_generated` stamp file to force regeneration:
  ```bash
  rm OpenMS-build/pyOpenMS/.cpp_extension_generated
  cmake --build OpenMS-build --target pyopenms -j4
  ```

### MSSpectrum Data Sources
- **FloatDataArrays**: Store per-peak float data like ion mobility values (use `getIMData()` to get the correct index)
- **StringDataArrays**: Store per-peak string data like ion annotations (look for name `'IonNames'`)
- **IntegerDataArrays**: Store per-peak integer data
- **MetaValues**: Store spectrum-level metadata (not per-peak)

### Testing pyOpenMS
```bash
PYTHONPATH=/path/to/OpenMS-build/pyOpenMS:$PYTHONPATH python3 -m pytest src/pyOpenMS/tests/ -v
```

### Common Patterns
- `__str__` should return simple/short output for user display
- `__repr__` should return detailed debug info with class name and key properties
- Use `self.size()` or wrap with `__len__` for container-like classes
- **Keep addons minimal** - only add methods that can't be auto-generated (like `__repr__`, `__str__`)
- **Don't add redundant aliases** - if `minX()` exists, don't add `getMin()` that just calls it

### Class-Specific Quirks

#### MetaInfoInterface inheritance
Not all classes inherit from MetaInfoInterface. Check before assuming meta value support:
- **Has meta values**: MSSpectrum, MSChromatogram, Feature, ConsensusFeature, PeptideHit, etc.
- **No meta values**: Mobilogram (no `getKeys()`, `getMetaValue()`, `metaValueExists()`)

#### Type requirements
- **PeptideIdentificationList**: `setPeptideIdentifications()` requires `PeptideIdentificationList`, not Python list:
  ```python
  pep_id_list = pyopenms.PeptideIdentificationList()
  pep_id_list.push_back(pep_id)
  cf.setPeptideIdentifications(pep_id_list)
  ```
- **AASequence validation**: `AASequence.fromString()` only accepts valid amino acid letters (A-Z except B, J, O, U, X, Z). Numbers and special chars fail.

#### Default/empty values
- `getDriftTimeUnitAsString()` returns `'<NONE>'` (not empty string) when no unit set
- Check for this literal string when testing for "no unit"

### Testing Quirks
- **Path shadowing**: Running pytest from source tree may pick up incomplete source `pyopenms/` instead of built version
- **Solution**: Run from `/tmp` or use explicit `PYTHONPATH=/path/to/OpenMS-build/pyOpenMS python3 -m pytest ...`
- for difficult wrapping issues also consult the autowrap implementation at https://github.com/OpenMS/autowrap
- pyopenms and openms can be build with multiple threads