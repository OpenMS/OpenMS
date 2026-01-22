# pyOpenMS Development

## Architecture

- **Autowrap** generates Cython bindings from `.pxd` files in `pxds/`
- **Addon files** (`.pyx`) in `addons/` inject custom Python methods into generated classes
- **DataFrame wrappers** in `pyopenms/_dataframes.py` add pandas methods (keeps pandas optional)
- **Type converters** in `converters/special_autowrap_conversionproviders.py` handle C++ to Python type conversion

## Addon Rules

- 4-space indent, no `cdef class` declaration - code is injected directly
- Don't add custom `__init__` unless necessary (overwrites autowrap's dispatcher)
- Don't add Python-only methods to `.pxd` files (tries to wrap non-existent C++)
- Don't use `cdef` typed variables for autowrap return values in `def` methods

## Type Handling

- Autowrap returns Python strings - never call `.decode('utf-8')` on them
- Use `cdef double/int/float` for C++ types, not for autowrap string returns
- `PeptideIdentificationList` required (not Python list) for `setPeptideIdentifications()`
- `AASequence.fromString()`: valid amino acids only (A-Z except B, J, O, U, X, Z)

## DataFrame Pattern

- `get_data_dict()` in Cython addon returns numpy arrays
- `get_df()` in `_dataframes.py` wraps with pandas
- Example: `_MSSpectrumDF` extends `_MSSpectrum`

## Type Converters

To make C++ types accept Python primitives (e.g., `DPosition<1>` from `float`):
1. Create converter in `converters/special_autowrap_conversionproviders.py`
2. Register in `converters/__init__.py`
3. Remove `wrap-ignore` from methods using that type

## Common Patterns

- `__str__`: short user display
- `__repr__`: detailed with class name and key properties
- Keep addons minimal - only add non-auto-generatable methods
- Don't add redundant aliases (if `minX()` exists, don't add `getMin()`)

## Data Sources (MSSpectrum)

- **FloatDataArrays**: per-peak float data (ion mobility via `getIMData()`)
- **StringDataArrays**: per-peak strings (ion annotations, name `'IonNames'`)
- **IntegerDataArrays**: per-peak integers
- **MetaValues**: spectrum-level metadata (not per-peak)

## Gotchas

- `getDriftTimeUnitAsString()` returns `'<NONE>'` not empty string
- Mobilogram has no MetaInfoInterface (no `getKeys()`, `getMetaValue()`)
- `# wrap-doc:` must be properly indented (strict whitespace parsing)
- Test from `/tmp` to avoid path shadowing with source `pyopenms/`

## Classes with MetaInfoInterface

**Has meta values**: MSSpectrum, MSChromatogram, Feature, ConsensusFeature, PeptideHit
**No meta values**: Mobilogram

## Reference

- Autowrap: https://github.com/OpenMS/autowrap
