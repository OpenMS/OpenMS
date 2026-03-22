# Arrow Schema Registry Design

**Date:** 2026-03-22
**Status:** Draft

## Problem

Arrow/Parquet schemas in OpenMS are defined inline across 10+ source files with ~117 unique column names as scattered string literals. There is no centralized schema registry, no validation on read or write, and no mechanism to ensure C++ and Python agree on schema structure. This creates risk for schema drift, inconsistent column naming, and silent compatibility breakage as the codebase evolves.

## Goals

In priority order:

1. **Developer consistency** — single source of truth for each schema; compile-time enforcement via constants; schema changes happen in one place
2. **Cross-language contract** — Python (pyOpenMS) and C++ share the same canonical schema definitions
3. **File compatibility** — validation at read and write boundaries catches schema mismatches (strict on write, subset on read)

## Non-Goals

- QPX file-level metadata (version, creator, UUID) — orthogonal, layered on later
- Bundle-level definitions (which tables belong together in a directory)
- Schema migration
- Refactoring IO class serialization logic — only schema definitions move

## Design

### File Structure

**New files:**

| File | Purpose |
|------|---------|
| `src/openms/include/OpenMS/FORMAT/ArrowSchemaRegistry.h` | All 18 schema structs + validation API (forward-declares Arrow types; implementations in .cpp) |
| `src/openms/source/FORMAT/ArrowSchemaRegistry.cpp` | Schema factory and validation implementation |
| `src/tests/class_tests/openms/source/ArrowSchemaRegistry_test.cpp` | Tests for all schemas and validation |

**Modified files:**

| File | Change |
|------|--------|
| `ProteinIdentificationArrowIO.cpp` | Replace inline schemas with registry calls, add validation |
| `FeatureMapArrowIO.cpp` | Replace inline schemas with registry calls, add validation |
| `ConsensusMapArrowIO.cpp` | Replace inline schemas with registry calls, add validation |
| `ConsensusMapArrowExport.cpp` | Replace inline schemas with registry calls, add validation |
| `QPXFile.cpp` | Replace inline schemas with registry calls, add validation |
| `MSExperimentArrowExport.cpp` | Replace inline schemas with registry calls, add validation |
| `TransitionParquetFile.cpp` | Replace inline schemas with registry calls, add validation |
| `OpenSwathOSWParquetWriter.cpp` | Replace inline schemas with registry calls, add validation |
| `MSChromatogramParquetConsumer.cpp` | Replace inline schemas with registry calls, add validation |
| `MobilogramParquetConsumer.cpp` | Replace inline schemas with registry calls, add validation |
| `arrow_zerocopy.cpp` | Expose schema structs to Python via nanobind (requires schema-only C Data Interface bridge + `ValidationResult` binding) |
| `msexperiment.py` (addon) | Replace hardcoded column inventory with registry constants |
| `consensusmap.py` (addon) | Replace hardcoded column inventory with registry constants |
| `peptideidentificationlist.py` (addon) | Replace hardcoded column inventory with registry constants |
| `sources.cmake` | Register new source files |
| `executables.cmake` | Register test executable |

All new files must be guarded with `#ifdef WITH_PARQUET` / `#endif`, consistent with existing Arrow/Parquet code. Schema structs with non-inline static methods require `OPENMS_DLLAPI` for Windows DLL builds.

### Schema Structs

Each table schema is a struct in `ArrowSchemaRegistry.h` containing:

1. **Column name constants** — `static constexpr const char*` for every field
2. **Nested type helpers** — static methods returning `std::shared_ptr<arrow::DataType>` for complex nested types (lists of structs, etc.)
3. **Schema factory** — `static std::shared_ptr<arrow::Schema> schema()` returning the canonical schema (the full superset of all possible columns)

Column names are scoped per schema (no shared namespace). Even when the string value is identical across schemas (e.g., `"rt"`), each schema owns its own constant independently.

Nullability for each field must match the existing inline definitions exactly (e.g., `unique_id` is non-null, `parent_feature_id` is nullable in FeatureSchema).

#### Example: FeatureSchema

```cpp
struct OPENMS_DLLAPI FeatureSchema
{
  // Column names
  static constexpr const char* UNIQUE_ID = "unique_id";
  static constexpr const char* PARENT_FEATURE_ID = "parent_feature_id";
  static constexpr const char* DEPTH = "depth";
  static constexpr const char* RT = "rt";
  static constexpr const char* MZ = "mz";
  static constexpr const char* INTENSITY = "intensity";
  static constexpr const char* CHARGE = "charge";
  static constexpr const char* QUALITY = "quality";
  static constexpr const char* QUALITY_RT = "quality_rt";
  static constexpr const char* QUALITY_MZ = "quality_mz";
  static constexpr const char* WIDTH = "width";
  static constexpr const char* RT_BB_MIN = "rt_bb_min";
  static constexpr const char* RT_BB_MAX = "rt_bb_max";
  static constexpr const char* MZ_BB_MIN = "mz_bb_min";
  static constexpr const char* MZ_BB_MAX = "mz_bb_max";
  static constexpr const char* CONVEX_HULLS = "convex_hulls";
  static constexpr const char* METAVALUES = "metavalues";

  /// Nested type: list<struct{hull_index:int32, points:list<struct{x:float64, y:float64}>}>
  static std::shared_ptr<arrow::DataType> convexHullType();

  /// Nested type: list<struct{name:utf8, value:utf8, value_type:utf8}>
  static std::shared_ptr<arrow::DataType> metavaluesType();

  /// Canonical Arrow schema for the features table (17 fields)
  static std::shared_ptr<arrow::Schema> schema();
};
```

#### Full Schema Inventory

**General format schemas:**

| Struct | Source Table | Columns | Write Validation | Read Validation |
|--------|-------------|---------|-----------------|-----------------|
| `ProteinSchema` | proteins.parquet | 10 | Strict | Subset |
| `ProteinGroupSchema` | protein_groups.parquet | 8 | Strict | Subset |
| `SearchParamsSchema` | search_params.parquet | 26 | Strict | Subset |
| `FeatureSchema` | features.parquet | 17 | Strict | Subset |
| `ConsensusFeatureSchema` | consensus_features.parquet | 9 | Strict | Subset |
| `ConsensusFeatureExportSchema` | consensus feature QPX export | 33 | Strict | Subset |
| `PSMSchema` | psms.parquet (QPX) | 25 | Strict | Subset |
| `SpectraLongSchema` | spectra (long format) | 12 (superset) | Subset | Subset |
| `SpectraSemiWideSchema` | spectra (semi-wide format) | 12 (superset) | Subset | Subset |
| `ChromatogramSchema` | chromatograms | 6 (superset) | Subset | Subset |

**OpenSWATH schemas:**

| Struct | Source Table | Columns | Write Validation | Read Validation |
|--------|-------------|---------|-----------------|-----------------|
| `OSWPrecursorSchema` | precursors.parquet | 10 | Strict | Subset |
| `OSWTransitionSchema` | transitions.parquet | 13 | Strict | Subset |
| `OSWFeatureSchema` | features.parquet | 65 | Strict | Subset |
| `OSWFeatureTransitionSchema` | feature_transition.parquet | 44 | Strict | Subset |
| `OSWFeaturePrecursorSchema` | feature_precursor.parquet | 5 | Strict | Subset |
| `OSWRunSchema` | runs.parquet | 2 | Strict | Subset |
| `XICSchema` | chromatogram .xic parquet | 18 | Strict | Subset |
| `XIMSchema` | mobilogram .xim parquet | 21 | Strict | Subset |

**Validation strategy:** Strict on write catches developer bugs (the primary goal). Subset on read ensures forward compatibility — current readers already accept partial tables (e.g., `ProteinIdentificationArrowIO` imports, PSM imports, OSW reader which discovers `ms1_`/`ms2_`/`var_` columns dynamically). Subset validation on read verifies that every column present has the correct name and type without requiring all columns to exist.

**Note:** `PSMSchema` defines the 25-column base PSM table from `QPXFile.cpp`. `FeatureMapArrowIO.cpp` and `ConsensusMapArrowIO.cpp` each add one extra ID column (`feature_unique_id` and `consensus_feature_unique_id` respectively) after calling `QPXFile::exportToArrow`. These extra columns are table-specific extensions and not part of `PSMSchema` itself.

#### Nested Type Helpers

Several schemas share structurally identical nested types. These are implemented as static methods on each schema that uses them, but the implementation can delegate to shared internal helpers to avoid duplication. The full inventory of nested types:

| Nested Type | Structure | Used By |
|-------------|-----------|---------|
| `metavaluesType()` | `list<struct{name:utf8, value:utf8, value_type:utf8}>` | FeatureSchema, ConsensusFeatureSchema, ConsensusFeatureExportSchema, PSMSchema, ProteinSchema, SearchParamsSchema |
| `modificationsType()` (protein) | `list<struct{position:int32, modification:utf8}>` | ProteinSchema |
| `modificationsType()` (PSM) | `list<struct{name:utf8, accession:utf8, positions:list<struct{position:utf8, scores:float64}>}>` | PSMSchema |
| `convexHullType()` | `list<struct{hull_index:int32, points:list<struct{x:float64, y:float64}>}>` | FeatureSchema |
| `handlesType()` | `list<struct{map_index:int64, unique_id:int64, rt:float64, mz:float64, intensity:float32, charge:int32, width:float32}>` | ConsensusFeatureSchema |
| `additionalScoresType()` | `list<struct{score_name:utf8, score_value:float64, higher_better:bool}>` | PSMSchema, ConsensusFeatureExportSchema |
| `intensitiesType()` | `list<struct{sample_accession:utf8, channel:utf8, intensity:float32}>` | ConsensusFeatureExportSchema |
| `floatDataType()` / `stringDataType()` / `integerDataType()` | `list<struct{name:utf8, values:list<T>}>` | ProteinGroupSchema |

### Validation API

```cpp
namespace OpenMS::ArrowSchemaValidation
{
  /// Validation mode
  enum class Mode
  {
    Strict,  ///< Exact match: all fields must be present with correct types, no extra fields
    Subset   ///< Every field present in the table must exist in the expected schema with correct type,
             ///< but not all expected fields need to be present. No unexpected fields allowed.
  };

  /// Result of schema validation
  struct OPENMS_DLLAPI ValidationResult
  {
    bool valid = true;
    std::vector<std::string> errors;

    /// Join all errors into a single string, separated by "; "
    std::string toString() const;
  };

  /// Validate an Arrow table's schema against the expected canonical schema.
  /// In Strict mode: checks field count, names, types (including nested), nullability, and order.
  /// In Subset mode: checks that every field present is a valid member of the expected schema
  /// with matching type and nullability. Missing fields are allowed; unexpected fields are errors.
  /// Schema-level metadata (key-value pairs) is ignored — only field structure is compared.
  /// Reports all mismatches (does not fail on first error).
  OPENMS_DLLAPI ValidationResult validate(
    const std::shared_ptr<arrow::Table>& table,
    const std::shared_ptr<arrow::Schema>& expected_schema,
    Mode mode = Mode::Strict);
}
```

**Error message examples:**

- `"Field count mismatch: got 15, expected 17"` (Strict only)
- `"Missing field 'convex_hulls'"` (Strict only)
- `"Unexpected field 'extra_column' not in expected schema"`
- `"Type mismatch for field 'mz': got int32, expected float64"`
- `"Nullability mismatch for field 'unique_id': got nullable, expected non-null"`

### Integration Points

#### Write Path

After constructing the Arrow table and before calling `parquet::arrow::WriteTable()`:

```cpp
// Fixed schema — strict validation
auto result = ArrowSchemaValidation::validate(table, FeatureSchema::schema());
if (!result.valid)
{
  throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
    "Schema validation failed: " + result.toString(), "");
}

// Configurable schema — subset validation
auto result = ArrowSchemaValidation::validate(table, SpectraLongSchema::schema(),
                                              ArrowSchemaValidation::Mode::Subset);
if (!result.valid)
{
  throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
    "Schema validation failed: " + result.toString(), "");
}
```

#### Read Path

After `parquet::arrow::ReadTable()` and before processing columns. Always uses Subset mode — current readers already accept partial tables, and the OSW reader discovers scoring columns dynamically by prefix:

```cpp
auto result = ArrowSchemaValidation::validate(table, FeatureSchema::schema(),
                                              ArrowSchemaValidation::Mode::Subset);
if (!result.valid)
{
  throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
    "Incompatible Parquet file schema: " + result.toString(), filename);
}
```

This catches type mismatches and unknown columns without rejecting files that are missing optional fields.

#### Column Access

Existing column access via string literals:

```cpp
// Before
auto col = table->GetColumnByName("rt");

// After
auto col = table->GetColumnByName(FeatureSchema::RT);
```

### Python Exposure

Schema structs are exposed via nanobind in the existing `arrow_zerocopy.cpp` module:

```python
import pyopenms

# Get canonical schema as pyarrow.Schema (via Arrow C Data Interface)
schema = pyopenms.FeatureSchema.schema()

# Access column name constants
col_name = pyopenms.FeatureSchema.RT  # "rt"

# Validate a pyarrow.Table
result = pyopenms.validate_arrow_schema(table, pyopenms.FeatureSchema.schema())
result.valid    # bool
result.errors   # list[str]

# Subset validation for configurable schemas
result = pyopenms.validate_arrow_schema(table, pyopenms.SpectraLongSchema.schema(), mode="subset")
```

The `schema()` method exports the `arrow::Schema` through the Arrow C Data Interface (same zero-copy path already used for table export). Column name constants are exposed as string attributes on the nanobind class.

### Testing Strategy

The test file `ArrowSchemaRegistry_test.cpp` covers:

1. **Schema construction** — each schema struct's `schema()` returns a non-null schema with the expected number of fields
2. **Column name consistency** — each constant matches the corresponding field name in the schema
3. **Nested type structure** — complex types (convex hulls, metavalues, modifications) have the correct nested structure
4. **Strict validation — valid table** — a table built from the canonical schema passes validation
5. **Strict validation — missing field** — removing a field produces the expected error
6. **Strict validation — extra field** — adding an unexpected field produces the expected error
7. **Strict validation — type mismatch** — changing a field type produces the expected error
8. **Strict validation — nullability mismatch** — changing nullability produces the expected error
9. **Subset validation — valid subset** — a table with a subset of columns passes subset validation
10. **Subset validation — unknown field** — a table with an unexpected column fails subset validation
11. **Subset validation — type mismatch** — a subset table with wrong type for a known column fails

Python-side tests in `test_arrow_zerocopy.py`:

1. **Schema round-trip** — `pyopenms.FeatureSchema.schema()` returns a `pyarrow.Schema` with correct fields
2. **Column constants** — `pyopenms.FeatureSchema.RT == "rt"`
3. **Strict validation from Python** — `validate_arrow_schema()` works on pyarrow tables
4. **Subset validation from Python** — subset mode works for configurable schemas
