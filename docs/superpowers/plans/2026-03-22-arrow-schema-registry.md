# Arrow Schema Registry Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Centralize all 18 Arrow/Parquet schema definitions into a single registry with validation, replacing scattered inline schema construction across 10+ source files.

**Architecture:** A single `ArrowSchemaRegistry.h/.cpp` containing 18 schema structs (each with column name constants, nested type helpers, and a schema factory). A `validate()` function checks tables against canonical schemas in Strict (write) or Subset (read) mode. Python gets access via nanobind in `arrow_zerocopy.cpp`.

**Tech Stack:** C++20, Apache Arrow, Parquet, nanobind, OpenMS test framework, pyarrow

**Spec:** `docs/superpowers/specs/2026-03-22-arrow-schema-registry-design.md`

---

### Task 1: Build system registration + skeleton files

**Files:**
- Create: `src/openms/include/OpenMS/FORMAT/ArrowSchemaRegistry.h`
- Create: `src/openms/source/FORMAT/ArrowSchemaRegistry.cpp`
- Create: `src/tests/class_tests/openms/source/ArrowSchemaRegistry_test.cpp`
- Modify: `src/openms/source/FORMAT/sources.cmake:109-122`
- Modify: `src/tests/class_tests/openms/executables.cmake:293-303`

- [ ] **Step 1: Create skeleton header**

Create `src/openms/include/OpenMS/FORMAT/ArrowSchemaRegistry.h` with:

```cpp
#pragma once

#ifdef WITH_PARQUET

#include <OpenMS/OpenMSConfig.h>
#include <memory>
#include <string>
#include <vector>

// Forward declarations
namespace arrow
{
  class Schema;
  class DataType;
  class Table;
}

namespace OpenMS
{

  namespace ArrowSchemaValidation
  {
    enum class Mode
    {
      Strict,
      Subset
    };

    struct OPENMS_DLLAPI ValidationResult
    {
      bool valid = true;
      std::vector<std::string> errors;
      std::string toString() const;
    };

    OPENMS_DLLAPI ValidationResult validate(
      const std::shared_ptr<arrow::Table>& table,
      const std::shared_ptr<arrow::Schema>& expected_schema,
      Mode mode = Mode::Strict);
  }

} // namespace OpenMS

#endif // WITH_PARQUET
```

- [ ] **Step 2: Create skeleton source**

Create `src/openms/source/FORMAT/ArrowSchemaRegistry.cpp` with:

```cpp
#ifdef WITH_PARQUET

#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <arrow/api.h>
#include <sstream>

namespace OpenMS
{

  namespace ArrowSchemaValidation
  {
    std::string ValidationResult::toString() const
    {
      std::ostringstream oss;
      for (size_t i = 0; i < errors.size(); ++i)
      {
        if (i > 0) oss << "; ";
        oss << errors[i];
      }
      return oss.str();
    }

    ValidationResult validate(
      const std::shared_ptr<arrow::Table>& table,
      const std::shared_ptr<arrow::Schema>& expected_schema,
      Mode mode)
    {
      ValidationResult result;
      // Implementation in Task 3
      return result;
    }
  }

} // namespace OpenMS

#endif // WITH_PARQUET
```

- [ ] **Step 3: Create skeleton test**

Create `src/tests/class_tests/openms/source/ArrowSchemaRegistry_test.cpp` with:

```cpp
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#ifdef WITH_PARQUET

#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <arrow/api.h>

using namespace OpenMS;

START_TEST(ArrowSchemaRegistry, "$Id$")

START_SECTION(ValidationResult::toString())
{
  ArrowSchemaValidation::ValidationResult result;
  result.valid = false;
  result.errors = {"error1", "error2"};
  TEST_STRING_EQUAL(result.toString(), "error1; error2")
}
END_SECTION

END_TEST

#else // WITH_PARQUET

START_TEST(ArrowSchemaRegistry, "$Id$")
END_TEST

#endif
```

- [ ] **Step 4: Register in build system**

In `src/openms/source/FORMAT/sources.cmake`, add after line 121 (`ConsensusMapArrowIO.cpp`):
```cmake
  list(APPEND sources_list ArrowSchemaRegistry.cpp)
```

In `src/tests/class_tests/openms/executables.cmake`, add `ArrowSchemaRegistry_test` to the list at line 294:
```cmake
  list(APPEND format_executables_list Arrow_test MSExperimentArrowExport_test ConsensusMapArrowExport_test QPXFile_test
    MSChromatogramParquetConsumer_test
    MobilogramParquetConsumer_test
    XICParquetFile_test
    XIMParquetFile_test
    OpenSwathOSWParquetRoundTrip_test
    ProteinIdentificationArrowIO_test
    FeatureMapArrowIO_test
    ConsensusMapArrowIO_test
    ArrowSchemaRegistry_test)
```

- [ ] **Step 5: Build and run skeleton test**

Run:
```bash
cd /home/sachsenb/Development/OpenMS && cmake --build OpenMS-build --target ArrowSchemaRegistry_test -j$(nproc)
```
Then:
```bash
cd /home/sachsenb/Development/OpenMS/OpenMS-build && ctest -R ArrowSchemaRegistry_test -V
```
Expected: PASS (skeleton test runs and passes)

- [ ] **Step 6: Commit**

```bash
git add src/openms/include/OpenMS/FORMAT/ArrowSchemaRegistry.h \
  src/openms/source/FORMAT/ArrowSchemaRegistry.cpp \
  src/tests/class_tests/openms/source/ArrowSchemaRegistry_test.cpp \
  src/openms/source/FORMAT/sources.cmake \
  src/tests/class_tests/openms/executables.cmake
git commit -m "feat: add ArrowSchemaRegistry skeleton with build system registration"
```

---

### Task 2: Validation logic (Strict + Subset)

**Files:**
- Modify: `src/openms/source/FORMAT/ArrowSchemaRegistry.cpp`
- Modify: `src/tests/class_tests/openms/source/ArrowSchemaRegistry_test.cpp`

- [ ] **Step 1: Write failing tests for validation**

Add to the test file before `END_TEST`:

```cpp
START_SECTION(validate — Strict mode — exact match passes)
{
  auto schema = arrow::schema({
    arrow::field("a", arrow::utf8(), false),
    arrow::field("b", arrow::float64(), true),
  });
  auto arr_a = arrow::ArrayFromJSON(arrow::utf8(), R"(["x","y"])");
  auto arr_b = arrow::ArrayFromJSON(arrow::float64(), "[1.0, 2.0]");
  auto table = arrow::Table::Make(schema, {arr_a, arr_b});

  auto result = ArrowSchemaValidation::validate(table, schema);
  TEST_EQUAL(result.valid, true)
  TEST_EQUAL(result.errors.size(), 0)
}
END_SECTION

START_SECTION(validate — Strict mode — missing field)
{
  auto expected = arrow::schema({
    arrow::field("a", arrow::utf8(), false),
    arrow::field("b", arrow::float64(), true),
  });
  auto actual_schema = arrow::schema({
    arrow::field("a", arrow::utf8(), false),
  });
  auto arr_a = arrow::ArrayFromJSON(arrow::utf8(), R"(["x"])");
  auto table = arrow::Table::Make(actual_schema, {arr_a});

  auto result = ArrowSchemaValidation::validate(table, expected);
  TEST_EQUAL(result.valid, false)
  TEST_EQUAL(result.errors.empty(), false)
}
END_SECTION

START_SECTION(validate — Strict mode — extra field)
{
  auto expected = arrow::schema({
    arrow::field("a", arrow::utf8(), false),
  });
  auto actual_schema = arrow::schema({
    arrow::field("a", arrow::utf8(), false),
    arrow::field("b", arrow::float64(), true),
  });
  auto arr_a = arrow::ArrayFromJSON(arrow::utf8(), R"(["x"])");
  auto arr_b = arrow::ArrayFromJSON(arrow::float64(), "[1.0]");
  auto table = arrow::Table::Make(actual_schema, {arr_a, arr_b});

  auto result = ArrowSchemaValidation::validate(table, expected);
  TEST_EQUAL(result.valid, false)
}
END_SECTION

START_SECTION(validate — Strict mode — type mismatch)
{
  auto expected = arrow::schema({
    arrow::field("a", arrow::float64(), false),
  });
  auto actual_schema = arrow::schema({
    arrow::field("a", arrow::int32(), false),
  });
  auto arr_a = arrow::ArrayFromJSON(arrow::int32(), "[1]");
  auto table = arrow::Table::Make(actual_schema, {arr_a});

  auto result = ArrowSchemaValidation::validate(table, expected);
  TEST_EQUAL(result.valid, false)
}
END_SECTION

START_SECTION(validate — Strict mode — nullability mismatch)
{
  auto expected = arrow::schema({
    arrow::field("a", arrow::utf8(), false),
  });
  auto actual_schema = arrow::schema({
    arrow::field("a", arrow::utf8(), true),
  });
  auto arr_a = arrow::ArrayFromJSON(arrow::utf8(), R"(["x"])");
  auto table = arrow::Table::Make(actual_schema, {arr_a});

  auto result = ArrowSchemaValidation::validate(table, expected);
  TEST_EQUAL(result.valid, false)
}
END_SECTION

START_SECTION(validate — Subset mode — valid subset passes)
{
  auto expected = arrow::schema({
    arrow::field("a", arrow::utf8(), false),
    arrow::field("b", arrow::float64(), true),
    arrow::field("c", arrow::int32(), true),
  });
  auto actual_schema = arrow::schema({
    arrow::field("a", arrow::utf8(), false),
    arrow::field("c", arrow::int32(), true),
  });
  auto arr_a = arrow::ArrayFromJSON(arrow::utf8(), R"(["x"])");
  auto arr_c = arrow::ArrayFromJSON(arrow::int32(), "[1]");
  auto table = arrow::Table::Make(actual_schema, {arr_a, arr_c});

  auto result = ArrowSchemaValidation::validate(table, expected,
    ArrowSchemaValidation::Mode::Subset);
  TEST_EQUAL(result.valid, true)
}
END_SECTION

START_SECTION(validate — Subset mode — unknown field fails)
{
  auto expected = arrow::schema({
    arrow::field("a", arrow::utf8(), false),
  });
  auto actual_schema = arrow::schema({
    arrow::field("a", arrow::utf8(), false),
    arrow::field("z", arrow::float64(), true),
  });
  auto arr_a = arrow::ArrayFromJSON(arrow::utf8(), R"(["x"])");
  auto arr_z = arrow::ArrayFromJSON(arrow::float64(), "[1.0]");
  auto table = arrow::Table::Make(actual_schema, {arr_a, arr_z});

  auto result = ArrowSchemaValidation::validate(table, expected,
    ArrowSchemaValidation::Mode::Subset);
  TEST_EQUAL(result.valid, false)
}
END_SECTION

START_SECTION(validate — Subset mode — type mismatch fails)
{
  auto expected = arrow::schema({
    arrow::field("a", arrow::float64(), false),
  });
  auto actual_schema = arrow::schema({
    arrow::field("a", arrow::int32(), false),
  });
  auto arr_a = arrow::ArrayFromJSON(arrow::int32(), "[1]");
  auto table = arrow::Table::Make(actual_schema, {arr_a});

  auto result = ArrowSchemaValidation::validate(table, expected,
    ArrowSchemaValidation::Mode::Subset);
  TEST_EQUAL(result.valid, false)
}
END_SECTION

START_SECTION(validate — metadata is ignored)
{
  auto metadata = arrow::key_value_metadata({"key"}, {"value"});
  auto expected = arrow::schema({
    arrow::field("a", arrow::utf8(), false),
  });
  auto actual_schema = arrow::schema({
    arrow::field("a", arrow::utf8(), false),
  })->WithMetadata(metadata);
  auto arr_a = arrow::ArrayFromJSON(arrow::utf8(), R"(["x"])");
  auto table = arrow::Table::Make(actual_schema, {arr_a});

  auto result = ArrowSchemaValidation::validate(table, expected);
  TEST_EQUAL(result.valid, true)
}
END_SECTION
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /home/sachsenb/Development/OpenMS && cmake --build OpenMS-build --target ArrowSchemaRegistry_test -j$(nproc) && cd OpenMS-build && ctest -R ArrowSchemaRegistry_test -V`

Expected: Tests compile but several fail (validate always returns valid=true)

- [ ] **Step 3: Implement validate()**

Replace the placeholder `validate()` body in `ArrowSchemaRegistry.cpp`:

```cpp
    ValidationResult validate(
      const std::shared_ptr<arrow::Table>& table,
      const std::shared_ptr<arrow::Schema>& expected_schema,
      Mode mode)
    {
      ValidationResult result;
      auto actual_schema = table->schema();
      int actual_n = actual_schema->num_fields();
      int expected_n = expected_schema->num_fields();

      if (mode == Mode::Strict)
      {
        if (actual_n != expected_n)
        {
          result.valid = false;
          result.errors.push_back("Field count mismatch: got " + std::to_string(actual_n) +
            ", expected " + std::to_string(expected_n));
        }

        for (int i = 0; i < std::min(actual_n, expected_n); ++i)
        {
          auto actual_field = actual_schema->field(i);
          auto expected_field = expected_schema->field(i);

          if (actual_field->name() != expected_field->name())
          {
            result.valid = false;
            result.errors.push_back("Field name mismatch at index " + std::to_string(i) +
              ": got '" + actual_field->name() + "', expected '" + expected_field->name() + "'");
            continue;
          }
          if (!actual_field->type()->Equals(expected_field->type()))
          {
            result.valid = false;
            result.errors.push_back("Type mismatch for field '" + actual_field->name() +
              "': got " + actual_field->type()->ToString() + ", expected " + expected_field->type()->ToString());
          }
          if (actual_field->nullable() != expected_field->nullable())
          {
            result.valid = false;
            result.errors.push_back("Nullability mismatch for field '" + actual_field->name() +
              "': got " + std::string(actual_field->nullable() ? "nullable" : "non-null") +
              ", expected " + std::string(expected_field->nullable() ? "nullable" : "non-null"));
          }
        }

        // Check for extra fields beyond expected count
        for (int i = expected_n; i < actual_n; ++i)
        {
          result.valid = false;
          result.errors.push_back("Unexpected field '" + actual_schema->field(i)->name() +
            "' not in expected schema");
        }

        // Check for missing fields beyond actual count
        for (int i = actual_n; i < expected_n; ++i)
        {
          result.valid = false;
          result.errors.push_back("Missing field '" + expected_schema->field(i)->name() + "'");
        }
      }
      else // Mode::Subset
      {
        for (int i = 0; i < actual_n; ++i)
        {
          auto actual_field = actual_schema->field(i);
          auto expected_field = expected_schema->GetFieldByName(actual_field->name());

          if (!expected_field)
          {
            result.valid = false;
            result.errors.push_back("Unexpected field '" + actual_field->name() +
              "' not in expected schema");
            continue;
          }
          if (!actual_field->type()->Equals(expected_field->type()))
          {
            result.valid = false;
            result.errors.push_back("Type mismatch for field '" + actual_field->name() +
              "': got " + actual_field->type()->ToString() + ", expected " + expected_field->type()->ToString());
          }
          if (actual_field->nullable() != expected_field->nullable())
          {
            result.valid = false;
            result.errors.push_back("Nullability mismatch for field '" + actual_field->name() +
              "': got " + std::string(actual_field->nullable() ? "nullable" : "non-null") +
              ", expected " + std::string(expected_field->nullable() ? "nullable" : "non-null"));
          }
        }
      }

      return result;
    }
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /home/sachsenb/Development/OpenMS && cmake --build OpenMS-build --target ArrowSchemaRegistry_test -j$(nproc) && cd OpenMS-build && ctest -R ArrowSchemaRegistry_test -V`

Expected: All tests PASS

- [ ] **Step 5: Commit**

```bash
git add src/openms/source/FORMAT/ArrowSchemaRegistry.cpp \
  src/tests/class_tests/openms/source/ArrowSchemaRegistry_test.cpp
git commit -m "feat: implement ArrowSchemaValidation::validate with Strict and Subset modes"
```

---

### Task 3: General format schema structs (ProteinSchema, ProteinGroupSchema, SearchParamsSchema)

**Files:**
- Modify: `src/openms/include/OpenMS/FORMAT/ArrowSchemaRegistry.h`
- Modify: `src/openms/source/FORMAT/ArrowSchemaRegistry.cpp`
- Modify: `src/tests/class_tests/openms/source/ArrowSchemaRegistry_test.cpp`

**IMPORTANT — Nullability:** `ProteinIdentificationArrowIO.cpp` uses explicit `/*nullable=*/true` or `/*nullable=*/false` per field. Copy these exactly. Other files (ConsensusMapArrowIO, QPXFile, OpenSWATH files) omit the nullable parameter, defaulting to `nullable=true` for all fields. Always match the source file's behavior.

**Reference files (read only — extract exact field definitions from these):**
- `src/openms/source/FORMAT/ProteinIdentificationArrowIO.cpp:593-604` (ProteinSchema — 10 fields)
- `src/openms/source/FORMAT/ProteinIdentificationArrowIO.cpp:819-828` (ProteinGroupSchema — 8 fields)
- `src/openms/source/FORMAT/ProteinIdentificationArrowIO.cpp:1241-1268` (SearchParamsSchema — 26 fields)

Nested types to extract from the same file:
- `metavaluesType()`: `list<struct{name:utf8, value:utf8, value_type:utf8}>` (used by ProteinSchema and SearchParamsSchema; ProteinGroupSchema does NOT have metavalues)
- `modificationsType()`: `list<struct{position:int32, modification:utf8}>` (ProteinSchema only)
- `floatDataType()`/`stringDataType()`/`integerDataType()`: `list<struct{name:utf8, values:list<T>}>` (ProteinGroupSchema only, defined around lines 645-682)

- [ ] **Step 1: Write failing tests for ProteinSchema**

Add tests that:
- `ProteinSchema::schema()` returns non-null with 10 fields
- `ProteinSchema::ACCESSION` equals `"accession"` (and similarly for all 10 column constants)
- A table built from `ProteinSchema::schema()` passes strict validation

- [ ] **Step 2: Run tests to verify they fail** (compile error — struct not defined yet)

- [ ] **Step 3: Add ProteinSchema struct to header and implement schema() in .cpp**

Read `ProteinIdentificationArrowIO.cpp:411-423` for the modification struct type and `ProteinIdentificationArrowIO.cpp:593-604` for the exact field definitions. Copy field names, types, and nullability exactly.

Header adds:
```cpp
struct OPENMS_DLLAPI ProteinSchema
{
  static constexpr const char* ACCESSION = "accession";
  static constexpr const char* SCORE = "score";
  static constexpr const char* RANK = "rank";
  static constexpr const char* COVERAGE = "coverage";
  static constexpr const char* SEQUENCE = "sequence";
  static constexpr const char* DESCRIPTION = "description";
  static constexpr const char* IS_DECOY = "is_decoy";
  static constexpr const char* RUN_IDENTIFIER = "run_identifier";
  static constexpr const char* MODIFICATIONS = "modifications";
  static constexpr const char* METAVALUES = "metavalues";

  static std::shared_ptr<arrow::DataType> modificationsType();
  static std::shared_ptr<arrow::DataType> metavaluesType();
  static std::shared_ptr<arrow::Schema> schema();
};
```

.cpp implements `schema()` matching the inline definition at line 593-604 exactly.

- [ ] **Step 4: Run tests to verify they pass**

- [ ] **Step 5: Repeat steps 1-4 for ProteinGroupSchema** (8 fields, read lines 819-828 and nested float/string/integer data types from lines 645-682)

- [ ] **Step 6: Repeat steps 1-4 for SearchParamsSchema** (26 fields, read lines 1241-1268)

- [ ] **Step 7: Build and run full test suite**

Run: `cmake --build OpenMS-build --target ArrowSchemaRegistry_test -j$(nproc) && cd OpenMS-build && ctest -R ArrowSchemaRegistry_test -V`

Expected: All tests PASS

- [ ] **Step 8: Commit**

```bash
git add src/openms/include/OpenMS/FORMAT/ArrowSchemaRegistry.h \
  src/openms/source/FORMAT/ArrowSchemaRegistry.cpp \
  src/tests/class_tests/openms/source/ArrowSchemaRegistry_test.cpp
git commit -m "feat: add ProteinSchema, ProteinGroupSchema, SearchParamsSchema to registry"
```

---

### Task 4: General format schema structs (FeatureSchema, ConsensusFeatureSchema, PSMSchema)

**Files:**
- Modify: `src/openms/include/OpenMS/FORMAT/ArrowSchemaRegistry.h`
- Modify: `src/openms/source/FORMAT/ArrowSchemaRegistry.cpp`
- Modify: `src/tests/class_tests/openms/source/ArrowSchemaRegistry_test.cpp`

**IMPORTANT — Nullability:** `FeatureMapArrowIO.cpp` uses explicit `/*nullable=*/true` or `/*nullable=*/false` per field — copy exactly. `ConsensusMapArrowIO.cpp` and `QPXFile.cpp` omit the nullable parameter on all fields, defaulting to `nullable=true`. Do NOT add explicit `nullable=false` to schemas that don't have it in the source.

**Reference files (read only):**
- `src/openms/source/FORMAT/FeatureMapArrowIO.cpp:829-847` (convexHullType), `1036-1054` (FeatureSchema — 17 fields)
- `src/openms/source/FORMAT/ConsensusMapArrowIO.cpp:1008-1025` (handlesType), `1130-1140` (ConsensusFeatureSchema — 9 fields)
- `src/openms/source/FORMAT/QPXFile.cpp:74-98` (PSM nested types), `644-670` (PSMSchema — 25 fields)

- [ ] **Step 1: Write failing tests for FeatureSchema** (17 fields, convexHullType, metavaluesType)

- [ ] **Step 2: Implement FeatureSchema** — read exact fields from `FeatureMapArrowIO.cpp:1036-1054`

- [ ] **Step 3: Run tests, verify pass**

- [ ] **Step 4: Write failing tests for ConsensusFeatureSchema** (9 fields, handlesType)

- [ ] **Step 5: Implement ConsensusFeatureSchema** — read exact fields from `ConsensusMapArrowIO.cpp:1130-1140` and handles struct from `1008-1025`

- [ ] **Step 6: Run tests, verify pass**

- [ ] **Step 7: Write failing tests for PSMSchema** (25 fields, modificationsType, additionalScoresType)

- [ ] **Step 8: Implement PSMSchema** — read exact fields from `QPXFile.cpp:644-670` and nested types from `74-98`

- [ ] **Step 9: Run tests, verify pass**

- [ ] **Step 10: Commit**

```bash
git add src/openms/include/OpenMS/FORMAT/ArrowSchemaRegistry.h \
  src/openms/source/FORMAT/ArrowSchemaRegistry.cpp \
  src/tests/class_tests/openms/source/ArrowSchemaRegistry_test.cpp
git commit -m "feat: add FeatureSchema, ConsensusFeatureSchema, PSMSchema to registry"
```

---

### Task 5: Export format schema structs (ConsensusFeatureExportSchema, Spectra, Chromatogram)

**Files:**
- Modify: `src/openms/include/OpenMS/FORMAT/ArrowSchemaRegistry.h`
- Modify: `src/openms/source/FORMAT/ArrowSchemaRegistry.cpp`
- Modify: `src/tests/class_tests/openms/source/ArrowSchemaRegistry_test.cpp`

**Reference files (read only):**
- `src/openms/source/FORMAT/ConsensusMapArrowExport.cpp:77-141` (nested types), `789-824` (ConsensusFeatureExportSchema — 33 fields)
- `src/openms/source/FORMAT/MSExperimentArrowExport.cpp:358-454` (SpectraLongSchema — 12 superset)
- `src/openms/source/FORMAT/MSExperimentArrowExport.cpp:646-738` (SpectraSemiWideSchema — 12 superset)
- `src/openms/source/FORMAT/MSExperimentArrowExport.cpp:907-912,988-995` (ChromatogramSchema — 6 superset)

- [ ] **Step 1: Write failing tests for ConsensusFeatureExportSchema** (33 fields, intensitiesType)

- [ ] **Step 2: Implement ConsensusFeatureExportSchema** — read exact fields from `ConsensusMapArrowExport.cpp:789-824`

- [ ] **Step 3: Run tests, verify pass**

- [ ] **Step 4: Write failing tests for SpectraLongSchema** (12 fields superset). Test with subset validation since this is a configurable schema.

- [ ] **Step 5: Implement SpectraLongSchema** — read exact fields from `MSExperimentArrowExport.cpp:358-454`

- [ ] **Step 6: Run tests, verify pass**

- [ ] **Step 7: Write failing tests for SpectraSemiWideSchema** (12 fields superset, note list types for mz/intensity/ion_mobility)

- [ ] **Step 8: Implement SpectraSemiWideSchema** — read exact fields from `MSExperimentArrowExport.cpp:646-738`

- [ ] **Step 9: Run tests, verify pass**

- [ ] **Step 10: Write failing tests for ChromatogramSchema** (6 fields superset)

- [ ] **Step 11: Implement ChromatogramSchema** — read exact fields from `MSExperimentArrowExport.cpp:907-912,988-995`

- [ ] **Step 12: Run tests, verify pass**

- [ ] **Step 13: Commit**

```bash
git add src/openms/include/OpenMS/FORMAT/ArrowSchemaRegistry.h \
  src/openms/source/FORMAT/ArrowSchemaRegistry.cpp \
  src/tests/class_tests/openms/source/ArrowSchemaRegistry_test.cpp
git commit -m "feat: add ConsensusFeatureExportSchema, SpectraLongSchema, SpectraSemiWideSchema, ChromatogramSchema"
```

---

### Task 6: OpenSWATH schema structs (OSWPrecursor, OSWTransition, OSWRun, OSWFeaturePrecursor)

**Files:**
- Modify: `src/openms/include/OpenMS/FORMAT/ArrowSchemaRegistry.h`
- Modify: `src/openms/source/FORMAT/ArrowSchemaRegistry.cpp`
- Modify: `src/tests/class_tests/openms/source/ArrowSchemaRegistry_test.cpp`

**IMPORTANT — Nullability:** All OpenSWATH source files omit the nullable parameter, defaulting to `nullable=true` for all fields. Do NOT add explicit `nullable=false`.

**Reference files (read only):**
- `src/openms/source/ANALYSIS/OPENSWATH/TransitionParquetFile.cpp:556-566` (OSWPrecursorSchema — 10 fields)
- `src/openms/source/ANALYSIS/OPENSWATH/TransitionParquetFile.cpp:655-668` (OSWTransitionSchema — 13 fields)
- `src/openms/source/ANALYSIS/OPENSWATH/OpenSwathOSWParquetWriter.cpp:1307-1312` (OSWFeaturePrecursorSchema — 5 fields)
- `src/openms/source/ANALYSIS/OPENSWATH/OpenSwathOSWParquetWriter.cpp:1351-1353` (OSWRunSchema — 2 fields)

- [ ] **Step 1: Write failing tests for OSWPrecursorSchema, OSWTransitionSchema, OSWFeaturePrecursorSchema, OSWRunSchema**

These are small schemas (2-13 fields, no complex nested types). Test field count and column name constants.

- [ ] **Step 2: Implement all four schemas** — read exact fields from the reference files above

- [ ] **Step 3: Run tests, verify pass**

- [ ] **Step 4: Commit**

```bash
git add src/openms/include/OpenMS/FORMAT/ArrowSchemaRegistry.h \
  src/openms/source/FORMAT/ArrowSchemaRegistry.cpp \
  src/tests/class_tests/openms/source/ArrowSchemaRegistry_test.cpp
git commit -m "feat: add OSWPrecursorSchema, OSWTransitionSchema, OSWFeaturePrecursorSchema, OSWRunSchema"
```

---

### Task 7: OpenSWATH schema structs (OSWFeature, OSWFeatureTransition, XIC, XIM)

**Files:**
- Modify: `src/openms/include/OpenMS/FORMAT/ArrowSchemaRegistry.h`
- Modify: `src/openms/source/FORMAT/ArrowSchemaRegistry.cpp`
- Modify: `src/tests/class_tests/openms/source/ArrowSchemaRegistry_test.cpp`

**Reference files (read only):**
- `src/openms/source/ANALYSIS/OPENSWATH/OpenSwathOSWParquetWriter.cpp:1066-1131` (OSWFeatureSchema — 65 fields)
- `src/openms/source/ANALYSIS/OPENSWATH/OpenSwathOSWParquetWriter.cpp:1206-1250` (OSWFeatureTransitionSchema — 44 fields)
- `src/openms/source/FORMAT/DATAACCESS/MSChromatogramParquetConsumer.cpp:195-213` (XICSchema — 18 fields)
- `src/openms/source/FORMAT/DATAACCESS/MobilogramParquetConsumer.cpp:493-514` (XIMSchema — 21 fields)

- [ ] **Step 1: Write failing tests for OSWFeatureSchema** (65 fields — test field count and a representative sample of column constants)

- [ ] **Step 2: Implement OSWFeatureSchema** — read exact fields from `OpenSwathOSWParquetWriter.cpp:1066-1131`. This is the largest schema; take care to copy all 65 field names, types, and nullabilities.

- [ ] **Step 3: Run tests, verify pass**

- [ ] **Step 4: Write failing tests for OSWFeatureTransitionSchema** (44 fields)

- [ ] **Step 5: Implement OSWFeatureTransitionSchema** — read exact fields from `OpenSwathOSWParquetWriter.cpp:1206-1250`

- [ ] **Step 6: Run tests, verify pass**

- [ ] **Step 7: Write failing tests for XICSchema** (18 fields)

- [ ] **Step 8: Implement XICSchema** — read exact fields from `MSChromatogramParquetConsumer.cpp:195-213`

- [ ] **Step 9: Run tests, verify pass**

- [ ] **Step 10: Write failing tests for XIMSchema** (21 fields)

- [ ] **Step 11: Implement XIMSchema** — read exact fields from `MobilogramParquetConsumer.cpp:493-514`

- [ ] **Step 12: Run tests, verify pass**

- [ ] **Step 13: Commit**

```bash
git add src/openms/include/OpenMS/FORMAT/ArrowSchemaRegistry.h \
  src/openms/source/FORMAT/ArrowSchemaRegistry.cpp \
  src/tests/class_tests/openms/source/ArrowSchemaRegistry_test.cpp
git commit -m "feat: add OSWFeatureSchema, OSWFeatureTransitionSchema, XICSchema, XIMSchema"
```

---

### Task 8: Integrate registry into ProteinIdentificationArrowIO

**Files:**
- Modify: `src/openms/source/FORMAT/ProteinIdentificationArrowIO.cpp`

**Strategy:** Replace inline schema construction with registry calls, replace string literals with constants, add validation on write and read paths. This is the first integration — do it carefully as a pattern for subsequent tasks.

- [ ] **Step 1: Run existing ProteinIdentificationArrowIO tests to establish baseline**

Run: `cd OpenMS-build && ctest -R ProteinIdentificationArrowIO_test -V`

Expected: All tests PASS (baseline)

- [ ] **Step 2: Add `#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>` to ProteinIdentificationArrowIO.cpp**

- [ ] **Step 3: Replace proteins table inline schema (lines 593-604)**

Replace the `arrow::schema({...})` block with `ProteinSchema::schema()`. The nested type builders (`mod_struct_type`, `metavalues_builder.type()`) are still used by the array builders — only the schema construction changes.

Before:
```cpp
auto schema = arrow::schema({
  arrow::field("accession", arrow::utf8(), /*nullable=*/false),
  // ... 10 fields
});
```

After:
```cpp
auto schema = ProteinSchema::schema();
```

- [ ] **Step 4: Replace string literals in column access on the read path**

Find all `getColumn_(table, "accession")` style calls for proteins import (around lines 1484-1493) and replace with `getColumn_(table, ProteinSchema::ACCESSION)`. Do the same for all 10 protein columns.

- [ ] **Step 5: Add write-path validation** — add `ArrowSchemaValidation::validate(table, ProteinSchema::schema())` check after table construction (around line 606), before `WriteTable`.

- [ ] **Step 6: Add read-path validation** — add `ArrowSchemaValidation::validate(table, ProteinSchema::schema(), ArrowSchemaValidation::Mode::Subset)` check after `ReadTable` in the protein import function.

- [ ] **Step 7: Repeat steps 3-6 for protein groups table** (lines 819-828 schema, read path around line 1611+)

- [ ] **Step 8: Repeat steps 3-6 for search params table** (lines 1241-1268 schema, read path around line 1332+)

- [ ] **Step 9: Run existing tests to verify nothing broke**

Run: `cd OpenMS-build && ctest -R ProteinIdentificationArrowIO_test -V`

Expected: All tests still PASS

- [ ] **Step 10: Commit**

```bash
git add src/openms/source/FORMAT/ProteinIdentificationArrowIO.cpp
git commit -m "refactor: use ArrowSchemaRegistry in ProteinIdentificationArrowIO"
```

---

### Task 9: Integrate registry into FeatureMapArrowIO and ConsensusMapArrowIO

**Files:**
- Modify: `src/openms/source/FORMAT/FeatureMapArrowIO.cpp`
- Modify: `src/openms/source/FORMAT/ConsensusMapArrowIO.cpp`

Follow the same pattern as Task 8:

- [ ] **Step 1: Run baseline tests**

Run: `cd OpenMS-build && ctest -R "FeatureMapArrowIO_test|ConsensusMapArrowIO_test" -V`

- [ ] **Step 2: Integrate FeatureSchema into FeatureMapArrowIO.cpp**

Replace inline schema at `FeatureMapArrowIO.cpp:1036-1054` with `FeatureSchema::schema()`. Replace string literals in column access. Add write/read validation.

Note: this file also writes PSMs via `QPXFile::exportToArrow` — do NOT add PSM validation here (that belongs in QPXFile.cpp). Only add `FeatureSchema` validation for the features table.

- [ ] **Step 3: Run FeatureMapArrowIO tests**

- [ ] **Step 4: Integrate ConsensusFeatureSchema into ConsensusMapArrowIO.cpp**

Replace inline schema at `ConsensusMapArrowIO.cpp:1130-1140` with `ConsensusFeatureSchema::schema()`. Replace string literals. Add validation.

- [ ] **Step 5: Run ConsensusMapArrowIO tests**

- [ ] **Step 6: Commit**

```bash
git add src/openms/source/FORMAT/FeatureMapArrowIO.cpp \
  src/openms/source/FORMAT/ConsensusMapArrowIO.cpp
git commit -m "refactor: use ArrowSchemaRegistry in FeatureMapArrowIO and ConsensusMapArrowIO"
```

---

### Task 10: Integrate registry into QPXFile, ConsensusMapArrowExport, MSExperimentArrowExport

**Files:**
- Modify: `src/openms/source/FORMAT/QPXFile.cpp`
- Modify: `src/openms/source/FORMAT/ConsensusMapArrowExport.cpp`
- Modify: `src/openms/source/FORMAT/MSExperimentArrowExport.cpp`

- [ ] **Step 1: Run baseline tests**

Run: `cd OpenMS-build && ctest -R "QPXFile_test|ConsensusMapArrowExport_test|MSExperimentArrowExport_test" -V`

- [ ] **Step 2: Integrate PSMSchema into QPXFile.cpp**

Replace inline schema at `QPXFile.cpp:644-670` with `PSMSchema::schema()`. Replace string literals. Add strict write validation, subset read validation.

- [ ] **Step 3: Run QPXFile tests**

- [ ] **Step 4: Integrate ConsensusFeatureExportSchema into ConsensusMapArrowExport.cpp**

Replace inline schema at `ConsensusMapArrowExport.cpp:789-824` with `ConsensusFeatureExportSchema::schema()`. Add strict write validation.

- [ ] **Step 5: Run ConsensusMapArrowExport tests**

- [ ] **Step 6: Integrate SpectraLongSchema, SpectraSemiWideSchema, ChromatogramSchema into MSExperimentArrowExport.cpp**

These are configurable schemas — use **Subset validation on both write and read**. Replace column name string literals with schema constants. The dynamic column selection logic stays as-is; validation just checks that whatever columns are included are valid members of the superset.

- [ ] **Step 7: Run MSExperimentArrowExport tests**

- [ ] **Step 8: Commit**

```bash
git add src/openms/source/FORMAT/QPXFile.cpp \
  src/openms/source/FORMAT/ConsensusMapArrowExport.cpp \
  src/openms/source/FORMAT/MSExperimentArrowExport.cpp
git commit -m "refactor: use ArrowSchemaRegistry in QPXFile, ConsensusMapArrowExport, MSExperimentArrowExport"
```

---

### Task 11: Integrate registry into OpenSWATH files

**Files:**
- Modify: `src/openms/source/ANALYSIS/OPENSWATH/TransitionParquetFile.cpp`
- Modify: `src/openms/source/ANALYSIS/OPENSWATH/OpenSwathOSWParquetWriter.cpp`
- Modify: `src/openms/source/FORMAT/DATAACCESS/MSChromatogramParquetConsumer.cpp`
- Modify: `src/openms/source/FORMAT/DATAACCESS/MobilogramParquetConsumer.cpp`

- [ ] **Step 1: Run baseline tests**

Run: `cd OpenMS-build && ctest -R "OpenSwathOSWParquetRoundTrip_test|MSChromatogramParquetConsumer_test|MobilogramParquetConsumer_test" -V`

- [ ] **Step 2: Integrate OSWPrecursorSchema and OSWTransitionSchema into TransitionParquetFile.cpp**

Replace inline schemas at lines 556-566 and 655-668. Replace string literals. Add validation.

- [ ] **Step 3: Integrate OSWFeatureSchema, OSWFeatureTransitionSchema, OSWFeaturePrecursorSchema, OSWRunSchema into OpenSwathOSWParquetWriter.cpp**

Replace inline schemas at lines 1066-1131, 1206-1250, 1307-1312, 1351-1353. Add strict write validation.

- [ ] **Step 4: Integrate XICSchema into MSChromatogramParquetConsumer.cpp**

Replace inline schema at lines 195-213. Add validation.

- [ ] **Step 5: Integrate XIMSchema into MobilogramParquetConsumer.cpp**

Replace inline schema at lines 493-514. Add validation.

- [ ] **Step 6: Run all OpenSWATH tests**

Run: `cd OpenMS-build && ctest -R "OpenSwathOSWParquetRoundTrip_test|MSChromatogramParquetConsumer_test|MobilogramParquetConsumer_test|XICParquetFile_test|XIMParquetFile_test" -V`

Expected: All tests PASS

- [ ] **Step 7: Commit**

```bash
git add src/openms/source/ANALYSIS/OPENSWATH/TransitionParquetFile.cpp \
  src/openms/source/ANALYSIS/OPENSWATH/OpenSwathOSWParquetWriter.cpp \
  src/openms/source/FORMAT/DATAACCESS/MSChromatogramParquetConsumer.cpp \
  src/openms/source/FORMAT/DATAACCESS/MobilogramParquetConsumer.cpp
git commit -m "refactor: use ArrowSchemaRegistry in OpenSWATH and chromatogram/mobilogram consumers"
```

---

### Task 12: Python exposure via nanobind

**Files:**
- Modify: `src/pyOpenMS/bindings/arrow_zerocopy.cpp`
- Modify: `src/pyOpenMS/tests/unittests/test_arrow_zerocopy.py`

**Reference:** `src/pyOpenMS/bindings/arrow_zerocopy.cpp:157-164` (module definition)

- [ ] **Step 1: Add schema-to-pyarrow and pyarrow-to-schema bridge functions**

In `arrow_zerocopy.cpp`, add two helper functions using the C Data Interface. Use the existing `ArrowGuard` RAII pattern (already in the file around line 34) instead of raw `malloc`/`free`:

```cpp
// Export C++ arrow::Schema → pyarrow.Schema
nb::object schema_to_pyarrow(const std::shared_ptr<arrow::Schema>& schema)
{
  ArrowGuard guard;  // RAII wrapper managing ArrowSchema C struct
  auto status = arrow::ExportSchema(*schema, guard.schema);
  if (!status.ok()) { throw std::runtime_error(status.ToString()); }

  nb::module_ pa = nb::module_::import_("pyarrow");
  nb::object Schema = pa.attr("Schema");
  return Schema.attr("_import_from_c")(reinterpret_cast<uintptr_t>(guard.schema));
}

// Import pyarrow.Schema → C++ arrow::Schema
std::shared_ptr<arrow::Schema> pyarrow_schema_to_arrow(nb::object py_schema)
{
  ArrowGuard guard;  // RAII wrapper managing ArrowSchema C struct
  py_schema.attr("_export_to_c")(reinterpret_cast<uintptr_t>(guard.schema));
  auto result = arrow::ImportSchema(guard.schema);
  if (!result.ok()) { throw std::runtime_error(result.status().ToString()); }
  return result.MoveValueUnsafe();
}
```

Note: Check the existing `ArrowGuard` implementation to confirm it has a `schema()` accessor. If it only wraps `ArrowArray`, extend it or create a separate RAII wrapper for `ArrowSchema`.

- [ ] **Step 2: Add nanobind bindings for all 18 schema structs**

For each schema struct, expose column name constants as class attributes and `schema()` as a static method. Example for FeatureSchema:

```cpp
nb::class_<FeatureSchema>(m, "FeatureSchema")
  .def_prop_ro_static("UNIQUE_ID", [](nb::handle) { return FeatureSchema::UNIQUE_ID; })
  .def_prop_ro_static("RT", [](nb::handle) { return FeatureSchema::RT; })
  // ... all column constants
  .def_static("schema", []() { return schema_to_pyarrow(FeatureSchema::schema()); });
```

Repeat for all 18 structs.

- [ ] **Step 3: Add ValidationResult binding and validate_arrow_schema function**

```cpp
nb::class_<ArrowSchemaValidation::ValidationResult>(m, "ValidationResult")
  .def_ro("valid", &ArrowSchemaValidation::ValidationResult::valid)
  .def_ro("errors", &ArrowSchemaValidation::ValidationResult::errors)
  .def("__str__", &ArrowSchemaValidation::ValidationResult::toString);

m.def("validate_arrow_schema",
  [](nb::object py_table, nb::object py_schema, const std::string& mode) {
    // Convert pyarrow table to C++ via C Data Interface
    auto table = pyarrow_to_table(py_table);
    auto schema = pyarrow_schema_to_arrow(py_schema);  // implemented in Step 1 above
    auto m = (mode == "subset") ? ArrowSchemaValidation::Mode::Subset
                                : ArrowSchemaValidation::Mode::Strict;
    return ArrowSchemaValidation::validate(table, schema, m);
  },
  nb::arg("table"), nb::arg("schema"), nb::arg("mode") = "strict");
```

- [ ] **Step 4: Build pyopenms**

Run: `cmake --build OpenMS-build --target pyopenms -j$(nproc)`

Expected: Build succeeds

- [ ] **Step 5: Write Python tests**

Add to `src/pyOpenMS/tests/unittests/test_arrow_zerocopy.py`:

```python
def test_feature_schema_round_trip():
    schema = pyopenms.FeatureSchema.schema()
    assert isinstance(schema, pa.Schema)
    assert len(schema) == 17
    assert "rt" in schema.names

def test_feature_schema_constants():
    assert pyopenms.FeatureSchema.RT == "rt"
    assert pyopenms.FeatureSchema.MZ == "mz"
    assert pyopenms.FeatureSchema.UNIQUE_ID == "unique_id"

def test_validate_arrow_schema_strict():
    schema = pyopenms.FeatureSchema.schema()
    table = pa.table({name: pa.array([0], type=field.type)
                      for name, field in zip(schema.names, schema)})
    result = pyopenms.validate_arrow_schema(table, schema)
    assert result.valid

def test_validate_arrow_schema_subset():
    schema = pyopenms.SpectraLongSchema.schema()
    # Build a subset table with just 2 columns (rt is float32 in schema, mz is float64)
    subset = pa.table({
        "rt": pa.array([1.0], type=pa.float32()),
        "mz": pa.array([500.0])
    })
    result = pyopenms.validate_arrow_schema(subset, schema, mode="subset")
    assert result.valid
```

- [ ] **Step 6: Run Python tests**

Run: `PYTHONPATH=OpenMS-build/pyOpenMS python3 -m pytest src/pyOpenMS/tests/unittests/test_arrow_zerocopy.py -v -k "schema"`

Expected: All new tests PASS

- [ ] **Step 7: Commit**

```bash
git add src/pyOpenMS/bindings/arrow_zerocopy.cpp \
  src/pyOpenMS/tests/unittests/test_arrow_zerocopy.py
git commit -m "feat: expose ArrowSchemaRegistry to Python via nanobind"
```

---

### Task 13: Update Python addons to use registry constants

**Files:**
- Modify: `src/pyOpenMS/pyopenms/addons/msexperiment.py`
- Modify: `src/pyOpenMS/pyopenms/addons/consensusmap.py`
- Modify: `src/pyOpenMS/pyopenms/addons/peptideidentificationlist.py`

- [ ] **Step 1: Update msexperiment.py**

Find the hardcoded column inventory (around line 67+) and replace string literals with `pyopenms.SpectraLongSchema.RT`, etc. The addon should import column names from the registry rather than defining its own.

- [ ] **Step 2: Update consensusmap.py**

Find the hardcoded column inventory (around line 201+) and replace with registry constants from `ConsensusFeatureExportSchema`.

- [ ] **Step 3: Update peptideidentificationlist.py**

Find the hardcoded column inventory (around line 181+) and replace with registry constants from `PSMSchema`.

- [ ] **Step 4: Run Python tests**

Run: `PYTHONPATH=OpenMS-build/pyOpenMS python3 -m pytest src/pyOpenMS/tests/unittests/ -v -k "arrow or parquet or feature_dataframe or psm_dataframe"`

Expected: All tests PASS

- [ ] **Step 5: Commit**

```bash
git add src/pyOpenMS/pyopenms/addons/msexperiment.py \
  src/pyOpenMS/pyopenms/addons/consensusmap.py \
  src/pyOpenMS/pyopenms/addons/peptideidentificationlist.py
git commit -m "refactor: use ArrowSchemaRegistry constants in Python addons"
```

---

### Task 14: Final integration test

- [ ] **Step 1: Run the full Parquet test suite**

Run: `cd OpenMS-build && ctest -R "Parquet|Arrow|QPX|OSW|Chromatogram|Mobilogram|FeatureMap|ConsensusMap|ProteinIdentification" -V`

Expected: All tests PASS

- [ ] **Step 2: Run full Python test suite**

Run: `PYTHONPATH=OpenMS-build/pyOpenMS python3 -m pytest src/pyOpenMS/tests/unittests/ -v`

Expected: All tests PASS

- [ ] **Step 3: Verify no inline schema construction remains**

Search for `arrow::schema({` in the modified source files — should only appear in ArrowSchemaRegistry.cpp and test files:

```bash
grep -rn "arrow::schema({" src/openms/source/FORMAT/ src/openms/source/ANALYSIS/OPENSWATH/
```

Expected: Only hits in `ArrowSchemaRegistry.cpp`

- [ ] **Step 4: Commit any final fixes**
