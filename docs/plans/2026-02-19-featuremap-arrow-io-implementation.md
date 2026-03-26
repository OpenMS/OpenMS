# FeatureMapArrowIO Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Full lossless round-trip export/import of FeatureMap data to Arrow/Parquet, replacing FeatureXML.

**Architecture:** Multi-file directory (`.fmd/`) containing 5 Parquet files: features, PSMs, proteins, protein_groups, search_params. Reuses existing `ProteinIdentificationArrowIO` and `QPXFile` with a `feature_id` linkage column. Features table flattens subordinates with `parent_feature_id` + `depth`.

**Tech Stack:** Apache Arrow C++ API, Parquet, OpenMS FeatureMap/Feature/ConvexHull2D, existing Arrow IO helpers.

**Design doc:** `docs/plans/2026-02-19-featuremap-arrow-io-design.md`

---

### Task 1: Create header file with class declaration

**Files:**
- Create: `src/openms/include/OpenMS/FORMAT/FeatureMapArrowIO.h`

**Step 1: Write the header file**

```cpp
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#ifdef WITH_PARQUET

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/MSExperimentArrowExport.h>

#include <memory>
#include <string>

// Forward declarations
namespace arrow
{
  class Table;
}

namespace OpenMS
{

/**
  @brief Import and export FeatureMap data to/from Apache Arrow format

  This class provides static methods to export and import FeatureMap data
  to/from Apache Arrow Tables and Parquet files. The data is stored as a
  directory of Parquet files:
  - features.parquet: Feature data (top-level + subordinates, flat with parent_feature_id)
  - psms.parquet: Peptide identifications linked to features via feature_id
  - proteins.parquet: Protein hits (via ProteinIdentificationArrowIO)
  - protein_groups.parquet: Protein groups (via ProteinIdentificationArrowIO)
  - search_params.parquet: Search parameters (via ProteinIdentificationArrowIO)

  @experimental This API is experimental and may change in future versions.

  @ingroup FileIO
*/
class OPENMS_DLLAPI FeatureMapArrowIO
{
public:
  // ==================== Export methods ====================

  /**
    @brief Export features (top-level + subordinates) to Arrow Table

    Each Feature becomes one row. Subordinate features are flattened with
    parent_feature_id and depth columns for hierarchy reconstruction.

    @param[in] feature_map The FeatureMap to export
    @return Shared pointer to Arrow Table, or nullptr on error
  */
  static std::shared_ptr<arrow::Table> exportFeaturesToArrow(
    const FeatureMap& feature_map);

  /**
    @brief Export feature-linked and unassigned peptide IDs to Arrow Table

    Extends QPXFile PSM schema with a feature_id column. Feature-linked PSMs
    have feature_id set; unassigned PSMs have feature_id = null.

    @param[in] feature_map The FeatureMap to export
    @return Shared pointer to Arrow Table, or nullptr on error
  */
  static std::shared_ptr<arrow::Table> exportPSMsToArrow(
    const FeatureMap& feature_map);

  /**
    @brief Export all FeatureMap data to a directory of Parquet files

    Creates the output directory and writes features.parquet, psms.parquet,
    proteins.parquet, protein_groups.parquet, and search_params.parquet.

    @param[in] feature_map The FeatureMap to export
    @param[in] directory Output directory path (will be created)
    @param[in] config Parquet writing options
    @return true on success, false on error
  */
  static bool exportToParquet(
    const FeatureMap& feature_map,
    const String& directory,
    const ParquetWriteConfig& config = ParquetWriteConfig{});

  // ==================== Import methods ====================

  /**
    @brief Import features from Arrow Table into FeatureMap

    Reconstructs feature hierarchy from flat rows using parent_feature_id
    and depth columns.

    @param[in] table Arrow Table with feature data
    @param[out] feature_map FeatureMap to populate
    @return true on success, false on error
  */
  static bool importFeaturesFromArrow(
    const std::shared_ptr<arrow::Table>& table,
    FeatureMap& feature_map);

  /**
    @brief Import PSMs from Arrow Table, linking to features and unassigned

    PSMs with non-null feature_id are attached to matching features.
    PSMs with null feature_id become unassigned peptide identifications.

    @param[in] table Arrow Table with PSM data
    @param[out] feature_map FeatureMap with features already populated
    @return true on success, false on error
  */
  static bool importPSMsFromArrow(
    const std::shared_ptr<arrow::Table>& table,
    FeatureMap& feature_map);

  /**
    @brief Import all data from a directory of Parquet files

    Reads features, PSMs, proteins, protein groups, and search parameters
    from the directory and reconstructs a complete FeatureMap.

    @param[in] directory Input directory path containing Parquet files
    @param[out] feature_map FeatureMap to populate
    @return true on success, false on error
  */
  static bool importFromParquet(
    const String& directory,
    FeatureMap& feature_map);
};

} // namespace OpenMS

#endif // WITH_PARQUET
```

**Step 2: Commit**

```bash
git add src/openms/include/OpenMS/FORMAT/FeatureMapArrowIO.h
git commit -m "feat: add FeatureMapArrowIO header with class declaration"
```

---

### Task 2: Create empty implementation file and register in build

**Files:**
- Create: `src/openms/source/FORMAT/FeatureMapArrowIO.cpp`
- Modify: `src/openms/source/FORMAT/sources.cmake` (line ~115, add to WITH_PARQUET block)

**Step 1: Write the stub implementation file**

```cpp
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FeatureMapArrowIO.h>

#ifdef WITH_PARQUET

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/ProteinIdentificationArrowIO.h>

#include <arrow/api.h>
#include <arrow/builder.h>
#include <arrow/io/file.h>
#include <parquet/arrow/writer.h>
#include <parquet/arrow/reader.h>
#include <parquet/properties.h>

#include <filesystem>

namespace OpenMS
{

std::shared_ptr<arrow::Table> FeatureMapArrowIO::exportFeaturesToArrow(
  const FeatureMap& /*feature_map*/)
{
  // TODO: implement
  return nullptr;
}

std::shared_ptr<arrow::Table> FeatureMapArrowIO::exportPSMsToArrow(
  const FeatureMap& /*feature_map*/)
{
  // TODO: implement
  return nullptr;
}

bool FeatureMapArrowIO::exportToParquet(
  const FeatureMap& /*feature_map*/,
  const String& /*directory*/,
  const ParquetWriteConfig& /*config*/)
{
  // TODO: implement
  return false;
}

bool FeatureMapArrowIO::importFeaturesFromArrow(
  const std::shared_ptr<arrow::Table>& /*table*/,
  FeatureMap& /*feature_map*/)
{
  // TODO: implement
  return false;
}

bool FeatureMapArrowIO::importPSMsFromArrow(
  const std::shared_ptr<arrow::Table>& /*table*/,
  FeatureMap& /*feature_map*/)
{
  // TODO: implement
  return false;
}

bool FeatureMapArrowIO::importFromParquet(
  const String& /*directory*/,
  FeatureMap& /*feature_map*/)
{
  // TODO: implement
  return false;
}

} // namespace OpenMS

#endif // WITH_PARQUET
```

**Step 2: Add to sources.cmake**

In `src/openms/source/FORMAT/sources.cmake`, inside the `if (WITH_PARQUET)` block (after the `ProteinIdentificationArrowIO.cpp` line), add:

```cmake
  list(APPEND sources_list FeatureMapArrowIO.cpp)
```

**Step 3: Build to verify it compiles**

Run: `cmake --build OpenMS-build -j$(nproc) --target OpenMS 2>&1 | tail -5`
Expected: Successful build (no errors)

**Step 4: Commit**

```bash
git add src/openms/source/FORMAT/FeatureMapArrowIO.cpp src/openms/source/FORMAT/sources.cmake
git commit -m "feat: add FeatureMapArrowIO stub implementation and register in build"
```

---

### Task 3: Create test file and register in build

**Files:**
- Create: `src/tests/class_tests/openms/source/FeatureMapArrowIO_test.cpp`
- Modify: `src/tests/class_tests/openms/executables.cmake` (line ~293, add to WITH_PARQUET block)

**Step 1: Write the initial test file with an empty FeatureMap export test**

```cpp
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/FeatureMapArrowIO.h>
///////////////////////////

#include <OpenMS/config.h>

#ifdef WITH_PARQUET

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/DATASTRUCTURES/ConvexHull2D.h>

#include <arrow/api.h>

using namespace OpenMS;
using namespace std;

START_TEST(FeatureMapArrowIO, "$Id$")

START_SECTION(exportFeaturesToArrow - empty FeatureMap)
{
  FeatureMap fm;
  auto table = FeatureMapArrowIO::exportFeaturesToArrow(fm);
  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 0)
  // Should have 17 columns per design
  TEST_EQUAL(table->num_columns(), 17)
}
END_SECTION

END_TEST

#else // WITH_PARQUET

START_TEST(FeatureMapArrowIO, "$Id$")
END_TEST

#endif // WITH_PARQUET
```

**Step 2: Register the test in executables.cmake**

In `src/tests/class_tests/openms/executables.cmake`, inside the `if(WITH_PARQUET)` block (line ~293), add `FeatureMapArrowIO_test` to the list.

**Step 3: Build the test**

Run: `cmake --build OpenMS-build -j$(nproc) --target FeatureMapArrowIO_test 2>&1 | tail -5`
Expected: Build succeeds

**Step 4: Run the test to verify it fails (since export returns nullptr)**

Run: `ctest --test-dir OpenMS-build -R FeatureMapArrowIO_test -V 2>&1 | tail -20`
Expected: FAIL (table is nullptr, empty export not yet implemented)

**Step 5: Commit**

```bash
git add src/tests/class_tests/openms/source/FeatureMapArrowIO_test.cpp src/tests/class_tests/openms/executables.cmake
git commit -m "test: add FeatureMapArrowIO test scaffold with empty export test"
```

---

### Task 4: Implement exportFeaturesToArrow

**Files:**
- Modify: `src/openms/source/FORMAT/FeatureMapArrowIO.cpp`

**Step 1: Implement the feature export**

Replace the stub `exportFeaturesToArrow` method. The implementation must:

1. Count total features (top-level + all subordinates recursively) for capacity reservation.
2. Use a recursive helper `flattenFeatures_` that walks the feature tree depth-first, writing one row per feature with `parent_feature_id` and `depth`.
3. Build Arrow columns for all 17 fields from the design doc.
4. For convex hulls: use nested `list<struct{hull_index: int32, points: list<struct{x: float64, y: float64}>}>`.
5. For bounding box: compute from `feature.getConvexHull().getBoundingBox()`.
6. For metavalues: reuse the `appendMetaValues_` pattern from `ProteinIdentificationArrowIO.cpp`.

Key patterns to follow (from `ProteinIdentificationArrowIO.cpp`):
- Reserve capacity with `builder.Reserve(num_rows)`
- Check status on every builder operation
- Use `(void)` cast for Append results
- Use the `appendMetaValues_` helper pattern for metavalue serialization
- Finalize with `builder.Finish(&array)`

The anonymous namespace helpers (`appendMetaValues_`, `writeArrowTableToParquet_`, etc.) from `ProteinIdentificationArrowIO.cpp` need to be either:
- Duplicated in the new file (simplest, follows existing pattern), OR
- Extracted to a shared header (better but separate refactoring task)

Follow the existing pattern: duplicate the helpers needed in an anonymous namespace.

**Step 2: Build**

Run: `cmake --build OpenMS-build -j$(nproc) --target FeatureMapArrowIO_test 2>&1 | tail -5`
Expected: Build succeeds

**Step 3: Run the empty export test**

Run: `ctest --test-dir OpenMS-build -R FeatureMapArrowIO_test -V 2>&1 | tail -20`
Expected: PASS (empty FeatureMap returns 0-row, 17-column table)

**Step 4: Commit**

```bash
git add src/openms/source/FORMAT/FeatureMapArrowIO.cpp
git commit -m "feat: implement FeatureMapArrowIO::exportFeaturesToArrow"
```

---

### Task 5: Add feature export tests with real data

**Files:**
- Modify: `src/tests/class_tests/openms/source/FeatureMapArrowIO_test.cpp`

**Step 1: Add test for single feature with convex hulls and metavalues**

Add a test section that:
1. Creates a FeatureMap with one feature (RT=100.0, MZ=500.5, intensity=1000, charge=2, quality=0.95).
2. Adds a convex hull with 4 points.
3. Adds metavalues (int, float, string).
4. Exports and verifies all column values, types, and nullability.

**Step 2: Add test for subordinate features**

Add a test section that:
1. Creates a feature with 2 subordinates (one of which has its own subordinate = depth 2).
2. Exports and verifies:
   - Total rows = 4 (1 top-level + 2 depth-1 + 1 depth-2)
   - `parent_feature_id` is null for top-level, correct for subordinates
   - `depth` is 0, 1, 1, 2 respectively

**Step 3: Run tests**

Run: `ctest --test-dir OpenMS-build -R FeatureMapArrowIO_test -V 2>&1 | tail -30`
Expected: All PASS

**Step 4: Commit**

```bash
git add src/tests/class_tests/openms/source/FeatureMapArrowIO_test.cpp
git commit -m "test: add feature export tests with convex hulls, metavalues, subordinates"
```

---

### Task 6: Implement importFeaturesFromArrow

**Files:**
- Modify: `src/openms/source/FORMAT/FeatureMapArrowIO.cpp`

**Step 1: Implement the feature import**

Replace the stub `importFeaturesFromArrow` method. The implementation must:

1. Use the helper functions (`getColumn_`, `getStringValue_`, `getDoubleValue_`, etc.) pattern from `ProteinIdentificationArrowIO.cpp`.
2. First pass: read all rows and create Feature objects in a `map<int64_t, Feature>` keyed by `unique_id`.
3. Second pass: for each feature with `parent_feature_id != null`, attach as subordinate to parent. Process in order of ascending depth to ensure parents exist before children.
4. Third pass: features with `parent_feature_id == null` (depth 0) are added to the FeatureMap.
5. Reconstruct convex hulls from the nested list structure.
6. Read metavalues using the `readMetaValues_` pattern.

**Step 2: Build and run tests**

Run: `cmake --build OpenMS-build -j$(nproc) --target FeatureMapArrowIO_test && ctest --test-dir OpenMS-build -R FeatureMapArrowIO_test -V 2>&1 | tail -30`
Expected: Build succeeds, existing tests pass

**Step 3: Commit**

```bash
git add src/openms/source/FORMAT/FeatureMapArrowIO.cpp
git commit -m "feat: implement FeatureMapArrowIO::importFeaturesFromArrow"
```

---

### Task 7: Add feature round-trip tests

**Files:**
- Modify: `src/tests/class_tests/openms/source/FeatureMapArrowIO_test.cpp`

**Step 1: Add round-trip test for features**

Add a test section that:
1. Creates a FeatureMap with features including subordinates, convex hulls, bounding boxes, and metavalues.
2. Exports to Arrow -> imports back.
3. Compares all field values between original and imported FeatureMap:
   - RT, MZ, intensity, charge, quality, quality_rt, quality_mz, width
   - unique_id, subordinate count, subordinate values
   - Convex hull point counts and coordinates
   - Bounding box values
   - Metavalue names, values, and types (int/float/string preservation)

**Step 2: Run tests**

Run: `cmake --build OpenMS-build -j$(nproc) --target FeatureMapArrowIO_test && ctest --test-dir OpenMS-build -R FeatureMapArrowIO_test -V 2>&1 | tail -30`
Expected: All PASS

**Step 3: Commit**

```bash
git add src/tests/class_tests/openms/source/FeatureMapArrowIO_test.cpp
git commit -m "test: add feature round-trip tests for FeatureMapArrowIO"
```

---

### Task 8: Implement exportPSMsToArrow

**Files:**
- Modify: `src/openms/source/FORMAT/FeatureMapArrowIO.cpp`

**Step 1: Implement the PSM export**

Replace the stub `exportPSMsToArrow` method. The approach:

1. Collect all PeptideIdentifications: iterate features (and their subordinates recursively) to get feature-linked PSMs, then collect unassigned PSMs.
2. For each PeptideIdentification, record the associated `feature_id` (or null for unassigned).
3. Reuse the column-building logic from `QPXFile::exportToArrow` (the PSM schema). Since `QPXFile` takes `vector<ProteinIdentification>` and `PeptideIdentificationList`, we need to adapt.
4. Add a `feature_id` (int64, nullable) column at the beginning of the schema.

**Important**: The simplest approach is to build the complete PSM table inline (following the QPXFile pattern) rather than calling QPXFile directly, since QPXFile doesn't have a feature_id column. Copy the relevant builder setup and schema from `QPXFile.cpp`, adding the feature_id column.

Alternative simpler approach: collect all PeptideIdentifications with their feature_ids, then call a shared helper. Given the complexity of PSM schema, this task may require significant code. Focus on getting the column builders right by following QPXFile.cpp exactly.

**Step 2: Build**

Run: `cmake --build OpenMS-build -j$(nproc) --target FeatureMapArrowIO_test 2>&1 | tail -5`
Expected: Build succeeds

**Step 3: Commit**

```bash
git add src/openms/source/FORMAT/FeatureMapArrowIO.cpp
git commit -m "feat: implement FeatureMapArrowIO::exportPSMsToArrow"
```

---

### Task 9: Implement importPSMsFromArrow

**Files:**
- Modify: `src/openms/source/FORMAT/FeatureMapArrowIO.cpp`

**Step 1: Implement the PSM import**

Replace the stub `importPSMsFromArrow` method. The approach:

1. Read the PSM Arrow table.
2. Group rows by `feature_id` and `P_ID` (peptide identification index) to reconstruct PeptideIdentifications.
3. For each group with non-null `feature_id`: find the matching feature in the FeatureMap (build a `unique_id -> Feature*` index first), reconstruct PeptideIdentification, and attach it.
4. For each group with null `feature_id`: reconstruct PeptideIdentification and add to `feature_map.getUnassignedPeptideIdentifications()`.
5. Reconstruct PeptideHit fields (sequence, modifications, scores, metavalues) following the QPXFile pattern in reverse.

**Step 2: Build**

Run: `cmake --build OpenMS-build -j$(nproc) --target FeatureMapArrowIO_test 2>&1 | tail -5`
Expected: Build succeeds

**Step 3: Commit**

```bash
git add src/openms/source/FORMAT/FeatureMapArrowIO.cpp
git commit -m "feat: implement FeatureMapArrowIO::importPSMsFromArrow"
```

---

### Task 10: Add PSM round-trip tests

**Files:**
- Modify: `src/tests/class_tests/openms/source/FeatureMapArrowIO_test.cpp`

**Step 1: Add PSM export and round-trip tests**

Add test sections that:
1. Create a FeatureMap with features that have PeptideIdentifications attached.
2. Add unassigned PeptideIdentifications.
3. Export PSMs to Arrow, verify `feature_id` column values.
4. Import back, verify feature-linked and unassigned PSMs survive round-trip.

**Step 2: Run tests**

Run: `cmake --build OpenMS-build -j$(nproc) --target FeatureMapArrowIO_test && ctest --test-dir OpenMS-build -R FeatureMapArrowIO_test -V 2>&1 | tail -30`
Expected: All PASS

**Step 3: Commit**

```bash
git add src/tests/class_tests/openms/source/FeatureMapArrowIO_test.cpp
git commit -m "test: add PSM round-trip tests for FeatureMapArrowIO"
```

---

### Task 11: Implement exportToParquet and importFromParquet

**Files:**
- Modify: `src/openms/source/FORMAT/FeatureMapArrowIO.cpp`

**Step 1: Implement exportToParquet**

1. Create output directory with `std::filesystem::create_directories`.
2. Call `exportFeaturesToArrow()` -> write to `<dir>/features.parquet` using `writeArrowTableToParquet_` with `file_type = "features"`.
3. Call `exportPSMsToArrow()` -> write to `<dir>/psms.parquet` with `file_type = "psms"`.
4. Delegate to `ProteinIdentificationArrowIO::exportProteinsToParquet(fm.getProteinIdentifications(), <dir>/proteins.parquet)`.
5. Same for protein_groups and search_params.

**Step 2: Implement importFromParquet**

1. Read search_params, proteins, protein_groups via `ProteinIdentificationArrowIO::importFromParquet`.
2. Set `feature_map.getProteinIdentifications()` from the result.
3. Read `features.parquet` via `readParquetTable_` -> `importFeaturesFromArrow`.
4. Read `psms.parquet` via `readParquetTable_` -> `importPSMsFromArrow`.

**Step 3: Build**

Run: `cmake --build OpenMS-build -j$(nproc) --target FeatureMapArrowIO_test 2>&1 | tail -5`
Expected: Build succeeds

**Step 4: Commit**

```bash
git add src/openms/source/FORMAT/FeatureMapArrowIO.cpp
git commit -m "feat: implement FeatureMapArrowIO exportToParquet and importFromParquet"
```

---

### Task 12: Add full Parquet directory round-trip test

**Files:**
- Modify: `src/tests/class_tests/openms/source/FeatureMapArrowIO_test.cpp`

**Step 1: Add comprehensive round-trip test**

Create a test that:
1. Builds a rich FeatureMap with:
   - 2 top-level features (one with subordinates, one without)
   - Convex hulls on features
   - PeptideIdentifications linked to features
   - Unassigned PeptideIdentifications
   - ProteinIdentifications with hits and groups
   - Metavalues on features
2. Exports to a temp directory via `exportToParquet`.
3. Imports back via `importFromParquet`.
4. Verifies all data survives: features, subordinates, hulls, PSMs, proteins, groups, metavalues.

**Step 2: Run tests**

Run: `cmake --build OpenMS-build -j$(nproc) --target FeatureMapArrowIO_test && ctest --test-dir OpenMS-build -R FeatureMapArrowIO_test -V 2>&1 | tail -30`
Expected: All PASS

**Step 3: Commit**

```bash
git add src/tests/class_tests/openms/source/FeatureMapArrowIO_test.cpp
git commit -m "test: add full Parquet directory round-trip test for FeatureMapArrowIO"
```

---

### Task 13: Add pyOpenMS nanobind binding

**Files:**
- Modify: `src/pyOpenMS/bindings/bind_format.cpp`

**Step 1: Add nanobind binding**

Add the binding to `bind_format.cpp` following the existing patterns (see `src/pyOpenMS/CLAUDE.md` for wrapping instructions):

```cpp
#include <OpenMS/FORMAT/FeatureMapArrowIO.h>

// Inside NB_MODULE:
nb::class_<OpenMS::FeatureMapArrowIO>(m, "FeatureMapArrowIO",
    "Import and export FeatureMap data to/from Apache Arrow Parquet format.")
    .def_static("exportToParquet", &OpenMS::FeatureMapArrowIO::exportToParquet,
        "feature_map"_a, "directory"_a, "config"_a,
        "Export FeatureMap to directory of Parquet files")
    .def_static("importFromParquet", &OpenMS::FeatureMapArrowIO::importFromParquet,
        "directory"_a, "feature_map"_a,
        "Import FeatureMap from directory of Parquet files")
    ;
```

Note: Only expose the high-level `exportToParquet`/`importFromParquet` to Python (not the Arrow Table methods, since Arrow C Data Interface is handled separately).

**Step 2: Build pyOpenMS to verify**

Run: `cmake --build OpenMS-build --target pyopenms -j$(nproc) 2>&1 | tail -10`

**Step 3: Commit**

```bash
git add src/pyOpenMS/bindings/bind_format.cpp
git commit -m "feat: add pyOpenMS nanobind binding for FeatureMapArrowIO"
```

---

### Task 14: Final review and cleanup

**Files:**
- Review all modified files

**Step 1: Run all Parquet-related tests**

Run: `ctest --test-dir OpenMS-build -R "Arrow|Parquet|QPX|FeatureMap" -V 2>&1 | tail -30`
Expected: All PASS

**Step 2: Verify no regressions in existing Arrow tests**

Run: `ctest --test-dir OpenMS-build -R "ProteinIdentificationArrowIO|ConsensusMapArrowExport|QPXFile" -V 2>&1 | tail -20`
Expected: All PASS

**Step 3: Final commit if any cleanup needed**

```bash
git add -A && git commit -m "chore: final cleanup for FeatureMapArrowIO"
```
