# FeatureMapArrowIO Design

## Goal

Full lossless round-trip export/import of `FeatureMap` data to Apache Arrow / Parquet, replacing FeatureXML as the storage format.

## Architecture

Multi-file directory approach. A `.fmd/` (feature map directory) contains 5 Parquet files:

```
<name>.fmd/
  features.parquet         # One row per feature (top-level + subordinates)
  psms.parquet             # QPXFile schema + feature_id linkage
  proteins.parquet         # ProteinIdentificationArrowIO (reused)
  protein_groups.parquet   # ProteinIdentificationArrowIO (reused)
  search_params.parquet    # ProteinIdentificationArrowIO (reused)
```

Reuses existing `ProteinIdentificationArrowIO` for protein/group/search_params and extends `QPXFile` PSM schema with a `feature_id` column.

## Features Table Schema

One row per feature. Top-level and subordinate features are in the same table, linked by `parent_feature_id`.

| Column | Type | Nullable | Description |
|--------|------|----------|-------------|
| `unique_id` | `int64` | no | Feature unique ID (`UniqueIdInterface`) |
| `parent_feature_id` | `int64` | yes | Null for top-level; parent's unique_id for subordinates |
| `depth` | `int32` | no | 0 = top-level, 1+ = subordinate nesting level |
| `rt` | `float64` | no | Retention time |
| `mz` | `float64` | no | m/z value |
| `intensity` | `float32` | no | Peak intensity |
| `charge` | `int32` | no | Charge state |
| `overall_quality` | `float32` | no | Overall quality measure |
| `quality_rt` | `float32` | no | Quality in RT dimension |
| `quality_mz` | `float32` | no | Quality in MZ dimension |
| `width` | `float32` | yes | FWHM (null if 0) |
| `rt_bb_min` | `float64` | yes | Bounding box RT min (from overall convex hull) |
| `rt_bb_max` | `float64` | yes | Bounding box RT max |
| `mz_bb_min` | `float64` | yes | Bounding box MZ min |
| `mz_bb_max` | `float64` | yes | Bounding box MZ max |
| `convex_hulls` | `list<struct{hull_index: int32, points: list<struct{x: float64, y: float64}>}>` | no | Per-isotope convex hull points |
| `metavalues` | `list<struct{name: utf8, value: utf8, value_type: utf8}>` | no | Arbitrary metadata |

### Design notes

- `float64` for RT/MZ coordinates (full precision for lossless round-trip).
- `float32` for intensity/quality (matches C++ `float` storage).
- Bounding box columns computed from overall convex hull on export.
- Convex hulls stored inline as nested list (avoids separate file, data is per-feature).
- Subordinate features flattened with `parent_feature_id` + `depth` columns for queryability and arbitrary nesting depth support.

## PSM Table Schema

Extends the existing QPXFile PSM schema with one column:

| Column | Type | Nullable | Description |
|--------|------|----------|-------------|
| `feature_id` | `int64` | yes | Links to `features.unique_id`. Null = unassigned peptide ID |
| *(all QPXFile PSM columns)* | | | sequence, peptidoform, modifications, scores, etc. |

### Export logic

1. For each feature, iterate `getPeptideIdentifications()` and export PSMs with `feature_id` set.
2. Then iterate `FeatureMap::getUnassignedPeptideIdentifications()` with `feature_id = null`.

### Import logic

1. Read PSM table, group by `feature_id`.
2. Non-null `feature_id` rows: reconstruct `PeptideIdentification` and attach to matching feature.
3. Null `feature_id` rows: add to `FeatureMap::unassigned_peptide_identifications_`.

## Protein-Level Data

Reuse `ProteinIdentificationArrowIO` unchanged:
- `exportProteinsToParquet` / `importProteinsFromParquet`
- `exportProteinGroupsToParquet` / `importProteinGroupsFromParquet`
- `exportSearchParamsToParquet` / `importSearchParamsFromParquet`

Joined via `run_identifier` as established.

## FeatureMap-Level Metadata

Stored as Parquet file-level key-value metadata in `features.parquet`:
- `document_id`, `file_path` (from `DocumentIdentifier`)
- `data_processing` (serialized as JSON, small and rarely queried)
- Standard QPX metadata (`qpx_version`, `creator`, `file_type`, `creation_date`, `uuid`, `software_provider`)

## Class API

```cpp
class OPENMS_DLLAPI FeatureMapArrowIO
{
public:
  // ==================== Export ====================

  /// Export features (top-level + subordinates) to Arrow Table
  static std::shared_ptr<arrow::Table> exportFeaturesToArrow(
    const FeatureMap& feature_map);

  /// Export feature-linked + unassigned PSMs to Arrow Table
  static std::shared_ptr<arrow::Table> exportPSMsToArrow(
    const FeatureMap& feature_map);

  /// Export all tables to a directory of Parquet files
  static bool exportToParquet(
    const FeatureMap& feature_map,
    const String& directory,
    const ParquetWriteConfig& config = ParquetWriteConfig{});

  // ==================== Import ====================

  /// Import features from Arrow Table
  static bool importFeaturesFromArrow(
    const std::shared_ptr<arrow::Table>& table,
    FeatureMap& feature_map);

  /// Import PSMs from Arrow Table, linking to features + unassigned
  static bool importPSMsFromArrow(
    const std::shared_ptr<arrow::Table>& table,
    FeatureMap& feature_map);

  /// Import all from a directory of Parquet files
  static bool importFromParquet(
    const String& directory,
    FeatureMap& feature_map);
};
```

### Export flow (`exportToParquet`)

1. Create output directory.
2. Call `exportFeaturesToArrow()` -> write `features.parquet`.
3. Call `exportPSMsToArrow()` -> write `psms.parquet`.
4. Delegate to `ProteinIdentificationArrowIO` for proteins, protein_groups, search_params.

### Import flow (`importFromParquet`)

1. Read `search_params.parquet` -> `ProteinIdentificationArrowIO::importSearchParamsFromParquet`.
2. Read `proteins.parquet` -> `ProteinIdentificationArrowIO::importProteinsFromParquet`.
3. Read `protein_groups.parquet` -> `ProteinIdentificationArrowIO::importProteinGroupsFromParquet`.
4. Set `feature_map.getProteinIdentifications()` from the above.
5. Read `features.parquet` -> `importFeaturesFromArrow` (reconstruct hierarchy from flat rows).
6. Read `psms.parquet` -> `importPSMsFromArrow` (link to features and unassigned).

### Feature hierarchy reconstruction (import)

1. Read all rows, build `unique_id -> Feature` map.
2. Sort by depth (ascending) to process parents before children.
3. For each row with `parent_feature_id != null`, add as subordinate to parent.
4. Rows with `parent_feature_id == null` (depth 0) become top-level features in the FeatureMap.

## Testing Strategy

- Round-trip test: create FeatureMap with subordinates, convex hulls, peptide IDs, metavalues -> export -> import -> compare.
- Edge cases: empty FeatureMap, features with no subordinates, deeply nested subordinates, features with no peptide IDs, only unassigned peptide IDs.
- Verify bounding box values match convex hull bounds.
- Verify `run_identifier` joins work across PSM and protein tables.
