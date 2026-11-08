// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FeatureMapArrowIO.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/ArrowIOHelpers.h>
#include "ArrowFileIOHelpers.h"
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/FORMAT/ProteinIdentificationArrowIO.h>
#include <OpenMS/FORMAT/ModificationDefinitionIO.h>
#include <OpenMS/FORMAT/QPXFile.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/CHEMISTRY/ProForma.h>
#include <OpenMS/METADATA/ProteinRunReferences.h>

#include <arrow/api.h>
#include <arrow/builder.h>
#include <arrow/io/file.h>
#include <parquet/arrow/writer.h>
#include <parquet/arrow/reader.h>
#include <parquet/properties.h>

#include <filesystem>
#include <unordered_set>

namespace OpenMS
{

namespace // anonymous
{
  /// Count total features recursively (top-level + all subordinates).
  Size countFeaturesRecursive_(const std::vector<Feature>& features)
  {
    Size count = 0;
    for (const auto& f : features)
    {
      ++count;
      count += countFeaturesRecursive_(f.getSubordinates());
    }
    return count;
  }

  /// Read convex hulls from a list<struct{hull_index, points: list<struct{x, y}>}> column at a given row.
  void readConvexHulls_(
    const std::shared_ptr<arrow::Array>& array,
    int64_t row,
    Feature& feature)
  {
    if (!array || array->IsNull(row)) return;
    auto list_arr = std::static_pointer_cast<arrow::ListArray>(array);
    auto hull_struct_arr = std::static_pointer_cast<arrow::StructArray>(list_arr->value_slice(row));
    if (!hull_struct_arr || hull_struct_arr->length() == 0) return;

    auto points_list_arr = std::static_pointer_cast<arrow::ListArray>(hull_struct_arr->field(1));

    for (int64_t h = 0; h < hull_struct_arr->length(); ++h)
    {
      ConvexHull2D hull;
      ConvexHull2D::PointArrayType hull_points;

      auto point_struct_arr = std::static_pointer_cast<arrow::StructArray>(points_list_arr->value_slice(h));
      if (point_struct_arr && point_struct_arr->length() > 0)
      {
        auto x_arr = std::static_pointer_cast<arrow::DoubleArray>(point_struct_arr->field(0));
        auto y_arr = std::static_pointer_cast<arrow::DoubleArray>(point_struct_arr->field(1));

        hull_points.reserve(point_struct_arr->length());
        for (int64_t p = 0; p < point_struct_arr->length(); ++p)
        {
          hull_points.emplace_back(x_arr->Value(p), y_arr->Value(p));
        }
      }

      hull.setHullPoints(hull_points);
      feature.getConvexHulls().push_back(hull);
    }
  }

} // anonymous namespace


std::shared_ptr<arrow::Table> FeatureMapArrowIO::exportFeaturesToArrow(
  const FeatureMap& feature_map)
{
  // -- Scalar column builders --
  arrow::Int64Builder unique_id_builder;
  arrow::Int64Builder parent_feature_id_builder;
  arrow::Int32Builder depth_builder;
  arrow::DoubleBuilder rt_builder, mz_builder;
  arrow::FloatBuilder intensity_builder;
  arrow::Int32Builder charge_builder;
  arrow::FloatBuilder overall_quality_builder, quality_rt_builder, quality_mz_builder;
  arrow::FloatBuilder width_builder;
  arrow::DoubleBuilder rt_bb_min_builder, rt_bb_max_builder, mz_bb_min_builder, mz_bb_max_builder;

  // -- convex_hulls: list<struct{hull_index: int32, points: list<struct{x: float64, y: float64}>}> --
  auto point_x_b = std::make_shared<arrow::DoubleBuilder>();
  auto point_y_b = std::make_shared<arrow::DoubleBuilder>();
  auto point_struct_type = arrow::struct_({
    arrow::field("x", arrow::float64()),
    arrow::field("y", arrow::float64())
  });
  auto point_struct_b = std::make_shared<arrow::StructBuilder>(
    point_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{point_x_b, point_y_b});
  auto points_list_b = std::make_shared<arrow::ListBuilder>(arrow::default_memory_pool(), point_struct_b);

  auto hull_index_b = std::make_shared<arrow::Int32Builder>();
  auto hull_struct_type = arrow::struct_({
    arrow::field("hull_index", arrow::int32()),
    arrow::field("points", arrow::list(point_struct_type))
  });
  auto hull_struct_b = std::make_shared<arrow::StructBuilder>(
    hull_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{hull_index_b, points_list_b});
  arrow::ListBuilder convex_hulls_builder(arrow::default_memory_pool(), hull_struct_b);

  // -- metavalues: list<struct{name: utf8, value: utf8, value_type: utf8}> --
  auto mv_name_b = std::make_shared<arrow::StringBuilder>();
  auto mv_value_b = std::make_shared<arrow::StringBuilder>();
  auto mv_type_b = std::make_shared<arrow::StringBuilder>();
  auto mv_struct_type = arrow::struct_({
    arrow::field("name", arrow::utf8()),
    arrow::field("value", arrow::utf8()),
    arrow::field("value_type", arrow::utf8())
  });
  auto mv_struct_b = std::make_shared<arrow::StructBuilder>(
    mv_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{mv_name_b, mv_value_b, mv_type_b});
  arrow::ListBuilder metavalues_builder(arrow::default_memory_pool(), mv_struct_b);

  // Count total features (top-level + subordinates) for capacity reservation
  Size num_rows = countFeaturesRecursive_(feature_map.getData());

  // Reserve capacity for scalar builders
  arrow::Status status;

  #define RESERVE_OR_RETURN(builder_name) \
    status = builder_name.Reserve(num_rows); \
    if (!status.ok()) { OPENMS_LOG_ERROR << "FeatureMapArrowIO: " #builder_name " Reserve failed: " << status.ToString() << std::endl; return nullptr; }

  RESERVE_OR_RETURN(unique_id_builder)
  RESERVE_OR_RETURN(parent_feature_id_builder)
  RESERVE_OR_RETURN(depth_builder)
  RESERVE_OR_RETURN(rt_builder)
  RESERVE_OR_RETURN(mz_builder)
  RESERVE_OR_RETURN(intensity_builder)
  RESERVE_OR_RETURN(charge_builder)
  RESERVE_OR_RETURN(overall_quality_builder)
  RESERVE_OR_RETURN(quality_rt_builder)
  RESERVE_OR_RETURN(quality_mz_builder)
  RESERVE_OR_RETURN(width_builder)
  RESERVE_OR_RETURN(rt_bb_min_builder)
  RESERVE_OR_RETURN(rt_bb_max_builder)
  RESERVE_OR_RETURN(mz_bb_min_builder)
  RESERVE_OR_RETURN(mz_bb_max_builder)

  #undef RESERVE_OR_RETURN

  // Exclude FWHM from metavalues (it is already stored as a dedicated column or derived from convex hulls)
  static const std::unordered_set<std::string> excluded_mvs = {"FWHM"};

  // Recursive lambda to flatten features depth-first
  std::function<void(const Feature&, int64_t, int32_t)> appendFeature =
    [&](const Feature& feature, int64_t parent_id, int32_t depth)
  {
    // === unique_id (not nullable) ===
    (void)unique_id_builder.Append(static_cast<int64_t>(feature.getUniqueId()));

    // === parent_feature_id (nullable: null for top-level) ===
    if (depth == 0)
    {
      (void)parent_feature_id_builder.AppendNull();
    }
    else
    {
      (void)parent_feature_id_builder.Append(parent_id);
    }

    // === depth (not nullable) ===
    (void)depth_builder.Append(depth);

    // === rt (not nullable) ===
    (void)rt_builder.Append(feature.getRT());

    // === mz (not nullable) ===
    (void)mz_builder.Append(feature.getMZ());

    // === intensity (not nullable) ===
    (void)intensity_builder.Append(feature.getIntensity());

    // === charge (not nullable) ===
    (void)charge_builder.Append(static_cast<int32_t>(feature.getCharge()));

    // === overall_quality (not nullable) ===
    (void)overall_quality_builder.Append(feature.getOverallQuality());

    // === quality_rt (not nullable) ===
    (void)quality_rt_builder.Append(feature.getQuality(0));

    // === quality_mz (not nullable) ===
    (void)quality_mz_builder.Append(feature.getQuality(1));

    // === width (nullable: null if 0.0, since 0.0 is the unset default; round-trip correct) ===
    float w = feature.getWidth();
    if (w == 0.0f)
    {
      (void)width_builder.AppendNull();
    }
    else
    {
      (void)width_builder.Append(w);
    }

    // === bounding box (nullable: null if no convex hulls) ===
    const auto& hulls = feature.getConvexHulls();
    if (hulls.empty())
    {
      (void)rt_bb_min_builder.AppendNull();
      (void)rt_bb_max_builder.AppendNull();
      (void)mz_bb_min_builder.AppendNull();
      (void)mz_bb_max_builder.AppendNull();
    }
    else
    {
      auto bb = feature.getConvexHull().getBoundingBox();
      (void)rt_bb_min_builder.Append(bb.minPosition()[0]);
      (void)rt_bb_max_builder.Append(bb.maxPosition()[0]);
      (void)mz_bb_min_builder.Append(bb.minPosition()[1]);
      (void)mz_bb_max_builder.Append(bb.maxPosition()[1]);
    }

    // === convex_hulls (not nullable) ===
    (void)convex_hulls_builder.Append(); // begin list for this feature
    for (int32_t hull_idx = 0; hull_idx < static_cast<int32_t>(hulls.size()); ++hull_idx)
    {
      const auto& hull = hulls[hull_idx];
      const auto& points = hull.getHullPoints();

      (void)hull_struct_b->Append(); // begin struct for this hull
      (void)hull_index_b->Append(hull_idx);
      (void)points_list_b->Append(); // begin list of points
      for (const auto& pt : points)
      {
        (void)point_struct_b->Append();
        (void)point_x_b->Append(pt[0]); // RT
        (void)point_y_b->Append(pt[1]); // MZ
      }
    }

    // === metavalues (not nullable) ===
    (void)metavalues_builder.Append();
    Internal::ArrowFileIO::appendMetaValues(feature, mv_name_b, mv_value_b, mv_type_b, mv_struct_b, excluded_mvs);

    // Recurse into subordinates
    int64_t this_id = static_cast<int64_t>(feature.getUniqueId());
    for (const auto& sub : feature.getSubordinates())
    {
      appendFeature(sub, this_id, depth + 1);
    }
  };

  // Walk all top-level features
  for (const auto& feature : feature_map)
  {
    appendFeature(feature, 0, 0);
  }

  // Finalize all arrays
  std::shared_ptr<arrow::Array> arr_unique_id, arr_parent_id, arr_depth;
  std::shared_ptr<arrow::Array> arr_rt, arr_mz, arr_intensity;
  std::shared_ptr<arrow::Array> arr_charge, arr_overall_quality;
  std::shared_ptr<arrow::Array> arr_quality_rt, arr_quality_mz;
  std::shared_ptr<arrow::Array> arr_width;
  std::shared_ptr<arrow::Array> arr_rt_bb_min, arr_rt_bb_max, arr_mz_bb_min, arr_mz_bb_max;
  std::shared_ptr<arrow::Array> arr_convex_hulls, arr_metavalues;

  #define FINISH_OR_RETURN(builder_name, arr_name) \
    status = builder_name.Finish(&arr_name); \
    if (!status.ok()) { OPENMS_LOG_ERROR << "FeatureMapArrowIO: " #builder_name " Finish failed: " << status.ToString() << std::endl; return nullptr; }

  FINISH_OR_RETURN(unique_id_builder, arr_unique_id)
  FINISH_OR_RETURN(parent_feature_id_builder, arr_parent_id)
  FINISH_OR_RETURN(depth_builder, arr_depth)
  FINISH_OR_RETURN(rt_builder, arr_rt)
  FINISH_OR_RETURN(mz_builder, arr_mz)
  FINISH_OR_RETURN(intensity_builder, arr_intensity)
  FINISH_OR_RETURN(charge_builder, arr_charge)
  FINISH_OR_RETURN(overall_quality_builder, arr_overall_quality)
  FINISH_OR_RETURN(quality_rt_builder, arr_quality_rt)
  FINISH_OR_RETURN(quality_mz_builder, arr_quality_mz)
  FINISH_OR_RETURN(width_builder, arr_width)
  FINISH_OR_RETURN(rt_bb_min_builder, arr_rt_bb_min)
  FINISH_OR_RETURN(rt_bb_max_builder, arr_rt_bb_max)
  FINISH_OR_RETURN(mz_bb_min_builder, arr_mz_bb_min)
  FINISH_OR_RETURN(mz_bb_max_builder, arr_mz_bb_max)
  FINISH_OR_RETURN(convex_hulls_builder, arr_convex_hulls)
  FINISH_OR_RETURN(metavalues_builder, arr_metavalues)

  #undef FINISH_OR_RETURN

  // Build schema from registry (17 columns)
  auto schema = FeatureSchema::schema();

  auto table = arrow::Table::Make(schema, {
    arr_unique_id, arr_parent_id, arr_depth,
    arr_rt, arr_mz, arr_intensity,
    arr_charge, arr_overall_quality,
    arr_quality_rt, arr_quality_mz,
    arr_width,
    arr_rt_bb_min, arr_rt_bb_max, arr_mz_bb_min, arr_mz_bb_max,
    arr_convex_hulls, arr_metavalues
  });

  // Validate table against registry schema
  auto validation = ArrowSchemaValidation::validate(table, FeatureSchema::schema());
  if (!validation.valid)
  {
    OPENMS_LOG_ERROR << "Schema validation failed: " << validation.toString() << "\n";
    return nullptr;
  }

  return table;
}

std::shared_ptr<arrow::Table> FeatureMapArrowIO::exportPSMsToArrow(
  const FeatureMap& feature_map)
{
  // 1. Collect all PeptideIdentifications with their associated feature_ids.
  //    For each PeptideIdentification we store (feature_unique_id, is_null).
  PeptideIdentificationList all_pep_ids;
  std::vector<std::pair<int64_t, bool>> feature_ids_per_pep_id; // (feature_id, is_null)

  // Helper: recursively collect PeptideIdentifications from a Feature and its subordinates.
  std::function<void(const Feature&)> collectFromFeature = [&](const Feature& f)
  {
    for (const auto& pep_id : f.getPeptideIdentifications())
    {
      all_pep_ids.push_back(pep_id);
      feature_ids_per_pep_id.push_back({static_cast<int64_t>(f.getUniqueId()), false});
    }
    for (const auto& sub : f.getSubordinates())
    {
      collectFromFeature(sub);
    }
  };

  for (const auto& feature : feature_map)
  {
    collectFromFeature(feature);
  }

  // Unassigned peptide identifications get feature_id = null
  for (const auto& pep_id : feature_map.getUnassignedPeptideIdentifications())
  {
    all_pep_ids.push_back(pep_id);
    feature_ids_per_pep_id.push_back({0, true}); // null
  }

  // 2. Call QPXFile::exportToArrow to produce the base PSM table.
  //    We use export_all_psms=true so all hits per PeptideIdentification are exported.
  auto base_table = QPXFile::exportToArrow(
    feature_map.getProteinIdentifications(), all_pep_ids, true);
  if (!base_table) { return nullptr; }

  // 3. Build the feature_id column.
  //    QPXFile::exportToArrow (with export_all_psms=true) produces one row per
  //    PeptideHit across all PeptideIdentifications. PeptideIdentifications
  //    with empty hits are skipped. We emit one feature_id per PeptideHit.
  arrow::Int64Builder feature_id_builder;

  for (size_t i = 0; i < all_pep_ids.size(); ++i)
  {
    // Skip PeptideIdentifications with no hits (QPXFile skips them too)
    if (all_pep_ids[i].getHits().empty())
    {
      continue;
    }

    // Emit one feature_id per PeptideHit in this PeptideIdentification
    for (size_t j = 0; j < all_pep_ids[i].getHits().size(); ++j)
    {
      if (feature_ids_per_pep_id[i].second) // is_null
      {
        (void)feature_id_builder.AppendNull();
      }
      else
      {
        (void)feature_id_builder.Append(feature_ids_per_pep_id[i].first);
      }
    }
  }

  std::shared_ptr<arrow::Array> feature_id_array;
  auto status = feature_id_builder.Finish(&feature_id_array);
  if (!status.ok())
  {
    OPENMS_LOG_ERROR << "FeatureMapArrowIO: feature_id_builder Finish failed: " << status.ToString() << std::endl;
    return nullptr;
  }

  // Verify row count matches
  if (feature_id_array->length() != base_table->num_rows())
  {
    OPENMS_LOG_ERROR << "FeatureMapArrowIO: feature_id column length (" << feature_id_array->length()
                     << ") does not match table row count (" << base_table->num_rows() << ")" << std::endl;
    return nullptr;
  }

  // 4. Add feature_id as the first column in the table.
  auto chunked_feature_id = std::make_shared<arrow::ChunkedArray>(feature_id_array);
  auto result = base_table->AddColumn(0, arrow::field("feature_unique_id", arrow::int64()), chunked_feature_id);
  if (!result.ok())
  {
    OPENMS_LOG_ERROR << "FeatureMapArrowIO: AddColumn failed: " << result.status().ToString() << std::endl;
    return nullptr;
  }
  return *result;
}

bool FeatureMapArrowIO::exportToParquet(
  const FeatureMap& feature_map,
  const std::string& directory,
  const ParquetWriteConfig& config)
{
  // Mirror XMLHandler::checkUniqueIdentifiers_ — fail before any file is opened
  // so we never leave a partial .featureparquet behind. Throws Exception::InvalidValue.
  ProteinIdentificationArrowIO::checkUniqueIdentifiers(feature_map.getProteinIdentifications());
  // every peptide identification needs its protein identification run
  ProteinRunReferences::check(feature_map);

  // 1. Create output directory
  try
  {
    std::filesystem::create_directories(std::string(directory));
  }
  catch (const std::filesystem::filesystem_error& e)
  {
    OPENMS_LOG_ERROR << "FeatureMapArrowIO: Failed to create directory: " << e.what() << std::endl;
    return false;
  }

  // 2. Export features table (with FeatureMap-level metadata)
  auto features_table = exportFeaturesToArrow(feature_map);
  if (!features_table)
  {
    OPENMS_LOG_ERROR << "FeatureMapArrowIO: Failed to create features Arrow table" << std::endl;
    return false;
  }

  // Collect FeatureMap-level metadata (DocumentIdentifier + DataProcessing + MetaValues).
  // The fmap_metavalues entry carries every FeatureMap-level MetaValue (notably
  // `spectra_data`, set by FeatureMap::setPrimaryMSRunPath at FeatureMap.cpp:415).
  // Mirrors the cmap_metavalues key in ConsensusMapArrowIO.cpp.
  std::unordered_map<std::string, std::string> feature_map_metadata;
  feature_map_metadata["document_id"] = feature_map.getIdentifier();
  feature_map_metadata["loaded_file_path"] = feature_map.getLoadedFilePath();
  feature_map_metadata["loaded_file_type"] = FileTypes::typeToName(feature_map.getLoadedFileType());
  feature_map_metadata["data_processing"] = Internal::ArrowFileIO::serializeDataProcessing(feature_map.getDataProcessing());
  feature_map_metadata["fmap_metavalues"] = Internal::ArrowFileIO::serializeMetaValues(feature_map);
  // The map's own UniqueIdInterface value -- distinct from the per-feature unique ids in the
  // table and from DocumentIdentifier above. Consumers key on it: ProteomicsLFQ copies it into
  // the consensus column header of the run the map came from.
  feature_map_metadata["map_unique_id"] = std::to_string(feature_map.getUniqueId());

  if (!Internal::ArrowFileIO::writeParquetTable(features_table, directory + "/features.parquet", "features", config, "FeatureMapArrowIO", feature_map_metadata))
  {
    return false;
  }

  // 3. Export PSMs table
  auto psms_table = exportPSMsToArrow(feature_map);
  if (!psms_table)
  {
    OPENMS_LOG_ERROR << "FeatureMapArrowIO: Failed to create PSMs Arrow table" << std::endl;
    return false;
  }
  if (!Internal::ArrowFileIO::writeParquetTable(psms_table, directory + "/psms.parquet", "psms", config, "FeatureMapArrowIO"))
  {
    return false;
  }

  // 4. Delegate protein data to ProteinIdentificationArrowIO
  const auto& prot_ids = feature_map.getProteinIdentifications();
  if (!ProteinIdentificationArrowIO::exportProteinsToParquet(
          prot_ids, directory + "/proteins.parquet", config))
  {
    return false;
  }
  if (!ProteinIdentificationArrowIO::exportProteinGroupsToParquet(
          prot_ids, directory + "/protein_groups.parquet", config))
  {
    return false;
  }
  if (!ProteinIdentificationArrowIO::exportSearchParamsToParquet(
          prot_ids, directory + "/search_params.parquet", config,
          ModificationDefinitionIO::encodeByRun(prot_ids, ModificationDefinitionIO::collect(feature_map))))
  {
    return false;
  }

  return true;
}

bool FeatureMapArrowIO::importFeaturesFromArrow(
  const std::shared_ptr<arrow::Table>& table,
  FeatureMap& feature_map)
{
  if (!table)
  {
    OPENMS_LOG_ERROR << "FeatureMapArrowIO: null table passed to importFeaturesFromArrow" << std::endl;
    return false;
  }

  // Combine chunks for reliable single-chunk access
  auto combined_result = table->CombineChunks(arrow::default_memory_pool());
  if (!combined_result.ok())
  {
    OPENMS_LOG_ERROR << "FeatureMapArrowIO: Failed to combine chunks" << std::endl;
    return false;
  }
  const auto& tbl = *combined_result;

  int64_t num_rows = tbl->num_rows();
  if (num_rows == 0)
  {
    return true;
  }

  // Validate table schema against registry (subset mode — file may have extra columns)
  auto validation = ArrowSchemaValidation::validate(tbl, FeatureSchema::schema(), ArrowSchemaValidation::Mode::Subset);
  if (!validation.valid)
  {
    OPENMS_LOG_ERROR << "Incompatible schema: " << validation.toString() << "\n";
    return false;
  }

  // Get all columns
  auto col_unique_id = ArrowIOHelpers::getColumn(tbl, FeatureSchema::UNIQUE_ID, true, "FeatureMapArrowIO");
  auto col_parent_id = ArrowIOHelpers::getColumn(tbl, FeatureSchema::PARENT_FEATURE_ID, true, "FeatureMapArrowIO");
  auto col_depth = ArrowIOHelpers::getColumn(tbl, FeatureSchema::DEPTH, true, "FeatureMapArrowIO");
  auto col_rt = ArrowIOHelpers::getColumn(tbl, FeatureSchema::RT, true, "FeatureMapArrowIO");
  auto col_mz = ArrowIOHelpers::getColumn(tbl, FeatureSchema::MZ, true, "FeatureMapArrowIO");
  auto col_intensity = ArrowIOHelpers::getColumn(tbl, FeatureSchema::INTENSITY, true, "FeatureMapArrowIO");
  auto col_charge = ArrowIOHelpers::getColumn(tbl, FeatureSchema::CHARGE, true, "FeatureMapArrowIO");
  auto col_overall_quality = ArrowIOHelpers::getColumn(tbl, FeatureSchema::QUALITY, true, "FeatureMapArrowIO");
  auto col_quality_rt = ArrowIOHelpers::getColumn(tbl, FeatureSchema::QUALITY_RT, true, "FeatureMapArrowIO");
  auto col_quality_mz = ArrowIOHelpers::getColumn(tbl, FeatureSchema::QUALITY_MZ, true, "FeatureMapArrowIO");
  auto col_width = ArrowIOHelpers::getColumn(tbl, FeatureSchema::WIDTH, /*required=*/false, "FeatureMapArrowIO");
  auto col_convex_hulls = ArrowIOHelpers::getColumn(tbl, FeatureSchema::CONVEX_HULLS, /*required=*/false, "FeatureMapArrowIO");
  auto col_metavalues = ArrowIOHelpers::getColumn(tbl, FeatureSchema::METAVALUES, /*required=*/false, "FeatureMapArrowIO");

  if (!col_unique_id || !col_parent_id || !col_depth || !col_rt || !col_mz ||
      !col_intensity || !col_charge || !col_overall_quality ||
      !col_quality_rt || !col_quality_mz)
  {
    OPENMS_LOG_ERROR << "FeatureMapArrowIO: Missing one or more required columns" << std::endl;
    return false;
  }

  // Build a vector of {depth, row_index, unique_id, parent_id, Feature} entries
  struct FeatureEntry
  {
    int32_t depth;
    int64_t row_index;
    int64_t unique_id;
    int64_t parent_id;  // -1 if top-level
    Feature feature;
  };

  std::vector<FeatureEntry> entries;
  entries.reserve(num_rows);

  // No metavalue keys to exclude for features
  static const std::unordered_set<std::string> excluded_mvs = {};

  for (int64_t i = 0; i < num_rows; ++i)
  {
    FeatureEntry entry;
    entry.row_index = i;
    entry.depth = ArrowIOHelpers::getInt32Value(col_depth, i, 0);
    entry.unique_id = ArrowIOHelpers::getInt64Value(col_unique_id, i, 0);
    entry.parent_id = ArrowIOHelpers::isNull(col_parent_id, i) ? -1 : ArrowIOHelpers::getInt64Value(col_parent_id, i, 0);

    Feature& f = entry.feature;
    f.setUniqueId(static_cast<UInt64>(entry.unique_id));
    f.setRT(ArrowIOHelpers::getDoubleValue(col_rt, i));
    f.setMZ(ArrowIOHelpers::getDoubleValue(col_mz, i));
    f.setIntensity(ArrowIOHelpers::getFloatValue(col_intensity, i));
    f.setCharge(static_cast<int>(ArrowIOHelpers::getInt32Value(col_charge, i)));
    f.setOverallQuality(ArrowIOHelpers::getFloatValue(col_overall_quality, i));
    f.setQuality(0, ArrowIOHelpers::getFloatValue(col_quality_rt, i));
    f.setQuality(1, ArrowIOHelpers::getFloatValue(col_quality_mz, i));

    // width: null means unset (default 0.0); only set if non-null and non-zero
    float w = ArrowIOHelpers::getFloatValue(col_width, i, 0.0f);
    if (w != 0.0f)
    {
      f.setWidth(w);
    }

    // Read convex hulls
    if (col_convex_hulls)
    {
      readConvexHulls_(col_convex_hulls, i, f);
    }

    // Read metavalues
    if (col_metavalues)
    {
      ArrowIOHelpers::readMetaValues(col_metavalues, i, f, excluded_mvs);
    }

    entries.push_back(std::move(entry));
  }

  // Sort by depth descending so children are processed before parents.
  // This allows us to build the tree bottom-up: attach children to parents
  // in the entries vector first, then add fully assembled top-level features
  // to the FeatureMap.
  std::stable_sort(entries.begin(), entries.end(),
    [](const FeatureEntry& a, const FeatureEntry& b) { return a.depth > b.depth; });

  // Build unique_id -> index-in-entries lookup map
  std::unordered_map<int64_t, size_t> id_to_index;
  for (size_t i = 0; i < entries.size(); ++i)
  {
    id_to_index[entries[i].unique_id] = i;
  }

  // Bottom-up assembly: for each entry (deepest first), attach it to its parent
  for (auto& entry : entries)
  {
    if (entry.parent_id != -1)
    {
      auto it = id_to_index.find(entry.parent_id);
      if (it != id_to_index.end())
      {
        entries[it->second].feature.getSubordinates().push_back(std::move(entry.feature));
      }
      else
      {
        OPENMS_LOG_WARN << "FeatureMapArrowIO: Could not find parent feature with id "
                        << entry.parent_id << " for feature " << entry.unique_id
                        << ". Will be added as top-level." << std::endl;
        // Mark as top-level by clearing parent_id
        entry.parent_id = -1;
      }
    }
  }

  // Add top-level features to the FeatureMap (preserving original order by row_index)
  // First, collect top-level entries sorted by row_index
  std::vector<size_t> top_level_indices;
  for (size_t i = 0; i < entries.size(); ++i)
  {
    if (entries[i].parent_id == -1)
    {
      top_level_indices.push_back(i);
    }
  }
  std::sort(top_level_indices.begin(), top_level_indices.end(),
    [&](size_t a, size_t b) { return entries[a].row_index < entries[b].row_index; });

  for (size_t idx : top_level_indices)
  {
    feature_map.push_back(std::move(entries[idx].feature));
  }

  return true;
}

bool FeatureMapArrowIO::importPSMsFromArrow(
  const std::shared_ptr<arrow::Table>& table,
  FeatureMap& feature_map)
{
  if (!table || table->num_rows() == 0) { return true; }

  // Combine chunks for reliable single-chunk access
  auto combined_result = table->CombineChunks(arrow::default_memory_pool());
  if (!combined_result.ok())
  {
    OPENMS_LOG_ERROR << "FeatureMapArrowIO: Failed to combine chunks in importPSMsFromArrow" << std::endl;
    return false;
  }
  const auto& tbl = *combined_result;
  int64_t num_rows = tbl->num_rows();

  auto psm_validation = ArrowSchemaValidation::validate(tbl, PSMSchema::schema(), ArrowSchemaValidation::Mode::Subset);
  if (!psm_validation.valid)
  {
    OPENMS_LOG_ERROR << "Incompatible PSM schema: " << psm_validation.toString() << "\n";
    return false;
  }

  // Build feature lookup: unique_id -> Feature* (recursively includes subordinates)
  std::unordered_map<int64_t, Feature*> feature_lookup;
  std::function<void(Feature&)> buildLookup = [&](Feature& f)
  {
    feature_lookup[static_cast<int64_t>(f.getUniqueId())] = &f;
    for (auto& sub : f.getSubordinates())
    {
      buildLookup(sub);
    }
  };
  for (auto& f : feature_map)
  {
    buildLookup(f);
  }

  // Read columns
  auto col_feature_id = ArrowIOHelpers::getColumn(tbl, "feature_unique_id", true, "FeatureMapArrowIO");
  auto col_p_id = ArrowIOHelpers::getColumn(tbl, PSMSchema::PEPTIDE_IDENTIFICATION_INDEX, true, "FeatureMapArrowIO");
  auto col_peptidoform = ArrowIOHelpers::getColumn(tbl, PSMSchema::PEPTIDOFORM, /*required=*/false, "FeatureMapArrowIO");
  auto col_sequence = ArrowIOHelpers::getColumn(tbl, PSMSchema::SEQUENCE, /*required=*/false, "FeatureMapArrowIO");
  auto col_charge = ArrowIOHelpers::getColumn(tbl, PSMSchema::PRECURSOR_CHARGE, true, "FeatureMapArrowIO");
  auto col_score = ArrowIOHelpers::getColumn(tbl, PSMSchema::SCORE, true, "FeatureMapArrowIO");
  auto col_score_type = ArrowIOHelpers::getColumn(tbl, PSMSchema::SCORE_TYPE, true, "FeatureMapArrowIO");
  // hit_index is positional (analytics view); rank semantics travel via psm_metavalues.
  auto col_rt = ArrowIOHelpers::getColumn(tbl, PSMSchema::RT, /*required=*/false, "FeatureMapArrowIO");
  auto col_mz = ArrowIOHelpers::getColumn(tbl, PSMSchema::OBSERVED_MZ, /*required=*/false, "FeatureMapArrowIO");
  auto col_spec_ref = ArrowIOHelpers::getColumn(tbl, PSMSchema::SPECTRUM_REFERENCE, /*required=*/false, "FeatureMapArrowIO");
  auto col_run_id = ArrowIOHelpers::getColumn(tbl, PSMSchema::RUN_IDENTIFIER, /*required=*/false, "FeatureMapArrowIO");
  auto col_is_decoy = ArrowIOHelpers::getColumn(tbl, PSMSchema::IS_DECOY, /*required=*/false, "FeatureMapArrowIO");
  auto col_protein_accs = ArrowIOHelpers::getColumn(tbl, PSMSchema::PROTEIN_ACCESSIONS, /*required=*/false, "FeatureMapArrowIO");
  auto col_additional_scores = ArrowIOHelpers::getColumn(tbl, PSMSchema::ADDITIONAL_SCORES, /*required=*/false, "FeatureMapArrowIO");
  auto col_psm_metavalues = ArrowIOHelpers::getColumn(tbl, PSMSchema::PSM_METAVALUES, /*required=*/false, "FeatureMapArrowIO");
  auto col_spectrum_metavalues = ArrowIOHelpers::getColumn(tbl, PSMSchema::SPECTRUM_METAVALUES, /*required=*/false, "FeatureMapArrowIO");
  auto col_predicted_rt = ArrowIOHelpers::getColumn(tbl, PSMSchema::PREDICTED_RT, /*required=*/false, "FeatureMapArrowIO");
  auto col_ion_mobility = ArrowIOHelpers::getColumn(tbl, PSMSchema::ION_MOBILITY, /*required=*/false, "FeatureMapArrowIO");
  auto col_hsb = ArrowIOHelpers::getColumn(tbl, PSMSchema::HIGHER_SCORE_BETTER, /*required=*/false, "FeatureMapArrowIO");
  auto col_scan = ArrowIOHelpers::getColumn(tbl, PSMSchema::SCAN, /*required=*/false, "FeatureMapArrowIO");
  auto col_ref_file = ArrowIOHelpers::getColumn(tbl, PSMSchema::REFERENCE_FILE_NAME, /*required=*/false, "FeatureMapArrowIO");

  if (!col_feature_id || !col_p_id || !col_charge || !col_score || !col_score_type)
  {
    OPENMS_LOG_ERROR << "FeatureMapArrowIO: Missing required columns for PSM import" << std::endl;
    return false;
  }

  // Group rows by P_ID to reconstruct PeptideIdentifications.
  // Each unique P_ID value corresponds to one PeptideIdentification.
  // Rows with the same P_ID share the same PeptideIdentification-level data
  // (rt, mz, score_type, spectrum_reference, run_identifier).
  //
  // We also track the feature_id for each group so we know where to attach it.
  // All rows in the same P_ID group should have the same feature_id.

  // Build a lookup for higher_score_better from ProteinIdentifications
  std::unordered_map<std::string, bool> higher_score_better_lookup;
  for (const auto& prot_id : feature_map.getProteinIdentifications())
  {
    higher_score_better_lookup[prot_id.getIdentifier()] = prot_id.isHigherScoreBetter();
  }

  // Process rows in order, grouping consecutive rows with the same P_ID.
  struct PepIdGroup
  {
    int64_t feature_id;
    bool feature_id_is_null;
    PeptideIdentification pep_id;
  };

  std::vector<PepIdGroup> groups;
  std::unordered_map<int32_t, size_t> p_id_to_idx;

  for (int64_t row = 0; row < num_rows; ++row)
  {
    int32_t p_id = ArrowIOHelpers::getInt32Value(col_p_id, row, -1);

    // Group PSM rows by P_ID using a first-seen lookup rather than assuming the
    // rows are contiguous. A standard OpenMS export writes rows of one P_ID
    // contiguously, but an externally sorted/filtered Parquet file may interleave
    // them; grouping by value (as QPXFile and ConsensusMapArrowIO already do) keeps
    // all hits of one logical PeptideIdentification together in either case.
    auto [pid_it, pid_inserted] = p_id_to_idx.try_emplace(p_id, groups.size());
    if (pid_inserted)
    {
      PepIdGroup group;
      group.feature_id_is_null = ArrowIOHelpers::isNull(col_feature_id, row);
      group.feature_id = group.feature_id_is_null ? 0 : ArrowIOHelpers::getInt64Value(col_feature_id, row, 0);

      // Set PeptideIdentification-level fields
      PeptideIdentification& pid = group.pep_id;
      pid.setScoreType(ArrowIOHelpers::getStringValue(col_score_type, row));
      // Run identifier (links to ProteinIdentification)
      if (col_run_id && !ArrowIOHelpers::isNull(col_run_id, row))
      {
        pid.setIdentifier(ArrowIOHelpers::getStringValue(col_run_id, row));
      }

      // higher_score_better: prefer per-PSM column, fall back to ProteinIdentification lookup
      if (col_hsb && !ArrowIOHelpers::isNull(col_hsb, row))
      {
        pid.setHigherScoreBetter(ArrowIOHelpers::getBoolValue(col_hsb, row, true));
      }
      else if (col_run_id && !ArrowIOHelpers::isNull(col_run_id, row))
      {
        auto hsb_it = higher_score_better_lookup.find(pid.getIdentifier());
        if (hsb_it != higher_score_better_lookup.end())
        {
          pid.setHigherScoreBetter(hsb_it->second);
        }
        else
        {
          pid.setHigherScoreBetter(true);
        }
      }
      else
      {
        pid.setHigherScoreBetter(true);
      }

      // RT
      if (col_rt && !ArrowIOHelpers::isNull(col_rt, row))
      {
        pid.setRT(ArrowIOHelpers::getDoubleValue(col_rt, row));
      }

      // MZ
      if (col_mz && !ArrowIOHelpers::isNull(col_mz, row))
      {
        pid.setMZ(ArrowIOHelpers::getDoubleValue(col_mz, row));
      }

      // Spectrum reference
      if (col_spec_ref && !ArrowIOHelpers::isNull(col_spec_ref, row))
      {
        pid.setSpectrumReference(ArrowIOHelpers::getStringValue(col_spec_ref, row));
      }

      // spectrum_metavalues -> PeptideIdentification metavalues
      if (col_spectrum_metavalues)
      {
        ArrowIOHelpers::readMetaValues(col_spectrum_metavalues, row, pid);
      }

      groups.push_back(std::move(group));
    }
    const size_t gidx = pid_it->second;

    // Add a PeptideHit to the current group
    PeptideHit hit;

    // Reconstruct sequence from peptidoform (ProForma) if available, else from sequence column
    if (col_peptidoform && !ArrowIOHelpers::isNull(col_peptidoform, row))
    {
      std::string peptidoform_str = ArrowIOHelpers::getStringValue(col_peptidoform, row);
      if (!peptidoform_str.empty())
      {
        try
        {
          auto pf = ProForma::parse(peptidoform_str);
          hit.setSequence(ProForma::toAASequence(pf, ProForma::ConversionPolicy::BEST_EFFORT));
        }
        catch (...)
        {
          // Fallback: try unmodified sequence
          if (col_sequence && !ArrowIOHelpers::isNull(col_sequence, row))
          {
            hit.setSequence(AASequence::fromString(ArrowIOHelpers::getStringValue(col_sequence, row)));
          }
        }
      }
    }
    else if (col_sequence && !ArrowIOHelpers::isNull(col_sequence, row))
    {
      hit.setSequence(AASequence::fromString(ArrowIOHelpers::getStringValue(col_sequence, row)));
    }

    hit.setCharge(static_cast<Int>(ArrowIOHelpers::getInt32Value(col_charge, row, 0)));
    hit.setScore(ArrowIOHelpers::getDoubleValue(col_score, row, 0.0));

    // is_decoy -> target_decoy metavalue
    if (col_is_decoy && !ArrowIOHelpers::isNull(col_is_decoy, row))
    {
      bool is_decoy = ArrowIOHelpers::getBoolValue(col_is_decoy, row, false);
      hit.setMetaValue("target_decoy", is_decoy ? "decoy" : "target");
    }

    // protein_accessions -> PeptideEvidence
    if (col_protein_accs && !ArrowIOHelpers::isNull(col_protein_accs, row))
    {
      auto list_arr = std::static_pointer_cast<arrow::ListArray>(col_protein_accs);
      auto struct_arr = std::static_pointer_cast<arrow::StructArray>(list_arr->values());
      auto acc_arr    = std::static_pointer_cast<arrow::StringArray>(struct_arr->GetFieldByName("accession"));
      auto before_arr = std::static_pointer_cast<arrow::StringArray>(struct_arr->GetFieldByName("aa_before"));
      auto after_arr  = std::static_pointer_cast<arrow::StringArray>(struct_arr->GetFieldByName("aa_after"));
      auto start_arr  = std::static_pointer_cast<arrow::Int32Array>(struct_arr->GetFieldByName("start"));
      auto end_arr    = std::static_pointer_cast<arrow::Int32Array>(struct_arr->GetFieldByName("end"));
      int64_t lstart = list_arr->value_offset(row);
      int64_t lend = lstart + list_arr->value_length(row);
      for (int64_t k = lstart; k < lend; ++k)
      {
        PeptideEvidence ev;
        ev.setProteinAccession(acc_arr->GetString(k));
        const std::string before_s = before_arr->IsNull(k) ? std::string{} : before_arr->GetString(k);
        ev.setAABefore(before_s.empty() ? PeptideEvidence::UNKNOWN_AA : before_s[0]);
        const std::string after_s = after_arr->IsNull(k) ? std::string{} : after_arr->GetString(k);
        ev.setAAAfter(after_s.empty() ? PeptideEvidence::UNKNOWN_AA : after_s[0]);
        ev.setStart(start_arr->IsNull(k) ? PeptideEvidence::UNKNOWN_POSITION : start_arr->Value(k));
        ev.setEnd  (end_arr  ->IsNull(k) ? PeptideEvidence::UNKNOWN_POSITION : end_arr  ->Value(k));
        hit.addPeptideEvidence(ev);
      }
    }

    // additional_scores -> metavalues on PeptideHit
    if (col_additional_scores && !ArrowIOHelpers::isNull(col_additional_scores, row))
    {
      auto list_arr = std::static_pointer_cast<arrow::ListArray>(col_additional_scores);
      auto struct_arr = std::static_pointer_cast<arrow::StructArray>(list_arr->values());
      auto names_arr = std::static_pointer_cast<arrow::StringArray>(struct_arr->GetFieldByName("score_name"));
      auto values_arr = std::static_pointer_cast<arrow::DoubleArray>(struct_arr->GetFieldByName("score_value"));

      int64_t start = list_arr->value_offset(row);
      int64_t end = start + list_arr->value_length(row);
      for (int64_t k = start; k < end; ++k)
      {
        std::string name = names_arr->GetString(k);
        double value = values_arr->Value(k);
        hit.setMetaValue(name, value);
      }
    }

    // predicted_rt -> PeptideHit metavalue
    if (col_predicted_rt && !ArrowIOHelpers::isNull(col_predicted_rt, row))
    {
      hit.setMetaValue("predicted_RT", ArrowIOHelpers::getDoubleValue(col_predicted_rt, row));
    }

    // ion_mobility -> PeptideHit metavalue
    if (col_ion_mobility && !ArrowIOHelpers::isNull(col_ion_mobility, row))
    {
      hit.setMetaValue("ion_mobility", ArrowIOHelpers::getDoubleValue(col_ion_mobility, row));
    }

    // scan -> PeptideHit metavalue
    if (col_scan && !ArrowIOHelpers::isNull(col_scan, row))
    {
      hit.setMetaValue("scan", static_cast<int>(ArrowIOHelpers::getInt32Value(col_scan, row)));
    }

    // reference_file_name -> PeptideHit metavalue
    if (col_ref_file && !ArrowIOHelpers::isNull(col_ref_file, row))
    {
      hit.setMetaValue("reference_file_name", ArrowIOHelpers::getStringValue(col_ref_file, row));
    }

    // psm_metavalues -> PeptideHit metavalues (exclude keys that have dedicated columns)
    if (col_psm_metavalues)
    {
      // 'target_decoy' is read from is_decoy above, and is only here if that cannot express it (e.g. "target+decoy")
      static const std::unordered_set<std::string> psm_excluded_mvs =
        {"predicted_RT", "predicted_rt", "ion_mobility", "IM",
         "scan", "reference_file_name"};
      ArrowIOHelpers::readMetaValues(col_psm_metavalues, row, hit, psm_excluded_mvs);
    }

    groups[gidx].pep_id.getHits().push_back(std::move(hit));
  }

  // Assign PeptideIdentifications to features or as unassigned
  for (auto& group : groups)
  {
    if (group.feature_id_is_null)
    {
      // Unassigned
      feature_map.getUnassignedPeptideIdentifications().push_back(std::move(group.pep_id));
    }
    else
    {
      auto it = feature_lookup.find(group.feature_id);
      if (it != feature_lookup.end())
      {
        it->second->getPeptideIdentifications().push_back(std::move(group.pep_id));
      }
      else
      {
        // Refuse rather than re-file the identification as unassigned. The two identification
        // lists of a FeatureMap partition the PSMs, and which list a PSM is in is part of the
        // result: FDR is estimated over both together, and every exporter reports them
        // separately. Silently moving one across that boundary yields a map that parses
        // cleanly and describes a different experiment from the one that was written, which
        // is the failure mode that is hardest to notice afterwards.
        OPENMS_LOG_ERROR << "FeatureMapArrowIO: PSM references feature unique id "
                         << group.feature_id << ", which is not in features.parquet. The file is "
                         << "inconsistent." << std::endl;
        return false;
      }
    }
  }

  return true;
}

bool FeatureMapArrowIO::importFromParquet(
  const std::string& directory,
  FeatureMap& feature_map)
{
  feature_map = FeatureMap{};

  // 1. Import protein identification data
  std::vector<ProteinIdentification> prot_ids;
  if (!ProteinIdentificationArrowIO::importFromParquet(
          directory + "/proteins.parquet",
          directory + "/protein_groups.parquet",
          directory + "/search_params.parquet",
          prot_ids))
  {
    return false;
  }
  feature_map.setProteinIdentifications(prot_ids);

  // 2. Import features and FeatureMap-level metadata
  auto features_table = Internal::ArrowFileIO::readParquetTable(directory + "/features.parquet", "FeatureMapArrowIO");
  if (!features_table) { return false; }

  // Extract FeatureMap-level metadata from file-level key-value metadata
  auto schema_md = features_table->schema()->metadata();
  if (schema_md)
  {
    // DocumentIdentifier fields
    int idx = schema_md->FindKey("document_id");
    if (idx >= 0) feature_map.setIdentifier(schema_md->value(idx));

    idx = schema_md->FindKey("loaded_file_path");
    if (idx >= 0) feature_map.setLoadedFilePath(schema_md->value(idx));

    // loaded_file_type: stored as type name string (e.g. "featureXML").
    // Note: DocumentIdentifier::setLoadedFileType(String) expects a file path and
    // determines type by content/extension, so we cannot directly restore the type
    // from a type name. The type can be re-derived from loaded_file_path if needed.

    // DataProcessing
    idx = schema_md->FindKey("data_processing");
    if (idx >= 0)
    {
      feature_map.setDataProcessing(Internal::ArrowFileIO::deserializeDataProcessing(schema_md->value(idx), "FeatureMapArrowIO"));
    }

    // FeatureMap-level MetaValues (e.g. `spectra_data` populated by setPrimaryMSRunPath).
    // Missing key (pre-fix .featureparquet) → no meta-values restored, matching legacy
    // behavior. Mirrors the cmap_metavalues handling in ConsensusMapArrowIO::importFromParquet.
    idx = schema_md->FindKey("fmap_metavalues");
    if (idx >= 0)
    {
      Internal::ArrowFileIO::deserializeMetaValues(schema_md->value(idx), feature_map);
    }

    // The map's own unique id. Absent in files written before it was stored, which then keep
    // the invalid id the map starts with, exactly as they did before.
    idx = schema_md->FindKey("map_unique_id");
    if (idx >= 0)
    {
      try
      {
        feature_map.setUniqueId(static_cast<UInt64>(std::stoull(schema_md->value(idx))));
      }
      catch (const std::exception&)
      {
        OPENMS_LOG_WARN << "FeatureMapArrowIO: could not read map_unique_id '"
                        << schema_md->value(idx) << "'; leaving the map id unset." << std::endl;
      }
    }
  }

  if (!importFeaturesFromArrow(features_table, feature_map))
  {
    return false;
  }

  // 3. Import PSMs (links to features already in the map)
  auto psms_table = Internal::ArrowFileIO::readParquetTable(directory + "/psms.parquet", "FeatureMapArrowIO");
  if (!psms_table) { return false; }
  if (!importPSMsFromArrow(psms_table, feature_map))
  {
    return false;
  }
  // every peptide identification needs its protein identification run (checked on the stored identifiers)
  ProteinRunReferences::check(feature_map);

  // 4. Synthesize fresh ProtID identifiers + apply rename to every pep_id collection
  //    we own (per-feature + unassigned). Mirrors IdXMLFile.cpp:530 — the stored
  //    identifier becomes informational; the in-memory identifier downstream sees
  //    is freshly synthesized with a UniqueIdGenerator suffix.
  {
    auto& prot_ids = feature_map.getProteinIdentifications();
    auto rename = ProteinIdentificationArrowIO::synthesizeRunIdentifiers(prot_ids);

    // features and their subordinates (which can carry peptide identifications, too)
    std::function<void(Feature&)> renameFeature = [&](Feature& f)
    {
      ProteinIdentificationArrowIO::applyRunIdentifierRename(rename, f.getPeptideIdentifications());
      for (auto& sub : f.getSubordinates())
      {
        renameFeature(sub);
      }
    };
    for (auto& feature : feature_map)
    {
      renameFeature(feature);
    }
    ProteinIdentificationArrowIO::applyRunIdentifierRename(
        rename, feature_map.getUnassignedPeptideIdentifications());
  }

  return true;
}

} // namespace OpenMS
