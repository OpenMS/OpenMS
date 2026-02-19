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
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/FORMAT/ProteinIdentificationArrowIO.h>

#include <arrow/api.h>
#include <arrow/builder.h>
#include <arrow/io/file.h>
#include <parquet/arrow/writer.h>
#include <parquet/arrow/reader.h>
#include <parquet/properties.h>

#include <cstdio>
#include <cstring>
#include <filesystem>
#include <random>
#include <unordered_set>

namespace OpenMS
{

namespace // anonymous
{
  // ==================== Export helpers ====================

  /// Append all MetaValues (excluding specified keys) from a MetaInfoInterface to the struct builders.
  void appendMetaValues_(
    const MetaInfoInterface& mii,
    std::shared_ptr<arrow::StringBuilder>& name_b,
    std::shared_ptr<arrow::StringBuilder>& value_b,
    std::shared_ptr<arrow::StringBuilder>& type_b,
    std::shared_ptr<arrow::StructBuilder>& struct_b,
    const std::unordered_set<std::string>& excluded_keys)
  {
    std::vector<String> keys;
    mii.getKeys(keys);
    for (const auto& key : keys)
    {
      if (excluded_keys.count(key)) continue;
      const DataValue& val = mii.getMetaValue(key);
      (void)struct_b->Append();
      (void)name_b->Append(key);
      (void)value_b->Append(val.toString());
      switch (val.valueType())
      {
        case DataValue::INT_VALUE: (void)type_b->Append("int"); break;
        case DataValue::DOUBLE_VALUE: (void)type_b->Append("float"); break;
        case DataValue::STRING_VALUE: (void)type_b->Append("str"); break;
        default: (void)type_b->Append("str"); break;
      }
    }
  }

  /// Write an Arrow table to a Parquet file with QPX-style metadata.
  bool writeArrowTableToParquet_(
    std::shared_ptr<arrow::Table> table,
    const String& filename,
    const std::string& file_type,
    const ParquetWriteConfig& config)
  {
    if (!table)
    {
      OPENMS_LOG_ERROR << "FeatureMapArrowIO: null table passed to writeArrowTableToParquet_" << std::endl;
      return false;
    }

    // Generate RFC 4122 version-4 UUID
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint32_t> dist;
    uint8_t bytes[16];
    for (int i = 0; i < 4; ++i)
    {
      uint32_t r = dist(gen);
      std::memcpy(bytes + i * 4, &r, 4);
    }
    bytes[6] = (bytes[6] & 0x0F) | 0x40; // version 4
    bytes[8] = (bytes[8] & 0x3F) | 0x80; // variant 1
    char buf[37];
    std::snprintf(buf, sizeof(buf),
      "%02x%02x%02x%02x-%02x%02x-%02x%02x-%02x%02x-%02x%02x%02x%02x%02x%02x",
      bytes[0], bytes[1], bytes[2], bytes[3],
      bytes[4], bytes[5], bytes[6], bytes[7],
      bytes[8], bytes[9], bytes[10], bytes[11],
      bytes[12], bytes[13], bytes[14], bytes[15]);
    std::string uuid_str(buf);

    auto metadata = arrow::key_value_metadata({
      {"qpx_version", "1.0"},
      {"creator", "OpenMS"},
      {"file_type", file_type},
      {"creation_date", DateTime::nowUTC().toString("yyyy-MM-ddThh:mm:ssZ")},
      {"uuid", uuid_str},
      {"software_provider", "OpenMS"}
    });
    table = table->ReplaceSchemaMetadata(metadata);

    // Open output file
    auto result = arrow::io::FileOutputStream::Open(filename);
    if (!result.ok())
    {
      OPENMS_LOG_ERROR << "FeatureMapArrowIO: Failed to open file: " << filename << std::endl;
      return false;
    }
    auto outfile = *result;

    // Configure Parquet writer
    auto builder = parquet::WriterProperties::Builder();

    switch (config.compression)
    {
      case ParquetWriteConfig::Compression::NONE:
        builder.compression(arrow::Compression::UNCOMPRESSED);
        break;
      case ParquetWriteConfig::Compression::SNAPPY:
        builder.compression(arrow::Compression::SNAPPY);
        break;
      case ParquetWriteConfig::Compression::GZIP:
        builder.compression(arrow::Compression::GZIP);
        builder.compression_level(config.compression_level);
        break;
      case ParquetWriteConfig::Compression::LZ4:
        builder.compression(arrow::Compression::LZ4);
        break;
      case ParquetWriteConfig::Compression::ZSTD:
        builder.compression(arrow::Compression::ZSTD);
        builder.compression_level(config.compression_level);
        break;
    }

    builder.data_pagesize(config.data_page_size);

    if (config.write_statistics)
    {
      builder.enable_statistics();
    }
    else
    {
      builder.disable_statistics();
    }

    auto writer_props = builder.build();
    auto arrow_props = parquet::ArrowWriterProperties::Builder().store_schema()->build();

    // Write
    auto write_status = parquet::arrow::WriteTable(
      *table, arrow::default_memory_pool(), outfile,
      config.row_group_size, writer_props, arrow_props);

    if (!write_status.ok())
    {
      OPENMS_LOG_ERROR << "FeatureMapArrowIO: Failed to write Parquet: "
                       << write_status.ToString() << std::endl;
      return false;
    }

    return true;
  }

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

  // No metavalue keys to exclude for features
  static const std::unordered_set<std::string> excluded_mvs = {};

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

    // === width (nullable: null if 0) ===
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
    appendMetaValues_(feature, mv_name_b, mv_value_b, mv_type_b, mv_struct_b, excluded_mvs);

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

  // Build schema (17 columns)
  auto schema = arrow::schema({
    arrow::field("unique_id", arrow::int64(), /*nullable=*/false),
    arrow::field("parent_feature_id", arrow::int64(), /*nullable=*/true),
    arrow::field("depth", arrow::int32(), /*nullable=*/false),
    arrow::field("rt", arrow::float64(), /*nullable=*/false),
    arrow::field("mz", arrow::float64(), /*nullable=*/false),
    arrow::field("intensity", arrow::float32(), /*nullable=*/false),
    arrow::field("charge", arrow::int32(), /*nullable=*/false),
    arrow::field("overall_quality", arrow::float32(), /*nullable=*/false),
    arrow::field("quality_rt", arrow::float32(), /*nullable=*/false),
    arrow::field("quality_mz", arrow::float32(), /*nullable=*/false),
    arrow::field("width", arrow::float32(), /*nullable=*/true),
    arrow::field("rt_bb_min", arrow::float64(), /*nullable=*/true),
    arrow::field("rt_bb_max", arrow::float64(), /*nullable=*/true),
    arrow::field("mz_bb_min", arrow::float64(), /*nullable=*/true),
    arrow::field("mz_bb_max", arrow::float64(), /*nullable=*/true),
    arrow::field("convex_hulls", convex_hulls_builder.type(), /*nullable=*/false),
    arrow::field("metavalues", metavalues_builder.type(), /*nullable=*/false),
  });

  auto table = arrow::Table::Make(schema, {
    arr_unique_id, arr_parent_id, arr_depth,
    arr_rt, arr_mz, arr_intensity,
    arr_charge, arr_overall_quality,
    arr_quality_rt, arr_quality_mz,
    arr_width,
    arr_rt_bb_min, arr_rt_bb_max, arr_mz_bb_min, arr_mz_bb_max,
    arr_convex_hulls, arr_metavalues
  });

  return table;
}

std::shared_ptr<arrow::Table> FeatureMapArrowIO::exportPSMsToArrow(
  const FeatureMap& /*feature_map*/)
{
  return nullptr;
}

bool FeatureMapArrowIO::exportToParquet(
  const FeatureMap& /*feature_map*/,
  const String& /*directory*/,
  const ParquetWriteConfig& /*config*/)
{
  return false;
}

bool FeatureMapArrowIO::importFeaturesFromArrow(
  const std::shared_ptr<arrow::Table>& /*table*/,
  FeatureMap& /*feature_map*/)
{
  return false;
}

bool FeatureMapArrowIO::importPSMsFromArrow(
  const std::shared_ptr<arrow::Table>& /*table*/,
  FeatureMap& /*feature_map*/)
{
  return false;
}

bool FeatureMapArrowIO::importFromParquet(
  const String& /*directory*/,
  FeatureMap& /*feature_map*/)
{
  return false;
}

} // namespace OpenMS

#endif // WITH_PARQUET
