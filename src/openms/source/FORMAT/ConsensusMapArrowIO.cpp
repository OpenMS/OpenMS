// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ConsensusMapArrowIO.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/ArrowIOHelpers.h>
#include "ArrowFileIOHelpers.h"
#include <OpenMS/FORMAT/ParquetFile.h>
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/FORMAT/ProteinIdentificationArrowIO.h>
#include <OpenMS/FORMAT/ModificationDefinitionIO.h>
#include <OpenMS/FORMAT/QPXFile.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/CHEMISTRY/ProForma.h>
#include <OpenMS/METADATA/ProteinRunReferences.h>

#include <arrow/api.h>
#include <arrow/builder.h>
#include <arrow/io/file.h>
#include <parquet/arrow/writer.h>
#include <parquet/arrow/reader.h>
#include <parquet/properties.h>

#include <cctype>
#include <filesystem>
#include <unordered_set>

namespace OpenMS
{

namespace // anonymous
{
  // ==================== Column header JSON helpers ====================

  /// Serialize ConsensusMap column headers to a JSON string (including metavalues).
  std::string serializeColumnHeaders_(const ConsensusMap::ColumnHeaders& headers)
  {
    std::string json = "[";
    bool first = true;
    for (const auto& [map_index, header] : headers)
    {
      if (!first) json += ",";
      json += "{\"map_index\":" + std::to_string(map_index)
            + ",\"filename\":\"" + ParquetFile::jsonEscape(header.filename) + "\""
            + ",\"label\":\"" + ParquetFile::jsonEscape(header.label) + "\""
            + ",\"size\":" + std::to_string(header.size)
            + ",\"unique_id\":" + std::to_string(header.unique_id)
            + ",\"metavalues\":" + Internal::ArrowFileIO::serializeMetaValues(header)
            + "}";
      first = false;
    }
    json += "]";
    return json;
  }

  /// Deserialize a JSON string into ConsensusMap column headers.
  ConsensusMap::ColumnHeaders deserializeColumnHeaders_(const std::string& json)
  {
    ConsensusMap::ColumnHeaders result;
    if (json.empty()) return result;

    size_t pos = 0;
    Internal::ArrowFileIO::skipWhitespace(json, pos);
    if (pos >= json.size() || json[pos] != '[') return result;
    ++pos;

    while (pos < json.size())
    {
      Internal::ArrowFileIO::skipWhitespace(json, pos);
      if (pos >= json.size() || json[pos] == ']') break;
      if (json[pos] == ',') { ++pos; continue; }
      if (json[pos] != '{') break;
      ++pos;

      ConsensusMap::ColumnHeader header;
      UInt64 map_index = 0;

      while (pos < json.size())
      {
        Internal::ArrowFileIO::skipWhitespace(json, pos);
        if (pos >= json.size() || json[pos] == '}') { ++pos; break; }
        if (json[pos] == ',') { ++pos; continue; }

        const size_t key_start = pos;
        std::string key = Internal::ArrowFileIO::parseJsonString(json, pos);
        if (pos == key_start) { pos = json.size(); break; } // no string key: malformed, give up
        Internal::ArrowFileIO::skipWhitespace(json, pos);
        if (pos < json.size() && json[pos] == ':') ++pos;
        Internal::ArrowFileIO::skipWhitespace(json, pos);

        if (key == "map_index" || key == "size" || key == "unique_id")
        {
          std::string num_str;
          while (pos < json.size() && (std::isdigit(json[pos]) || json[pos] == '-'))
          {
            num_str += json[pos++];
          }
          if (num_str.empty()) // not a number (null, string, ...): ignore, keep parsing
          {
            Internal::ArrowFileIO::skipJsonValue(json, pos);
            continue;
          }
          try
          {
            if (key == "map_index") map_index = std::stoull(num_str);
            else if (key == "size") header.size = std::stoull(num_str);
            else if (key == "unique_id") header.unique_id = std::stoull(num_str);
          }
          catch (...) {} // out of range or just "-": ignore
          Internal::ArrowFileIO::skipJsonValue(json, pos); // rest of a non-integer number such as 1.5
        }
        else if (key == "metavalues")
        {
          Internal::ArrowFileIO::parseMetaValueArray(json, pos, header);
        }
        else if (key == "filename")
        {
          header.filename = Internal::ArrowFileIO::parseJsonStringOrSkip(json, pos);
        }
        else if (key == "label")
        {
          header.label = Internal::ArrowFileIO::parseJsonStringOrSkip(json, pos);
        }
        else
        {
          Internal::ArrowFileIO::skipJsonValue(json, pos); // unknown key: ignore its value
        }
      }
      result[map_index] = header;
    }
    return result;
  }

  /// Read FeatureHandles from a list<struct{...}> column at a given row into a ConsensusFeature.
  void readHandles_(
    const std::shared_ptr<arrow::Array>& array,
    int64_t row,
    ConsensusFeature& cf)
  {
    if (!array || array->IsNull(row)) return;
    auto list_arr = std::static_pointer_cast<arrow::ListArray>(array);
    auto struct_arr = std::static_pointer_cast<arrow::StructArray>(list_arr->value_slice(row));
    if (!struct_arr || struct_arr->length() == 0) return;

    auto map_index_arr = std::static_pointer_cast<arrow::Int64Array>(struct_arr->field(0));
    auto unique_id_arr = std::static_pointer_cast<arrow::Int64Array>(struct_arr->field(1));
    auto rt_arr = std::static_pointer_cast<arrow::DoubleArray>(struct_arr->field(2));
    auto mz_arr = std::static_pointer_cast<arrow::DoubleArray>(struct_arr->field(3));
    auto intensity_arr = std::static_pointer_cast<arrow::FloatArray>(struct_arr->field(4));
    auto charge_arr = std::static_pointer_cast<arrow::Int32Array>(struct_arr->field(5));
    auto width_arr = std::static_pointer_cast<arrow::FloatArray>(struct_arr->field(6));

    for (int64_t i = 0; i < struct_arr->length(); ++i)
    {
      FeatureHandle handle;
      handle.setMapIndex(static_cast<UInt64>(map_index_arr->Value(i)));
      handle.setUniqueId(static_cast<UInt64>(unique_id_arr->Value(i)));
      handle.setRT(rt_arr->Value(i));
      handle.setMZ(mz_arr->Value(i));
      handle.setIntensity(intensity_arr->Value(i));
      handle.setCharge(static_cast<Int>(charge_arr->Value(i)));
      handle.setWidth(width_arr->Value(i));
      cf.insert(handle);
    }
  }

} // anonymous namespace


// ==================== Export methods ====================

std::shared_ptr<arrow::Table> ConsensusMapArrowIO::exportFeaturesToArrow(
  const ConsensusMap& cmap)
{
  // -- Scalar column builders --
  arrow::Int64Builder unique_id_builder;
  arrow::DoubleBuilder rt_builder, mz_builder;
  arrow::FloatBuilder intensity_builder;
  arrow::Int32Builder charge_builder;
  arrow::FloatBuilder quality_builder;
  arrow::FloatBuilder width_builder;

  // -- handles: list<struct{map_index, unique_id, rt, mz, intensity, charge, width}> --
  auto handle_map_index_b = std::make_shared<arrow::Int64Builder>();
  auto handle_unique_id_b = std::make_shared<arrow::Int64Builder>();
  auto handle_rt_b = std::make_shared<arrow::DoubleBuilder>();
  auto handle_mz_b = std::make_shared<arrow::DoubleBuilder>();
  auto handle_intensity_b = std::make_shared<arrow::FloatBuilder>();
  auto handle_charge_b = std::make_shared<arrow::Int32Builder>();
  auto handle_width_b = std::make_shared<arrow::FloatBuilder>();

  auto handle_struct_type = arrow::struct_({
    arrow::field("map_index", arrow::int64()),
    arrow::field("unique_id", arrow::int64()),
    arrow::field("rt", arrow::float64()),
    arrow::field("mz", arrow::float64()),
    arrow::field("intensity", arrow::float32()),
    arrow::field("charge", arrow::int32()),
    arrow::field("width", arrow::float32())
  });

  auto handle_struct_b = std::make_shared<arrow::StructBuilder>(
    handle_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{
      handle_map_index_b, handle_unique_id_b, handle_rt_b, handle_mz_b,
      handle_intensity_b, handle_charge_b, handle_width_b});
  arrow::ListBuilder handles_builder(arrow::default_memory_pool(), handle_struct_b);

  // -- metavalues: list<struct{name, value, value_type}> --
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

  Size num_rows = cmap.size();

  arrow::Status status;

  #define RESERVE_OR_RETURN(builder_name) \
    status = builder_name.Reserve(num_rows); \
    if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowIO: " #builder_name " Reserve failed: " << status.ToString() << std::endl; return nullptr; }

  RESERVE_OR_RETURN(unique_id_builder)
  RESERVE_OR_RETURN(rt_builder)
  RESERVE_OR_RETURN(mz_builder)
  RESERVE_OR_RETURN(intensity_builder)
  RESERVE_OR_RETURN(charge_builder)
  RESERVE_OR_RETURN(quality_builder)
  RESERVE_OR_RETURN(width_builder)

  #undef RESERVE_OR_RETURN

  // Exclude FWHM from metavalues (it is already stored as the dedicated 'width' column)
  static const std::unordered_set<std::string> excluded_mvs = {"FWHM"};

  for (const auto& cf : cmap)
  {
    (void)unique_id_builder.Append(static_cast<int64_t>(cf.getUniqueId()));
    (void)rt_builder.Append(cf.getRT());
    (void)mz_builder.Append(cf.getMZ());
    (void)intensity_builder.Append(cf.getIntensity());
    (void)charge_builder.Append(static_cast<int32_t>(cf.getCharge()));
    (void)quality_builder.Append(cf.getQuality());

    // Width: 0.0 is the unset default, stored as null to distinguish from explicitly set values.
    // On import, null is left as default (0.0), so round-trip is numerically correct.
    float w = cf.getWidth();
    if (w == 0.0f)
    {
      (void)width_builder.AppendNull();
    }
    else
    {
      (void)width_builder.Append(w);
    }

    // handles
    (void)handles_builder.Append();
    for (const auto& handle : cf.getFeatures())
    {
      (void)handle_struct_b->Append();
      (void)handle_map_index_b->Append(static_cast<int64_t>(handle.getMapIndex()));
      (void)handle_unique_id_b->Append(static_cast<int64_t>(handle.getUniqueId()));
      (void)handle_rt_b->Append(handle.getRT());
      (void)handle_mz_b->Append(handle.getMZ());
      (void)handle_intensity_b->Append(handle.getIntensity());
      (void)handle_charge_b->Append(static_cast<int32_t>(handle.getCharge()));
      (void)handle_width_b->Append(handle.getWidth());
    }

    // metavalues
    (void)metavalues_builder.Append();
    Internal::ArrowFileIO::appendMetaValues(cf, mv_name_b, mv_value_b, mv_type_b, mv_struct_b, excluded_mvs);
  }

  // Finalize arrays
  std::shared_ptr<arrow::Array> arr_unique_id, arr_rt, arr_mz, arr_intensity;
  std::shared_ptr<arrow::Array> arr_charge, arr_quality, arr_width;
  std::shared_ptr<arrow::Array> arr_handles, arr_metavalues;

  #define FINISH_OR_RETURN(builder_name, arr_name) \
    status = builder_name.Finish(&arr_name); \
    if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowIO: " #builder_name " Finish failed: " << status.ToString() << std::endl; return nullptr; }

  FINISH_OR_RETURN(unique_id_builder, arr_unique_id)
  FINISH_OR_RETURN(rt_builder, arr_rt)
  FINISH_OR_RETURN(mz_builder, arr_mz)
  FINISH_OR_RETURN(intensity_builder, arr_intensity)
  FINISH_OR_RETURN(charge_builder, arr_charge)
  FINISH_OR_RETURN(quality_builder, arr_quality)
  FINISH_OR_RETURN(width_builder, arr_width)
  FINISH_OR_RETURN(handles_builder, arr_handles)
  FINISH_OR_RETURN(metavalues_builder, arr_metavalues)

  #undef FINISH_OR_RETURN

  // Build schema from registry
  auto schema = ConsensusFeatureSchema::schema();

  auto table = arrow::Table::Make(schema, {
    arr_unique_id, arr_rt, arr_mz, arr_intensity,
    arr_charge, arr_quality, arr_width,
    arr_handles, arr_metavalues
  });

  // Validate table against registry schema
  auto validation = ArrowSchemaValidation::validate(table, ConsensusFeatureSchema::schema());
  if (!validation.valid)
  {
    OPENMS_LOG_ERROR << "Schema validation failed: " << validation.toString() << "\n";
    return nullptr;
  }

  return table;
}

std::shared_ptr<arrow::Table> ConsensusMapArrowIO::exportPSMsToArrow(
  const ConsensusMap& cmap)
{
  // Collect all PeptideIdentifications with their associated consensus_feature_unique_ids.
  PeptideIdentificationList all_pep_ids;
  std::vector<std::pair<int64_t, bool>> feature_ids_per_pep_id; // (id, is_null)

  // Collect from consensus features (flat - no subordinates)
  for (const auto& cf : cmap)
  {
    for (const auto& pep_id : cf.getPeptideIdentifications())
    {
      all_pep_ids.push_back(pep_id);
      feature_ids_per_pep_id.push_back({static_cast<int64_t>(cf.getUniqueId()), false});
    }
  }

  // Unassigned peptide identifications get consensus_feature_unique_id = null
  for (const auto& pep_id : cmap.getUnassignedPeptideIdentifications())
  {
    all_pep_ids.push_back(pep_id);
    feature_ids_per_pep_id.push_back({0, true});
  }

  // Delegate to QPXFile to produce the base PSM table
  auto base_table = QPXFile::exportToArrow(
    cmap.getProteinIdentifications(), all_pep_ids, true);
  if (!base_table) { return nullptr; }

  // Build the consensus_feature_unique_id column (one value per PeptideHit)
  arrow::Int64Builder feature_id_builder;

  for (size_t i = 0; i < all_pep_ids.size(); ++i)
  {
    if (all_pep_ids[i].getHits().empty())
    {
      continue;
    }

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
    OPENMS_LOG_ERROR << "ConsensusMapArrowIO: feature_id_builder Finish failed: " << status.ToString() << std::endl;
    return nullptr;
  }

  if (feature_id_array->length() != base_table->num_rows())
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowIO: consensus_feature_unique_id column length (" << feature_id_array->length()
                     << ") does not match table row count (" << base_table->num_rows() << ")" << std::endl;
    return nullptr;
  }

  auto chunked_feature_id = std::make_shared<arrow::ChunkedArray>(feature_id_array);
  auto result = base_table->AddColumn(0, arrow::field("consensus_feature_unique_id", arrow::int64()), chunked_feature_id);
  if (!result.ok())
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowIO: AddColumn failed: " << result.status().ToString() << std::endl;
    return nullptr;
  }
  return *result;
}

bool ConsensusMapArrowIO::exportToParquet(
  const ConsensusMap& cmap,
  const std::string& directory,
  const ParquetWriteConfig& config)
{
  // Mirror XMLHandler::checkUniqueIdentifiers_ — fail before any file is opened
  // so we never leave a partial .consensusparquet behind. Throws Exception::InvalidValue.
  ProteinIdentificationArrowIO::checkUniqueIdentifiers(cmap.getProteinIdentifications());
  // every peptide identification needs its protein identification run
  ProteinRunReferences::check(cmap);

  // 1. Create output directory
  try
  {
    std::filesystem::create_directories(std::string(directory));
  }
  catch (const std::filesystem::filesystem_error& e)
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowIO: Failed to create directory: " << e.what() << std::endl;
    return false;
  }

  // 2. Export consensus features table (with ConsensusMap-level metadata)
  auto features_table = exportFeaturesToArrow(cmap);
  if (!features_table)
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowIO: Failed to create consensus features Arrow table" << std::endl;
    return false;
  }

  std::unordered_map<std::string, std::string> cmap_metadata;
  cmap_metadata["column_headers"] = serializeColumnHeaders_(cmap.getColumnHeaders());
  cmap_metadata["experiment_type"] = cmap.getExperimentType();
  cmap_metadata["document_id"] = cmap.getIdentifier();
  cmap_metadata["loaded_file_path"] = cmap.getLoadedFilePath();
  cmap_metadata["loaded_file_type"] = FileTypes::typeToName(cmap.getLoadedFileType());
  cmap_metadata["data_processing"] = Internal::ArrowFileIO::serializeDataProcessing(cmap.getDataProcessing());
  cmap_metadata["unique_id"] = std::to_string(cmap.getUniqueId());
  cmap_metadata["cmap_metavalues"] = Internal::ArrowFileIO::serializeMetaValues(cmap);

  if (!Internal::ArrowFileIO::writeParquetTable(features_table, directory + "/consensus_features.parquet",
                                  "consensus_features", config, "ConsensusMapArrowIO", cmap_metadata))
  {
    return false;
  }

  // 3. Export PSMs table
  auto psms_table = exportPSMsToArrow(cmap);
  if (!psms_table)
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowIO: Failed to create PSMs Arrow table" << std::endl;
    return false;
  }
  if (!Internal::ArrowFileIO::writeParquetTable(psms_table, directory + "/psms.parquet", "psms", config, "ConsensusMapArrowIO"))
  {
    return false;
  }

  // 4. Delegate protein data to ProteinIdentificationArrowIO
  const auto& prot_ids = cmap.getProteinIdentifications();
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
          ModificationDefinitionIO::encodeByRun(prot_ids, ModificationDefinitionIO::collect(cmap))))
  {
    return false;
  }

  return true;
}

// ==================== Import methods ====================

bool ConsensusMapArrowIO::importFeaturesFromArrow(
  const std::shared_ptr<arrow::Table>& table,
  ConsensusMap& cmap)
{
  if (!table)
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowIO: null table passed to importFeaturesFromArrow" << std::endl;
    return false;
  }

  auto combined_result = table->CombineChunks(arrow::default_memory_pool());
  if (!combined_result.ok())
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowIO: Failed to combine chunks" << std::endl;
    return false;
  }
  const auto& tbl = *combined_result;

  int64_t num_rows = tbl->num_rows();
  if (num_rows == 0)
  {
    return true;
  }

  // Validate table schema against registry (subset mode — file may have extra columns)
  auto validation = ArrowSchemaValidation::validate(tbl, ConsensusFeatureSchema::schema(), ArrowSchemaValidation::Mode::Subset);
  if (!validation.valid)
  {
    OPENMS_LOG_ERROR << "Incompatible schema: " << validation.toString() << "\n";
    return false;
  }

  auto col_unique_id = ArrowIOHelpers::getColumn(tbl, ConsensusFeatureSchema::UNIQUE_ID, true, "ConsensusMapArrowIO");
  auto col_rt = ArrowIOHelpers::getColumn(tbl, ConsensusFeatureSchema::RT, true, "ConsensusMapArrowIO");
  auto col_mz = ArrowIOHelpers::getColumn(tbl, ConsensusFeatureSchema::MZ, true, "ConsensusMapArrowIO");
  auto col_intensity = ArrowIOHelpers::getColumn(tbl, ConsensusFeatureSchema::INTENSITY, true, "ConsensusMapArrowIO");
  auto col_charge = ArrowIOHelpers::getColumn(tbl, ConsensusFeatureSchema::CHARGE, true, "ConsensusMapArrowIO");
  auto col_quality = ArrowIOHelpers::getColumn(tbl, ConsensusFeatureSchema::QUALITY, true, "ConsensusMapArrowIO");
  auto col_width = ArrowIOHelpers::getColumn(tbl, ConsensusFeatureSchema::WIDTH, /*required=*/false, "ConsensusMapArrowIO");
  auto col_handles = ArrowIOHelpers::getColumn(tbl, ConsensusFeatureSchema::HANDLES, /*required=*/false, "ConsensusMapArrowIO");
  auto col_metavalues = ArrowIOHelpers::getColumn(tbl, ConsensusFeatureSchema::METAVALUES, /*required=*/false, "ConsensusMapArrowIO");

  if (!col_unique_id || !col_rt || !col_mz || !col_intensity || !col_charge || !col_quality)
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowIO: Missing one or more required columns" << std::endl;
    return false;
  }

  for (int64_t i = 0; i < num_rows; ++i)
  {
    ConsensusFeature cf;
    cf.setUniqueId(static_cast<UInt64>(ArrowIOHelpers::getInt64Value(col_unique_id, i, 0)));
    cf.setRT(ArrowIOHelpers::getDoubleValue(col_rt, i));
    cf.setMZ(ArrowIOHelpers::getDoubleValue(col_mz, i));
    cf.setIntensity(ArrowIOHelpers::getFloatValue(col_intensity, i));
    cf.setCharge(static_cast<Int>(ArrowIOHelpers::getInt32Value(col_charge, i)));
    cf.setQuality(ArrowIOHelpers::getFloatValue(col_quality, i));

    // Width: null means unset (default 0.0), so we only call setWidth for non-null values.
    if (col_width && !ArrowIOHelpers::isNull(col_width, i))
    {
      cf.setWidth(ArrowIOHelpers::getFloatValue(col_width, i));
    }

    if (col_handles)
    {
      readHandles_(col_handles, i, cf);
    }

    if (col_metavalues)
    {
      ArrowIOHelpers::readMetaValues(col_metavalues, i, cf);
    }

    cmap.push_back(std::move(cf));
  }

  return true;
}

bool ConsensusMapArrowIO::importPSMsFromArrow(
  const std::shared_ptr<arrow::Table>& table,
  ConsensusMap& cmap)
{
  if (!table || table->num_rows() == 0) { return true; }

  auto combined_result = table->CombineChunks(arrow::default_memory_pool());
  if (!combined_result.ok())
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowIO: Failed to combine chunks in importPSMsFromArrow" << std::endl;
    return false;
  }
  const auto& tbl = *combined_result;

  // Read consensus_feature_unique_id and PEPTIDE_IDENTIFICATION_INDEX side-by-side
  // so we can route each constructed PeptideIdentification to the right ConsensusFeature.
  auto col_feature_id = ArrowIOHelpers::getColumn(tbl, "consensus_feature_unique_id", true, "ConsensusMapArrowIO");
  auto col_p_id = ArrowIOHelpers::getColumn(tbl, PSMSchema::PEPTIDE_IDENTIFICATION_INDEX, true, "ConsensusMapArrowIO");
  if (!col_feature_id || !col_p_id)
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowIO: Missing required columns for PSM import" << std::endl;
    return false;
  }
  if (col_feature_id->type_id() != arrow::Type::INT64 ||
      col_p_id->type_id() != arrow::Type::INT32)
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowIO: Invalid column types for PSM import "
                     << "(expected consensus_feature_unique_id=int64 and "
                     << PSMSchema::PEPTIDE_IDENTIFICATION_INDEX << "=int32)" << std::endl;
    return false;
  }

  // Build (p_id -> feature_unique_id-or-null) map by capturing the first row each p_id appears.
  // Do not assume rows of the same p_id are contiguous.
  std::unordered_map<int32_t, std::pair<bool, int64_t>> p_id_to_feature; // bool: is_null
  for (int64_t row = 0; row < tbl->num_rows(); ++row)
  {
    int32_t p_id = ArrowIOHelpers::getInt32Value(col_p_id, row, -1);
    bool is_null = ArrowIOHelpers::isNull(col_feature_id, row);
    int64_t fid = is_null ? 0 : ArrowIOHelpers::getInt64Value(col_feature_id, row, 0);
    auto [it, inserted] = p_id_to_feature.try_emplace(p_id, std::make_pair(is_null, fid));
    if (!inserted)
    {
      const auto& [prev_null, prev_fid] = it->second;
      if (prev_null != is_null || (!is_null && prev_fid != fid))
      {
        OPENMS_LOG_ERROR << "ConsensusMapArrowIO: inconsistent consensus_feature_unique_id "
                         << "mappings for p_id " << p_id << " (first: "
                         << (prev_null ? "null" : std::to_string(prev_fid))
                         << ", later: " << (is_null ? "null" : std::to_string(fid))
                         << ")" << std::endl;
        return false;
      }
    }
  }

  // Delegate PSM construction to QPXFile::importFromArrow.
  std::vector<ProteinIdentification> prot_ids = cmap.getProteinIdentifications();
  PeptideIdentificationList pep_ids;
  if (!QPXFile::importFromArrow(tbl, prot_ids, pep_ids)) { return false; }
  cmap.setProteinIdentifications(prot_ids);

  // Walk each constructed PeptideIdentification and route to features by feature_id.
  std::unordered_map<int64_t, ConsensusFeature*> feature_lookup;
  for (auto& cf : cmap)
  {
    feature_lookup[static_cast<int64_t>(cf.getUniqueId())] = &cf;
  }

  // Re-derive the p_id sequence in first-seen order (independent of adjacency).
  std::vector<int32_t> p_id_order;
  std::unordered_set<int32_t> seen;
  seen.reserve(static_cast<size_t>(tbl->num_rows()));
  for (int64_t row = 0; row < tbl->num_rows(); ++row)
  {
    int32_t p_id = ArrowIOHelpers::getInt32Value(col_p_id, row, -1);
    if (seen.insert(p_id).second) { p_id_order.push_back(p_id); }
  }
  if (p_id_order.size() != pep_ids.size())
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowIO: PSM count mismatch after import "
                     << "(rows=" << p_id_order.size() << ", pep_ids=" << pep_ids.size() << ")" << std::endl;
    return false;
  }

  for (size_t i = 0; i < pep_ids.size(); ++i)
  {
    auto fit = p_id_to_feature.find(p_id_order[i]);
    bool is_null = (fit == p_id_to_feature.end()) ? true : fit->second.first;
    int64_t fid = (fit == p_id_to_feature.end()) ? 0 : fit->second.second;

    if (is_null)
    {
      cmap.getUnassignedPeptideIdentifications().push_back(std::move(pep_ids[i]));
    }
    else
    {
      auto it = feature_lookup.find(fid);
      if (it != feature_lookup.end())
      {
        it->second->getPeptideIdentifications().push_back(std::move(pep_ids[i]));
      }
      else
      {
        OPENMS_LOG_WARN << "ConsensusMapArrowIO: Could not find consensus feature with id "
                        << fid << " for PSM. Adding as unassigned." << std::endl;
        cmap.getUnassignedPeptideIdentifications().push_back(std::move(pep_ids[i]));
      }
    }
  }

  return true;
}

bool ConsensusMapArrowIO::importFromParquet(
  const std::string& directory,
  ConsensusMap& cmap)
{
  cmap = ConsensusMap{};

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
  cmap.setProteinIdentifications(prot_ids);

  // 2. Import consensus features and ConsensusMap-level metadata
  auto features_table = Internal::ArrowFileIO::readParquetTable(directory + "/consensus_features.parquet", "ConsensusMapArrowIO");
  if (!features_table) { return false; }

  // Extract ConsensusMap-level metadata from file-level key-value metadata
  auto schema_md = features_table->schema()->metadata();
  if (schema_md)
  {
    int idx = schema_md->FindKey("column_headers");
    if (idx >= 0) cmap.setColumnHeaders(deserializeColumnHeaders_(schema_md->value(idx)));

    idx = schema_md->FindKey("experiment_type");
    if (idx >= 0) cmap.setExperimentType(schema_md->value(idx));

    idx = schema_md->FindKey("document_id");
    if (idx >= 0) cmap.setIdentifier(schema_md->value(idx));

    idx = schema_md->FindKey("loaded_file_path");
    if (idx >= 0) cmap.setLoadedFilePath(schema_md->value(idx));

    idx = schema_md->FindKey("data_processing");
    if (idx >= 0)
    {
      cmap.setDataProcessing(Internal::ArrowFileIO::deserializeDataProcessing(schema_md->value(idx), "ConsensusMapArrowIO"));
    }

    idx = schema_md->FindKey("unique_id");
    if (idx >= 0)
    {
      try { cmap.setUniqueId(std::stoull(schema_md->value(idx))); }
      catch (...) {}
    }

    idx = schema_md->FindKey("cmap_metavalues");
    if (idx >= 0)
    {
      Internal::ArrowFileIO::deserializeMetaValues(schema_md->value(idx), cmap);
    }
  }

  if (!importFeaturesFromArrow(features_table, cmap))
  {
    return false;
  }

  // 3. Import PSMs (links to consensus features already in the map)
  auto psms_table = Internal::ArrowFileIO::readParquetTable(directory + "/psms.parquet", "ConsensusMapArrowIO");
  if (!psms_table) { return false; }
  const Size stored_runs = cmap.getProteinIdentifications().size();
  if (!importPSMsFromArrow(psms_table, cmap))
  {
    return false;
  }
  // every peptide identification needs its protein identification run: the PSM table import appends a shell run
  // for each run identifier the stored runs lack, which is a dangling reference in this file
  if (cmap.getProteinIdentifications().size() > stored_runs)
  {
    throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      ProteinRunReferences::missingRunMessage(cmap.getProteinIdentifications()[stored_runs].getIdentifier()));
  }
  ProteinRunReferences::check(cmap);

  // 4. Synthesize fresh ProtID identifiers + apply rename to every pep_id collection
  //    we own (per-consensus-feature + unassigned). Mirrors IdXMLFile.cpp:530 — the
  //    stored identifier becomes informational; the in-memory identifier downstream
  //    sees is freshly synthesized with a UniqueIdGenerator suffix.
  {
    auto& prot_ids = cmap.getProteinIdentifications();
    auto rename = ProteinIdentificationArrowIO::synthesizeRunIdentifiers(prot_ids);

    for (auto& cf : cmap)
    {
      ProteinIdentificationArrowIO::applyRunIdentifierRename(
          rename, cf.getPeptideIdentifications());
    }
    ProteinIdentificationArrowIO::applyRunIdentifierRename(
        rename, cmap.getUnassignedPeptideIdentifications());
  }

  return true;
}

} // namespace OpenMS
