// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ProteinIdentificationArrowIO.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/FORMAT/ArrowIOHelpers.h>
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>

#include <arrow/api.h>
#include <arrow/builder.h>
#include <arrow/io/file.h>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>
#include <parquet/properties.h>

#include <cstdio>
#include <cstring>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace OpenMS
{

namespace // anonymous
{
  // ==================== Export helpers ====================

  /// Append all MetaValues (excluding specified keys) from a MetaInfoInterface to the struct builders.
  /// @param mii The MetaInfoInterface to read keys from
  /// @param name_b StringBuilder for the name field
  /// @param value_b StringBuilder for the value field (toString)
  /// @param type_b StringBuilder for the value_type field
  /// @param struct_b StructBuilder for the list element struct
  /// @param excluded_keys Set of key names to skip
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
        case DataValue::DOUBLE_VALUE: (void)type_b->Append("double"); break;
        case DataValue::STRING_VALUE: (void)type_b->Append("string"); break;
        case DataValue::INT_LIST: (void)type_b->Append("int_list"); break;
        case DataValue::DOUBLE_LIST: (void)type_b->Append("double_list"); break;
        case DataValue::STRING_LIST: (void)type_b->Append("string_list"); break;
        default: (void)type_b->Append("string"); break;
      }
    }
  }

  /// Write an Arrow table to a Parquet file with QPX-style metadata.
  /// @param table The Arrow table to write
  /// @param filename Output file path
  /// @param file_type The file_type metadata value (e.g. "proteins", "protein_groups", "search_params")
  /// @param config Parquet writing options
  /// @return true on success, false on error
  bool writeArrowTableToParquet_(
    std::shared_ptr<arrow::Table> table,
    const String& filename,
    const std::string& file_type,
    const ParquetWriteConfig& config)
  {
    if (!table)
    {
      OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: null table passed to writeArrowTableToParquet_" << std::endl;
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
      OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: Failed to open file: " << filename << std::endl;
      return false;
    }
    const auto& outfile = *result;

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
      OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: Failed to write Parquet: "
                       << write_status.ToString() << std::endl;
      return false;
    }

    return true;
  }

  // ==================== Import helpers ====================

  /// Read a single Parquet file into an Arrow table.
  std::shared_ptr<arrow::Table> readParquetTable_(const String& filename)
  {
    auto infile_result = arrow::io::ReadableFile::Open(std::string(filename));
    if (!infile_result.ok())
    {
      OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: Failed to open file: " << filename << std::endl;
      return nullptr;
    }
    const auto& infile = *infile_result;

    auto reader_result = parquet::arrow::OpenFile(infile, arrow::default_memory_pool());
    if (!reader_result.ok())
    {
      OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: Failed to create Parquet reader for: " << filename << std::endl;
      return nullptr;
    }
    auto reader = std::move(reader_result.ValueOrDie());

    std::shared_ptr<arrow::Table> table;
    auto read_status = reader->ReadTable(&table);
    if (!read_status.ok())
    {
      OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: Failed to read table: " << read_status.ToString() << std::endl;
      return nullptr;
    }

    auto combined = table->CombineChunks(arrow::default_memory_pool());
    if (!combined.ok())
    {
      OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: Failed to combine chunks" << std::endl;
      return nullptr;
    }

    return *combined;
  }

  /// Fetch a named column from a table, combining chunks if needed.
  /// Returns nullptr if column not found (logs error if required).
  std::shared_ptr<arrow::Array> getColumn_(
    const std::shared_ptr<arrow::Table>& table,
    const std::string& name,
    bool required = true)
  {
    auto column = table->GetColumnByName(name);
    if (!column)
    {
      if (required)
      {
        OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: Missing required column '" << name << "'" << std::endl;
      }
      return nullptr;
    }
    if (column->num_chunks() == 0)
    {
      OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: Column '" << name << "' has no chunks" << std::endl;
      return nullptr;
    }
    if (column->num_chunks() == 1)
    {
      return column->chunk(0);
    }
    auto combined = arrow::Concatenate(column->chunks(), arrow::default_memory_pool());
    if (!combined.ok())
    {
      OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: Failed to combine chunks for column '" << name << "'" << std::endl;
      return nullptr;
    }
    return *combined;
  }

  /// Get string value at a row, returning empty string if null.
  String getStringValue_(const std::shared_ptr<arrow::Array>& array, int64_t row)
  {
    if (!array || array->IsNull(row)) return "";
    return std::static_pointer_cast<arrow::StringArray>(array)->GetString(row);
  }

  /// Get double value at a row, returning default_val if null.
  double getDoubleValue_(const std::shared_ptr<arrow::Array>& array, int64_t row, double default_val = 0.0)
  {
    if (!array || array->IsNull(row)) return default_val;
    return std::static_pointer_cast<arrow::DoubleArray>(array)->Value(row);
  }

  /// Get int32 value at a row, returning default_val if null.
  int32_t getInt32Value_(const std::shared_ptr<arrow::Array>& array, int64_t row, int32_t default_val = 0)
  {
    if (!array || array->IsNull(row)) return default_val;
    return std::static_pointer_cast<arrow::Int32Array>(array)->Value(row);
  }

  /// Get boolean value at a row, returning default_val if null.
  bool getBoolValue_(const std::shared_ptr<arrow::Array>& array, int64_t row, bool default_val = false)
  {
    if (!array || array->IsNull(row)) return default_val;
    return std::static_pointer_cast<arrow::BooleanArray>(array)->Value(row);
  }

  /// Check if value at row is null.
  bool isNull_(const std::shared_ptr<arrow::Array>& array, int64_t row)
  {
    return !array || array->IsNull(row);
  }

  /// Read a list<utf8> column at a given row into a vector of Strings.
  std::vector<String> readStringList_(const std::shared_ptr<arrow::Array>& array, int64_t row)
  {
    std::vector<String> result;
    if (!array || array->IsNull(row)) return result;
    auto list_arr = std::static_pointer_cast<arrow::ListArray>(array);
    auto values = std::static_pointer_cast<arrow::StringArray>(list_arr->value_slice(row));
    result.reserve(values->length());
    for (int64_t i = 0; i < values->length(); ++i)
    {
      result.emplace_back(values->GetString(i));
    }
    return result;
  }

  /// Build a map from run_identifier to index in the protein_identifications vector.
  std::unordered_map<std::string, size_t> buildRunIdMap_(
    const std::vector<ProteinIdentification>& protein_identifications)
  {
    std::unordered_map<std::string, size_t> map;
    for (size_t i = 0; i < protein_identifications.size(); ++i)
    {
      map[protein_identifications[i].getIdentifier()] = i;
    }
    return map;
  }

} // anonymous namespace


// ==================== Export methods ====================

std::shared_ptr<arrow::Table> ProteinIdentificationArrowIO::exportProteinsToArrow(
  const std::vector<ProteinIdentification>& protein_identifications)
{
  // -- Simple column builders --
  arrow::StringBuilder accession_builder, sequence_builder;
  arrow::StringBuilder description_builder, run_identifier_builder;
  arrow::DoubleBuilder score_builder, coverage_builder;
  arrow::Int32Builder rank_builder;
  arrow::BooleanBuilder is_decoy_builder;

  // -- modifications list<struct{position: int32, modification: utf8}> --
  auto mod_position_b = std::make_shared<arrow::Int32Builder>();
  auto mod_name_b = std::make_shared<arrow::StringBuilder>();
  auto mod_struct_type = arrow::struct_({
    arrow::field("position", arrow::int32()),
    arrow::field("modification", arrow::utf8())
  });
  auto mod_struct_b = std::make_shared<arrow::StructBuilder>(
    mod_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{mod_position_b, mod_name_b});
  arrow::ListBuilder modifications_builder(arrow::default_memory_pool(), mod_struct_b);

  // -- metavalues list<struct{name: utf8, value: utf8, value_type: utf8}> --
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

  // Estimate total rows for capacity reservation
  Size num_rows = 0;
  for (const auto& prot_id : protein_identifications)
  {
    num_rows += prot_id.getHits().size();
  }

  // Reserve capacity for simple builders
  arrow::Status status;

  status = accession_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: accession_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = score_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: score_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = rank_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: rank_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = coverage_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: coverage_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = sequence_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: sequence_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = description_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: description_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = is_decoy_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: is_decoy_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = run_identifier_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: run_identifier_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }

  // Metavalue keys excluded from metavalues column (they have dedicated columns)
  static const std::unordered_set<std::string> excluded_hit_mvs = {
    "Description", "target_decoy"
  };

  for (const auto& prot_id : protein_identifications)
  {
    const String& run_id = prot_id.getIdentifier();

    for (const auto& hit : prot_id.getHits())
    {
      // === accession (not nullable) ===
      (void)accession_builder.Append(hit.getAccession());

      // === score (not nullable) ===
      (void)score_builder.Append(hit.getScore());

      // === rank (nullable) ===
      (void)rank_builder.Append(static_cast<int32_t>(hit.getRank()));

      // === coverage (nullable, null if < 0 i.e. COVERAGE_UNKNOWN) ===
      double cov = hit.getCoverage();
      if (cov < 0)
      {
        (void)coverage_builder.AppendNull();
      }
      else
      {
        (void)coverage_builder.Append(cov);
      }

      // === sequence (nullable, null if empty) ===
      const String& seq = hit.getSequence();
      if (seq.empty())
      {
        (void)sequence_builder.AppendNull();
      }
      else
      {
        (void)sequence_builder.Append(seq);
      }

      // === description (nullable, null if empty or not set) ===
      if (hit.metaValueExists("Description"))
      {
        String desc = hit.getDescription();
        if (desc.empty())
        {
          (void)description_builder.AppendNull();
        }
        else
        {
          (void)description_builder.Append(desc);
        }
      }
      else
      {
        (void)description_builder.AppendNull();
      }

      // === is_decoy (nullable) ===
      {
        auto td_type = hit.getTargetDecoyType();
        if (td_type == ProteinHit::TargetDecoyType::UNKNOWN)
        {
          (void)is_decoy_builder.AppendNull();
        }
        else
        {
          (void)is_decoy_builder.Append(td_type == ProteinHit::TargetDecoyType::DECOY);
        }
      }

      // === run_identifier (not nullable) ===
      (void)run_identifier_builder.Append(run_id);

      // === modifications (nullable) ===
      const auto& mods = hit.getModifications();
      if (mods.empty())
      {
        (void)modifications_builder.AppendNull();
      }
      else
      {
        (void)modifications_builder.Append();
        for (const auto& mod_pair : mods)
        {
          (void)mod_struct_b->Append();
          (void)mod_position_b->Append(static_cast<int32_t>(mod_pair.first));
          (void)mod_name_b->Append(mod_pair.second.getFullId());
        }
      }

      // === metavalues (not nullable) ===
      (void)metavalues_builder.Append();
      appendMetaValues_(hit, mv_name_b, mv_value_b, mv_type_b, mv_struct_b, excluded_hit_mvs);

    } // end hit loop
  } // end protein identification loop

  // Finalize all arrays
  std::shared_ptr<arrow::Array> arr_accession, arr_score;
  std::shared_ptr<arrow::Array> arr_rank, arr_coverage;
  std::shared_ptr<arrow::Array> arr_sequence, arr_description, arr_is_decoy;
  std::shared_ptr<arrow::Array> arr_run_id;
  std::shared_ptr<arrow::Array> arr_modifications, arr_metavalues;

  status = accession_builder.Finish(&arr_accession);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: accession_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = score_builder.Finish(&arr_score);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: score_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = rank_builder.Finish(&arr_rank);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: rank_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = coverage_builder.Finish(&arr_coverage);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: coverage_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = sequence_builder.Finish(&arr_sequence);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: sequence_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = description_builder.Finish(&arr_description);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: description_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = is_decoy_builder.Finish(&arr_is_decoy);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: is_decoy_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = run_identifier_builder.Finish(&arr_run_id);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: run_identifier_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = modifications_builder.Finish(&arr_modifications);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: modifications_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = metavalues_builder.Finish(&arr_metavalues);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: metavalues_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }

  // Build table using registry schema
  auto schema = ProteinSchema::schema();

  auto table = arrow::Table::Make(schema, {
    arr_accession, arr_score,
    arr_rank, arr_coverage,
    arr_sequence, arr_description, arr_is_decoy,
    arr_run_id,
    arr_modifications, arr_metavalues
  });

  // Validate table against registry schema
  auto validation = ArrowSchemaValidation::validate(table, ProteinSchema::schema());
  if (!validation.valid)
  {
    OPENMS_LOG_ERROR << "Schema validation failed: " << validation.toString() << "\n";
    return nullptr;
  }

  return table;
}


bool ProteinIdentificationArrowIO::exportProteinsToParquet(
  const std::vector<ProteinIdentification>& protein_identifications,
  const String& filename,
  const ParquetWriteConfig& config)
{
  auto table = exportProteinsToArrow(protein_identifications);
  if (!table)
  {
    OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: Failed to create Arrow table for proteins" << std::endl;
    return false;
  }
  return writeArrowTableToParquet_(table, filename, "proteins", config);
}


std::shared_ptr<arrow::Table> ProteinIdentificationArrowIO::exportProteinGroupsToArrow(
  const std::vector<ProteinIdentification>& protein_identifications)
{
  // -- Simple column builders --
  arrow::StringBuilder group_type_builder, run_identifier_builder;
  arrow::DoubleBuilder probability_builder;
  arrow::Int32Builder group_index_builder;

  // -- accessions: list<utf8> --
  auto acc_value_b = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder accessions_builder(arrow::default_memory_pool(), acc_value_b);

  // -- float_data: list<struct{name: utf8, values: list<float64>}> --
  auto fd_inner_value_b = std::make_shared<arrow::DoubleBuilder>();
  auto fd_inner_list_b = std::make_shared<arrow::ListBuilder>(arrow::default_memory_pool(), fd_inner_value_b);
  auto fd_name_b = std::make_shared<arrow::StringBuilder>();
  auto fd_struct_type = arrow::struct_({
    arrow::field("name", arrow::utf8()),
    arrow::field("values", arrow::list(arrow::float64()))
  });
  auto fd_struct_b = std::make_shared<arrow::StructBuilder>(
    fd_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{fd_name_b, fd_inner_list_b});
  arrow::ListBuilder float_data_builder(arrow::default_memory_pool(), fd_struct_b);

  // -- string_data: list<struct{name: utf8, values: list<utf8>}> --
  auto sd_inner_value_b = std::make_shared<arrow::StringBuilder>();
  auto sd_inner_list_b = std::make_shared<arrow::ListBuilder>(arrow::default_memory_pool(), sd_inner_value_b);
  auto sd_name_b = std::make_shared<arrow::StringBuilder>();
  auto sd_struct_type = arrow::struct_({
    arrow::field("name", arrow::utf8()),
    arrow::field("values", arrow::list(arrow::utf8()))
  });
  auto sd_struct_b = std::make_shared<arrow::StructBuilder>(
    sd_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{sd_name_b, sd_inner_list_b});
  arrow::ListBuilder string_data_builder(arrow::default_memory_pool(), sd_struct_b);

  // -- integer_data: list<struct{name: utf8, values: list<int64>}> --
  auto id_inner_value_b = std::make_shared<arrow::Int64Builder>();
  auto id_inner_list_b = std::make_shared<arrow::ListBuilder>(arrow::default_memory_pool(), id_inner_value_b);
  auto id_name_b = std::make_shared<arrow::StringBuilder>();
  auto id_struct_type = arrow::struct_({
    arrow::field("name", arrow::utf8()),
    arrow::field("values", arrow::list(arrow::int64()))
  });
  auto id_struct_b = std::make_shared<arrow::StructBuilder>(
    id_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{id_name_b, id_inner_list_b});
  arrow::ListBuilder integer_data_builder(arrow::default_memory_pool(), id_struct_b);

  // Lambda to append a single ProteinGroup row
  auto append_group = [&](const ProteinIdentification::ProteinGroup& group,
                          const std::string& type_str,
                          const String& run_id,
                          int32_t index) -> bool
  {
    // group_type
    (void)group_type_builder.Append(type_str);

    // probability
    (void)probability_builder.Append(group.probability);

    // accessions (list<utf8>)
    (void)accessions_builder.Append();
    for (const auto& acc : group.accessions)
    {
      (void)acc_value_b->Append(acc);
    }

    // run_identifier
    (void)run_identifier_builder.Append(run_id);

    // group_index
    (void)group_index_builder.Append(index);

    // float_data
    const auto& fda = group.getFloatDataArrays();
    if (fda.empty())
    {
      (void)float_data_builder.AppendNull();
    }
    else
    {
      (void)float_data_builder.Append();
      for (const auto& arr : fda)
      {
        (void)fd_struct_b->Append();
        (void)fd_name_b->Append(arr.getName());
        (void)fd_inner_list_b->Append();
        for (float val : arr)
        {
          (void)fd_inner_value_b->Append(static_cast<double>(val));
        }
      }
    }

    // string_data
    const auto& sda = group.getStringDataArrays();
    if (sda.empty())
    {
      (void)string_data_builder.AppendNull();
    }
    else
    {
      (void)string_data_builder.Append();
      for (const auto& arr : sda)
      {
        (void)sd_struct_b->Append();
        (void)sd_name_b->Append(arr.getName());
        (void)sd_inner_list_b->Append();
        for (const auto& val : arr)
        {
          (void)sd_inner_value_b->Append(val);
        }
      }
    }

    // integer_data
    const auto& ida = group.getIntegerDataArrays();
    if (ida.empty())
    {
      (void)integer_data_builder.AppendNull();
    }
    else
    {
      (void)integer_data_builder.Append();
      for (const auto& arr : ida)
      {
        (void)id_struct_b->Append();
        (void)id_name_b->Append(arr.getName());
        (void)id_inner_list_b->Append();
        for (Int val : arr)
        {
          (void)id_inner_value_b->Append(static_cast<int64_t>(val));
        }
      }
    }

    return true;
  };

  // Iterate over all ProteinIdentifications
  for (const auto& prot_id : protein_identifications)
  {
    const String& run_id = prot_id.getIdentifier();

    // Protein groups with separate 0-based index
    int32_t pg_index = 0;
    for (const auto& group : prot_id.getProteinGroups())
    {
      append_group(group, "protein_group", run_id, pg_index++);
    }

    // Indistinguishable proteins with separate 0-based index
    int32_t ind_index = 0;
    for (const auto& group : prot_id.getIndistinguishableProteins())
    {
      append_group(group, "indistinguishable", run_id, ind_index++);
    }
  }

  // Finalize all arrays
  std::shared_ptr<arrow::Array> arr_group_type, arr_probability, arr_accessions;
  std::shared_ptr<arrow::Array> arr_run_id, arr_group_index;
  std::shared_ptr<arrow::Array> arr_float_data, arr_string_data, arr_integer_data;

  arrow::Status status;
  status = group_type_builder.Finish(&arr_group_type);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: group_type_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = probability_builder.Finish(&arr_probability);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: probability_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = accessions_builder.Finish(&arr_accessions);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: accessions_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = run_identifier_builder.Finish(&arr_run_id);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: run_identifier_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = group_index_builder.Finish(&arr_group_index);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: group_index_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = float_data_builder.Finish(&arr_float_data);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: float_data_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = string_data_builder.Finish(&arr_string_data);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: string_data_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = integer_data_builder.Finish(&arr_integer_data);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: integer_data_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }

  // Build table using registry schema
  auto schema = ProteinGroupSchema::schema();

  auto table = arrow::Table::Make(schema, {
    arr_group_type, arr_probability, arr_accessions,
    arr_run_id, arr_group_index,
    arr_float_data, arr_string_data, arr_integer_data
  });

  // Validate table against registry schema
  auto validation = ArrowSchemaValidation::validate(table, ProteinGroupSchema::schema());
  if (!validation.valid)
  {
    OPENMS_LOG_ERROR << "Schema validation failed: " << validation.toString() << "\n";
    return nullptr;
  }

  return table;
}


bool ProteinIdentificationArrowIO::exportProteinGroupsToParquet(
  const std::vector<ProteinIdentification>& protein_identifications,
  const String& filename,
  const ParquetWriteConfig& config)
{
  auto table = exportProteinGroupsToArrow(protein_identifications);
  if (!table)
  {
    OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: Failed to create Arrow table for protein groups" << std::endl;
    return false;
  }
  return writeArrowTableToParquet_(table, filename, "protein_groups", config);
}


std::shared_ptr<arrow::Table> ProteinIdentificationArrowIO::exportSearchParamsToArrow(
  const std::vector<ProteinIdentification>& protein_identifications)
{
  // -- Simple column builders --
  arrow::StringBuilder run_identifier_builder, search_engine_builder;
  arrow::StringBuilder search_engine_version_builder;
  arrow::StringBuilder inference_engine_builder, inference_engine_version_builder;
  arrow::StringBuilder score_type_builder;
  arrow::TimestampBuilder date_builder(arrow::timestamp(arrow::TimeUnit::SECOND), arrow::default_memory_pool());
  arrow::BooleanBuilder higher_score_better_builder;
  arrow::DoubleBuilder significance_threshold_builder;
  arrow::StringBuilder db_builder, db_version_builder, taxonomy_builder, charges_builder;
  arrow::StringBuilder mass_type_builder;
  arrow::DoubleBuilder precursor_mass_tolerance_builder, fragment_mass_tolerance_builder;
  arrow::BooleanBuilder precursor_mass_tolerance_ppm_builder, fragment_mass_tolerance_ppm_builder;
  arrow::StringBuilder digestion_enzyme_builder, enzyme_term_specificity_builder;
  arrow::Int32Builder missed_cleavages_builder;

  // -- fixed_modifications: list<utf8> --
  auto fixed_mod_value_b = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder fixed_modifications_builder(arrow::default_memory_pool(), fixed_mod_value_b);

  // -- variable_modifications: list<utf8> --
  auto var_mod_value_b = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder variable_modifications_builder(arrow::default_memory_pool(), var_mod_value_b);

  // -- primary_ms_run_paths: list<utf8> --
  auto ms_run_value_b = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder primary_ms_run_paths_builder(arrow::default_memory_pool(), ms_run_value_b);

  // -- metavalues list<struct{name: utf8, value: utf8, value_type: utf8}> (from ProteinIdentification) --
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

  // -- sp_metavalues list<struct{name: utf8, value: utf8, value_type: utf8}> (from SearchParameters) --
  auto sp_mv_name_b = std::make_shared<arrow::StringBuilder>();
  auto sp_mv_value_b = std::make_shared<arrow::StringBuilder>();
  auto sp_mv_type_b = std::make_shared<arrow::StringBuilder>();
  auto sp_mv_struct_b = std::make_shared<arrow::StructBuilder>(
    mv_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{sp_mv_name_b, sp_mv_value_b, sp_mv_type_b});
  arrow::ListBuilder sp_metavalues_builder(arrow::default_memory_pool(), sp_mv_struct_b);

  Size num_rows = protein_identifications.size();

  // Reserve capacity for simple builders
  arrow::Status status;
  status = run_identifier_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: run_identifier_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = search_engine_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: search_engine_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = search_engine_version_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: search_engine_version_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = inference_engine_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: inference_engine_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = inference_engine_version_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: inference_engine_version_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = date_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: date_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = score_type_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: score_type_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = higher_score_better_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: higher_score_better_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = significance_threshold_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: significance_threshold_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = db_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: db_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = db_version_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: db_version_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = taxonomy_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: taxonomy_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = charges_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: charges_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = mass_type_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: mass_type_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = precursor_mass_tolerance_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: precursor_mass_tolerance_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = precursor_mass_tolerance_ppm_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: precursor_mass_tolerance_ppm_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = fragment_mass_tolerance_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: fragment_mass_tolerance_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = fragment_mass_tolerance_ppm_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: fragment_mass_tolerance_ppm_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = digestion_enzyme_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: digestion_enzyme_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = enzyme_term_specificity_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: enzyme_term_specificity_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = missed_cleavages_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: missed_cleavages_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }

  // No metavalue keys to exclude
  static const std::unordered_set<std::string> no_exclusions;

  for (const auto& prot_id : protein_identifications)
  {
    const auto& sp = prot_id.getSearchParameters();

    // === run_identifier (not nullable) ===
    (void)run_identifier_builder.Append(prot_id.getIdentifier());

    // === search_engine (not nullable) ===
    (void)search_engine_builder.Append(prot_id.getSearchEngine());

    // === search_engine_version (nullable, null if empty) ===
    const String& se_version = prot_id.getSearchEngineVersion();
    if (se_version.empty())
    {
      (void)search_engine_version_builder.AppendNull();
    }
    else
    {
      (void)search_engine_version_builder.Append(se_version);
    }

    // === inference_engine (nullable, null if empty) ===
    const String inference_engine = prot_id.getInferenceEngine();
    if (inference_engine.empty())
    {
      (void)inference_engine_builder.AppendNull();
    }
    else
    {
      (void)inference_engine_builder.Append(inference_engine);
    }

    // === inference_engine_version (nullable, null if empty) ===
    const String inference_engine_version = prot_id.getInferenceEngineVersion();
    if (inference_engine_version.empty())
    {
      (void)inference_engine_version_builder.AppendNull();
    }
    else
    {
      (void)inference_engine_version_builder.Append(inference_engine_version);
    }

    // === date (nullable, null if unset) ===
    {
      const DateTime& dt = prot_id.getDateTime();
      String date_str = dt.toString("yyyy-MM-ddThh:mm:ss");
      if (date_str.empty())
      {
        (void)date_builder.AppendNull();
      }
      else
      {
        UInt month, day, year, hour, minute, second;
        dt.get(month, day, year, hour, minute, second);
        std::tm tm{};
        tm.tm_year = static_cast<int>(year) - 1900;
        tm.tm_mon = static_cast<int>(month) - 1;
        tm.tm_mday = static_cast<int>(day);
        tm.tm_hour = static_cast<int>(hour);
        tm.tm_min = static_cast<int>(minute);
        tm.tm_sec = static_cast<int>(second);
#ifdef _WIN32
        int64_t epoch = static_cast<int64_t>(_mkgmtime(&tm));
#else
        int64_t epoch = static_cast<int64_t>(timegm(&tm));
#endif
        if (epoch < 0)
        {
          (void)date_builder.AppendNull(); // invalid date (e.g. default-constructed DateTime)
        }
        else
        {
          (void)date_builder.Append(epoch);
        }
      }
    }

    // === score_type (not nullable) ===
    (void)score_type_builder.Append(prot_id.getScoreType());

    // === higher_score_better (not nullable) ===
    (void)higher_score_better_builder.Append(prot_id.isHigherScoreBetter());

    // === significance_threshold (nullable) ===
    (void)significance_threshold_builder.Append(prot_id.getSignificanceThreshold());

    // === db (nullable, null if empty) ===
    if (sp.db.empty())
    {
      (void)db_builder.AppendNull();
    }
    else
    {
      (void)db_builder.Append(sp.db);
    }

    // === db_version (nullable, null if empty) ===
    if (sp.db_version.empty())
    {
      (void)db_version_builder.AppendNull();
    }
    else
    {
      (void)db_version_builder.Append(sp.db_version);
    }

    // === taxonomy (nullable, null if empty) ===
    if (sp.taxonomy.empty())
    {
      (void)taxonomy_builder.AppendNull();
    }
    else
    {
      (void)taxonomy_builder.Append(sp.taxonomy);
    }

    // === charges (nullable, null if empty) ===
    if (sp.charges.empty())
    {
      (void)charges_builder.AppendNull();
    }
    else
    {
      (void)charges_builder.Append(sp.charges);
    }

    // === mass_type (not nullable) ===
    if (sp.mass_type == ProteinIdentification::PeakMassType::MONOISOTOPIC)
    {
      (void)mass_type_builder.Append("MONOISOTOPIC");
    }
    else
    {
      (void)mass_type_builder.Append("AVERAGE");
    }

    // === precursor_mass_tolerance (not nullable) ===
    (void)precursor_mass_tolerance_builder.Append(sp.precursor_mass_tolerance);

    // === precursor_mass_tolerance_ppm (not nullable) ===
    (void)precursor_mass_tolerance_ppm_builder.Append(sp.precursor_mass_tolerance_ppm);

    // === fragment_mass_tolerance (not nullable) ===
    (void)fragment_mass_tolerance_builder.Append(sp.fragment_mass_tolerance);

    // === fragment_mass_tolerance_ppm (not nullable) ===
    (void)fragment_mass_tolerance_ppm_builder.Append(sp.fragment_mass_tolerance_ppm);

    // === digestion_enzyme (nullable, null if empty or "unknown_enzyme") ===
    const String& enzyme_name = sp.digestion_enzyme.getName();
    if (enzyme_name.empty() || enzyme_name == "unknown_enzyme")
    {
      (void)digestion_enzyme_builder.AppendNull();
    }
    else
    {
      (void)digestion_enzyme_builder.Append(enzyme_name);
    }

    // === enzyme_term_specificity (nullable) ===
    switch (sp.enzyme_term_specificity)
    {
      case EnzymaticDigestion::SPEC_FULL:
        (void)enzyme_term_specificity_builder.Append("FULL");
        break;
      case EnzymaticDigestion::SPEC_SEMI:
        (void)enzyme_term_specificity_builder.Append("SEMI");
        break;
      case EnzymaticDigestion::SPEC_NONE:
        (void)enzyme_term_specificity_builder.Append("NONE");
        break;
      default:
        (void)enzyme_term_specificity_builder.AppendNull();
        break;
    }

    // === missed_cleavages (not nullable) ===
    (void)missed_cleavages_builder.Append(static_cast<int32_t>(sp.missed_cleavages));

    // === fixed_modifications (not nullable, list<utf8>) ===
    (void)fixed_modifications_builder.Append();
    for (const auto& mod : sp.fixed_modifications)
    {
      (void)fixed_mod_value_b->Append(mod);
    }

    // === variable_modifications (not nullable, list<utf8>) ===
    (void)variable_modifications_builder.Append();
    for (const auto& mod : sp.variable_modifications)
    {
      (void)var_mod_value_b->Append(mod);
    }

    // === primary_ms_run_paths (not nullable, list<utf8>) ===
    StringList ms_runs;
    prot_id.getPrimaryMSRunPath(ms_runs);
    (void)primary_ms_run_paths_builder.Append();
    for (const auto& run : ms_runs)
    {
      (void)ms_run_value_b->Append(run);
    }

    // === metavalues (not nullable) - from ProteinIdentification ===
    (void)metavalues_builder.Append();
    appendMetaValues_(prot_id, mv_name_b, mv_value_b, mv_type_b, mv_struct_b, no_exclusions);

    // === sp_metavalues (not nullable) - from SearchParameters ===
    (void)sp_metavalues_builder.Append();
    appendMetaValues_(sp, sp_mv_name_b, sp_mv_value_b, sp_mv_type_b, sp_mv_struct_b, no_exclusions);

  } // end prot_id loop

  // Finalize all arrays
  std::shared_ptr<arrow::Array> arr_run_id, arr_search_engine, arr_se_version;
  std::shared_ptr<arrow::Array> arr_inf_engine, arr_inf_version;
  std::shared_ptr<arrow::Array> arr_date, arr_score_type, arr_higher_better;
  std::shared_ptr<arrow::Array> arr_sig_threshold;
  std::shared_ptr<arrow::Array> arr_db, arr_db_version, arr_taxonomy, arr_charges;
  std::shared_ptr<arrow::Array> arr_mass_type;
  std::shared_ptr<arrow::Array> arr_precursor_tol, arr_precursor_ppm;
  std::shared_ptr<arrow::Array> arr_fragment_tol, arr_fragment_ppm;
  std::shared_ptr<arrow::Array> arr_enzyme, arr_enzyme_spec;
  std::shared_ptr<arrow::Array> arr_missed_cleavages;
  std::shared_ptr<arrow::Array> arr_fixed_mods, arr_var_mods;
  std::shared_ptr<arrow::Array> arr_ms_run_paths;
  std::shared_ptr<arrow::Array> arr_metavalues;
  std::shared_ptr<arrow::Array> arr_sp_metavalues;

  status = run_identifier_builder.Finish(&arr_run_id);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: run_identifier_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = search_engine_builder.Finish(&arr_search_engine);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: search_engine_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = search_engine_version_builder.Finish(&arr_se_version);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: search_engine_version_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = inference_engine_builder.Finish(&arr_inf_engine);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: inference_engine_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = inference_engine_version_builder.Finish(&arr_inf_version);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: inference_engine_version_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = date_builder.Finish(&arr_date);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: date_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = score_type_builder.Finish(&arr_score_type);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: score_type_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = higher_score_better_builder.Finish(&arr_higher_better);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: higher_score_better_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = significance_threshold_builder.Finish(&arr_sig_threshold);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: significance_threshold_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = db_builder.Finish(&arr_db);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: db_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = db_version_builder.Finish(&arr_db_version);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: db_version_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = taxonomy_builder.Finish(&arr_taxonomy);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: taxonomy_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = charges_builder.Finish(&arr_charges);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: charges_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = mass_type_builder.Finish(&arr_mass_type);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: mass_type_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = precursor_mass_tolerance_builder.Finish(&arr_precursor_tol);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: precursor_mass_tolerance_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = precursor_mass_tolerance_ppm_builder.Finish(&arr_precursor_ppm);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: precursor_mass_tolerance_ppm_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = fragment_mass_tolerance_builder.Finish(&arr_fragment_tol);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: fragment_mass_tolerance_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = fragment_mass_tolerance_ppm_builder.Finish(&arr_fragment_ppm);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: fragment_mass_tolerance_ppm_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = digestion_enzyme_builder.Finish(&arr_enzyme);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: digestion_enzyme_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = enzyme_term_specificity_builder.Finish(&arr_enzyme_spec);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: enzyme_term_specificity_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = missed_cleavages_builder.Finish(&arr_missed_cleavages);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: missed_cleavages_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = fixed_modifications_builder.Finish(&arr_fixed_mods);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: fixed_modifications_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = variable_modifications_builder.Finish(&arr_var_mods);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: variable_modifications_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = primary_ms_run_paths_builder.Finish(&arr_ms_run_paths);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: primary_ms_run_paths_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = metavalues_builder.Finish(&arr_metavalues);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: metavalues_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = sp_metavalues_builder.Finish(&arr_sp_metavalues);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: sp_metavalues_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }

  // Build table using registry schema
  auto schema = SearchParamsSchema::schema();

  auto table = arrow::Table::Make(schema, {
    arr_run_id, arr_search_engine, arr_se_version,
    arr_inf_engine, arr_inf_version,
    arr_date, arr_score_type, arr_higher_better,
    arr_sig_threshold,
    arr_db, arr_db_version, arr_taxonomy, arr_charges,
    arr_mass_type,
    arr_precursor_tol, arr_precursor_ppm,
    arr_fragment_tol, arr_fragment_ppm,
    arr_enzyme, arr_enzyme_spec,
    arr_missed_cleavages,
    arr_fixed_mods, arr_var_mods,
    arr_ms_run_paths, arr_metavalues,
    arr_sp_metavalues
  });

  // Validate table against registry schema
  auto validation = ArrowSchemaValidation::validate(table, SearchParamsSchema::schema());
  if (!validation.valid)
  {
    OPENMS_LOG_ERROR << "Schema validation failed: " << validation.toString() << "\n";
    return nullptr;
  }

  return table;
}


bool ProteinIdentificationArrowIO::exportSearchParamsToParquet(
  const std::vector<ProteinIdentification>& protein_identifications,
  const String& filename,
  const ParquetWriteConfig& config)
{
  auto table = exportSearchParamsToArrow(protein_identifications);
  if (!table)
  {
    OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: Failed to create Arrow table for search params" << std::endl;
    return false;
  }
  return writeArrowTableToParquet_(table, filename, "search_params", config);
}


// ==================== Import methods ====================

bool ProteinIdentificationArrowIO::importSearchParamsFromArrow(
  const std::shared_ptr<arrow::Table>& table,
  std::vector<ProteinIdentification>& protein_identifications)
{
  if (!table || table->num_rows() == 0) return true;

  // Validate table schema against registry (subset mode — file may have extra columns)
  auto validation = ArrowSchemaValidation::validate(table, SearchParamsSchema::schema(), ArrowSchemaValidation::Mode::Subset);
  if (!validation.valid)
  {
    OPENMS_LOG_ERROR << "Incompatible schema: " << validation.toString() << "\n";
    return false;
  }

  // Get all columns
  auto col_run_id = getColumn_(table, SearchParamsSchema::RUN_IDENTIFIER);
  auto col_search_engine = getColumn_(table, SearchParamsSchema::SEARCH_ENGINE);
  auto col_se_version = getColumn_(table, SearchParamsSchema::SEARCH_ENGINE_VERSION, false);
  auto col_inf_engine = getColumn_(table, SearchParamsSchema::INFERENCE_ENGINE, false);
  auto col_inf_version = getColumn_(table, SearchParamsSchema::INFERENCE_ENGINE_VERSION, false);
  auto col_date = getColumn_(table, SearchParamsSchema::DATE, false);
  auto col_score_type = getColumn_(table, SearchParamsSchema::SCORE_TYPE);
  auto col_higher_better = getColumn_(table, SearchParamsSchema::HIGHER_SCORE_BETTER);
  auto col_sig_threshold = getColumn_(table, SearchParamsSchema::SIGNIFICANCE_THRESHOLD, false);
  auto col_db = getColumn_(table, SearchParamsSchema::DB, false);
  auto col_db_version = getColumn_(table, SearchParamsSchema::DB_VERSION, false);
  auto col_taxonomy = getColumn_(table, SearchParamsSchema::TAXONOMY, false);
  auto col_charges = getColumn_(table, SearchParamsSchema::CHARGES, false);
  auto col_mass_type = getColumn_(table, SearchParamsSchema::MASS_TYPE);
  auto col_precursor_tol = getColumn_(table, SearchParamsSchema::PRECURSOR_MASS_TOLERANCE);
  auto col_precursor_ppm = getColumn_(table, SearchParamsSchema::PRECURSOR_MASS_TOLERANCE_PPM);
  auto col_fragment_tol = getColumn_(table, SearchParamsSchema::FRAGMENT_MASS_TOLERANCE);
  auto col_fragment_ppm = getColumn_(table, SearchParamsSchema::FRAGMENT_MASS_TOLERANCE_PPM);
  auto col_enzyme = getColumn_(table, SearchParamsSchema::DIGESTION_ENZYME, false);
  auto col_enzyme_spec = getColumn_(table, SearchParamsSchema::ENZYME_TERM_SPECIFICITY, false);
  auto col_missed_cleavages = getColumn_(table, SearchParamsSchema::MISSED_CLEAVAGES);
  auto col_fixed_mods = getColumn_(table, SearchParamsSchema::FIXED_MODIFICATIONS, false);
  auto col_var_mods = getColumn_(table, SearchParamsSchema::VARIABLE_MODIFICATIONS, false);
  auto col_ms_run_paths = getColumn_(table, SearchParamsSchema::PRIMARY_MS_RUN_PATHS, false);
  auto col_metavalues = getColumn_(table, SearchParamsSchema::METAVALUES, false);
  auto col_sp_metavalues = getColumn_(table, SearchParamsSchema::SP_METAVALUES, false);

  if (!col_run_id || !col_search_engine || !col_score_type || !col_higher_better ||
      !col_mass_type || !col_precursor_tol || !col_precursor_ppm ||
      !col_fragment_tol || !col_fragment_ppm || !col_missed_cleavages)
  {
    OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: Missing required columns in search_params table" << std::endl;
    return false;
  }

  for (int64_t row = 0; row < table->num_rows(); ++row)
  {
    ProteinIdentification prot_id;

    // Run-level metadata
    prot_id.setIdentifier(getStringValue_(col_run_id, row));
    prot_id.setSearchEngine(getStringValue_(col_search_engine, row));
    prot_id.setSearchEngineVersion(getStringValue_(col_se_version, row));
    prot_id.setScoreType(getStringValue_(col_score_type, row));
    prot_id.setHigherScoreBetter(getBoolValue_(col_higher_better, row));
    prot_id.setSignificanceThreshold(getDoubleValue_(col_sig_threshold, row));

    // Date (timestamp; normalize to seconds based on the column's actual time unit).
    // The schema declares timestamp(SECOND), but Parquet may persist/restore the column
    // as nanoseconds, so we have to inspect the runtime unit rather than assume seconds.
    if (!isNull_(col_date, row))
    {
      auto ts_arr = std::static_pointer_cast<arrow::TimestampArray>(col_date);
      auto ts_type = std::static_pointer_cast<arrow::TimestampType>(ts_arr->type());
      int64_t raw = ts_arr->Value(row);
      // Floor-style division so a small negative sub-second value (e.g. raw = -500
      // ms) maps to -1 second instead of 0; with C++'s truncating integer division
      // it would otherwise pass the epoch_secs >= 0 check below and silently render
      // as 1970-01-01.  Quotient/remainder form avoids UB at INT64_MIN that the
      // simpler `-((-r + d - 1) / d)` formulation would hit by negating r.
      auto floor_div = [](int64_t r, int64_t d) -> int64_t {
        int64_t q = r / d;
        int64_t rem = r % d;
        if (rem != 0 && ((rem < 0) != (d < 0))) { --q; }
        return q;
      };
      int64_t epoch_secs = raw;
      switch (ts_type->unit())
      {
        case arrow::TimeUnit::SECOND: epoch_secs = raw; break;
        case arrow::TimeUnit::MILLI:  epoch_secs = floor_div(raw, 1000LL); break;
        case arrow::TimeUnit::MICRO:  epoch_secs = floor_div(raw, 1000000LL); break;
        case arrow::TimeUnit::NANO:   epoch_secs = floor_div(raw, 1000000000LL); break;
      }
      if (epoch_secs >= 0) // negative epochs are invalid on Windows
      {
        time_t t = static_cast<time_t>(epoch_secs);
        std::tm tm{};
#ifdef _WIN32
        errno_t err = gmtime_s(&tm, &t);
        if (err == 0)
#else
        if (gmtime_r(&t, &tm) != nullptr)
#endif
        {
          DateTime dt;
          dt.set(static_cast<UInt>(tm.tm_mon + 1),
                 static_cast<UInt>(tm.tm_mday),
                 static_cast<UInt>(tm.tm_year + 1900),
                 static_cast<UInt>(tm.tm_hour),
                 static_cast<UInt>(tm.tm_min),
                 static_cast<UInt>(tm.tm_sec));
          prot_id.setDateTime(dt);
        }
      }
    }

    // Primary MS run paths — only call when the column is present.
    // Missing column → no spectra_data UserParam (preserves round-trip semantics).
    // Present-but-empty column → empty spectra_data UserParam (matches idXML behaviour).
    if (col_ms_run_paths)
    {
      auto ms_run_paths = readStringList_(col_ms_run_paths, row);
      StringList sl(ms_run_paths.begin(), ms_run_paths.end());
      prot_id.setPrimaryMSRunPath(sl);
    }

    // SearchParameters
    ProteinIdentification::SearchParameters sp;
    sp.db = getStringValue_(col_db, row);
    sp.db_version = getStringValue_(col_db_version, row);
    sp.taxonomy = getStringValue_(col_taxonomy, row);
    sp.charges = getStringValue_(col_charges, row);

    String mass_type_str = getStringValue_(col_mass_type, row);
    sp.mass_type = (mass_type_str == "AVERAGE") ?
      ProteinIdentification::PeakMassType::AVERAGE :
      ProteinIdentification::PeakMassType::MONOISOTOPIC;

    sp.precursor_mass_tolerance = getDoubleValue_(col_precursor_tol, row);
    sp.precursor_mass_tolerance_ppm = getBoolValue_(col_precursor_ppm, row);
    sp.fragment_mass_tolerance = getDoubleValue_(col_fragment_tol, row);
    sp.fragment_mass_tolerance_ppm = getBoolValue_(col_fragment_ppm, row);
    sp.missed_cleavages = static_cast<UInt>(getInt32Value_(col_missed_cleavages, row));

    // Enzyme
    String enzyme_name = getStringValue_(col_enzyme, row);
    if (!enzyme_name.empty())
    {
      if (ProteaseDB::getInstance()->hasEnzyme(enzyme_name))
      {
        sp.digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme(enzyme_name);
      }
      else
      {
        OPENMS_LOG_WARN << "ProteinIdentificationArrowIO: Unknown enzyme '" << enzyme_name << "'" << std::endl;
      }
    }

    // Enzyme term specificity
    String spec_str = getStringValue_(col_enzyme_spec, row);
    if (spec_str == "FULL") sp.enzyme_term_specificity = EnzymaticDigestion::SPEC_FULL;
    else if (spec_str == "SEMI") sp.enzyme_term_specificity = EnzymaticDigestion::SPEC_SEMI;
    else if (spec_str == "NONE") sp.enzyme_term_specificity = EnzymaticDigestion::SPEC_NONE;
    else if (!spec_str.empty())
    {
      OPENMS_LOG_WARN << "ProteinIdentificationArrowIO: Unknown enzyme_term_specificity '"
                      << spec_str << "'" << std::endl;
    }

    // Modifications
    auto fixed_mods = readStringList_(col_fixed_mods, row);
    sp.fixed_modifications.assign(fixed_mods.begin(), fixed_mods.end());
    auto var_mods = readStringList_(col_var_mods, row);
    sp.variable_modifications.assign(var_mods.begin(), var_mods.end());

    // SearchParameters metavalues (restore before setSearchParameters)
    ArrowIOHelpers::readMetaValues(col_sp_metavalues, row, sp);

    prot_id.setSearchParameters(sp);

    // Inference engine (must be set after setSearchParameters)
    String inf_engine = getStringValue_(col_inf_engine, row);
    if (!inf_engine.empty())
    {
      prot_id.setInferenceEngine(inf_engine);
      String inf_version = getStringValue_(col_inf_version, row);
      if (!inf_version.empty())
      {
        prot_id.setInferenceEngineVersion(inf_version);
      }
    }

    // ProteinIdentification metavalues (exclude keys handled via dedicated columns)
    static const std::unordered_set<std::string> excluded_prot_id_mvs = {
      "InferenceEngine", "InferenceEngineVersion",
      "spectra_data", "spectra_data_raw"
    };
    ArrowIOHelpers::readMetaValues(col_metavalues, row, prot_id, excluded_prot_id_mvs);

    protein_identifications.push_back(std::move(prot_id));
  }

  return true;
}


bool ProteinIdentificationArrowIO::importProteinsFromArrow(
  const std::shared_ptr<arrow::Table>& table,
  std::vector<ProteinIdentification>& protein_identifications)
{
  if (!table || table->num_rows() == 0) return true;

  // Validate table schema against registry (subset mode — file may have extra columns)
  auto validation = ArrowSchemaValidation::validate(table, ProteinSchema::schema(), ArrowSchemaValidation::Mode::Subset);
  if (!validation.valid)
  {
    OPENMS_LOG_ERROR << "Incompatible schema: " << validation.toString() << "\n";
    return false;
  }

  // Get all columns
  auto col_accession = getColumn_(table, ProteinSchema::ACCESSION);
  auto col_score = getColumn_(table, ProteinSchema::SCORE);
  auto col_rank = getColumn_(table, ProteinSchema::RANK, false);
  auto col_coverage = getColumn_(table, ProteinSchema::COVERAGE, false);
  auto col_sequence = getColumn_(table, ProteinSchema::SEQUENCE, false);
  auto col_description = getColumn_(table, ProteinSchema::DESCRIPTION, false);
  auto col_is_decoy = getColumn_(table, ProteinSchema::IS_DECOY, false);
  auto col_run_id = getColumn_(table, ProteinSchema::RUN_IDENTIFIER);
  auto col_modifications = getColumn_(table, ProteinSchema::MODIFICATIONS, false);
  auto col_metavalues = getColumn_(table, ProteinSchema::METAVALUES, false);

  if (!col_accession || !col_score || !col_run_id)
  {
    OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: Missing required columns in proteins table" << std::endl;
    return false;
  }

  // Build map from run_identifier to index
  auto run_id_map = buildRunIdMap_(protein_identifications);

  // Metavalue keys excluded from hit metavalues (they have dedicated columns)
  static const std::unordered_set<std::string> excluded_hit_mvs = {
    "Description", "target_decoy"
  };

  for (int64_t row = 0; row < table->num_rows(); ++row)
  {
    String run_id = getStringValue_(col_run_id, row);

    // Find or create matching ProteinIdentification
    auto it = run_id_map.find(run_id);
    size_t prot_id_idx;
    if (it != run_id_map.end())
    {
      prot_id_idx = it->second;
    }
    else
    {
      // Create minimal ProteinIdentification with just run_identifier
      OPENMS_LOG_WARN << "ProteinIdentificationArrowIO: run_identifier '" << run_id
                      << "' not found in search_params; creating minimal ProteinIdentification" << std::endl;
      ProteinIdentification new_prot_id;
      new_prot_id.setIdentifier(run_id);

      prot_id_idx = protein_identifications.size();
      protein_identifications.push_back(std::move(new_prot_id));
      run_id_map[run_id] = prot_id_idx;
    }

    // Build ProteinHit
    ProteinHit hit;
    hit.setAccession(getStringValue_(col_accession, row));
    hit.setScore(getDoubleValue_(col_score, row));

    if (!isNull_(col_rank, row))
    {
      hit.setRank(static_cast<UInt>(getInt32Value_(col_rank, row)));
    }

    if (!isNull_(col_coverage, row))
    {
      hit.setCoverage(getDoubleValue_(col_coverage, row));
    }

    if (!isNull_(col_sequence, row))
    {
      hit.setSequence(getStringValue_(col_sequence, row));
    }

    if (!isNull_(col_description, row))
    {
      hit.setDescription(getStringValue_(col_description, row));
    }

    // is_decoy -> target_decoy metavalue
    if (!isNull_(col_is_decoy, row))
    {
      bool is_decoy = getBoolValue_(col_is_decoy, row);
      hit.setMetaValue("target_decoy", is_decoy ? "decoy" : "target");
    }

    // Modifications
    if (col_modifications && !col_modifications->IsNull(row))
    {
      auto list_arr = std::static_pointer_cast<arrow::ListArray>(col_modifications);
      auto struct_arr = std::static_pointer_cast<arrow::StructArray>(list_arr->value_slice(row));
      if (struct_arr && struct_arr->length() > 0)
      {
        auto pos_arr = std::static_pointer_cast<arrow::Int32Array>(struct_arr->field(0));
        auto mod_name_arr = std::static_pointer_cast<arrow::StringArray>(struct_arr->field(1));

        std::set<std::pair<Size, ResidueModification>> mods;
        for (int64_t i = 0; i < struct_arr->length(); ++i)
        {
          Size position = static_cast<Size>(pos_arr->Value(i));
          String mod_id = mod_name_arr->GetString(i);
          try
          {
            const ResidueModification* mod = ModificationsDB::getInstance()->getModification(mod_id);
            mods.insert(std::make_pair(position, *mod));
          }
          catch (...)
          {
            OPENMS_LOG_WARN << "ProteinIdentificationArrowIO: Could not find modification '" << mod_id << "'" << std::endl;
          }
        }
        hit.setModifications(mods);
      }
    }

    // MetaValues
    ArrowIOHelpers::readMetaValues(col_metavalues, row, hit, excluded_hit_mvs);

    protein_identifications[prot_id_idx].insertHit(std::move(hit));
  }

  return true;
}


bool ProteinIdentificationArrowIO::importProteinGroupsFromArrow(
  const std::shared_ptr<arrow::Table>& table,
  std::vector<ProteinIdentification>& protein_identifications)
{
  if (!table || table->num_rows() == 0) return true;

  // Validate table schema against registry (subset mode — file may have extra columns)
  auto validation = ArrowSchemaValidation::validate(table, ProteinGroupSchema::schema(), ArrowSchemaValidation::Mode::Subset);
  if (!validation.valid)
  {
    OPENMS_LOG_ERROR << "Incompatible schema: " << validation.toString() << "\n";
    return false;
  }

  // Get all columns
  auto col_group_type = getColumn_(table, ProteinGroupSchema::GROUP_TYPE);
  auto col_probability = getColumn_(table, ProteinGroupSchema::PROBABILITY);
  auto col_accessions = getColumn_(table, ProteinGroupSchema::ACCESSIONS);
  auto col_run_id = getColumn_(table, ProteinGroupSchema::RUN_IDENTIFIER);
  auto col_group_index = getColumn_(table, ProteinGroupSchema::GROUP_INDEX, false); // informational, not needed for reconstruction
  auto col_float_data = getColumn_(table, ProteinGroupSchema::FLOAT_DATA, false);
  auto col_string_data = getColumn_(table, ProteinGroupSchema::STRING_DATA, false);
  auto col_integer_data = getColumn_(table, ProteinGroupSchema::INTEGER_DATA, false);

  if (!col_group_type || !col_probability || !col_accessions || !col_run_id)
  {
    OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: Missing required columns in protein_groups table" << std::endl;
    return false;
  }

  auto run_id_map = buildRunIdMap_(protein_identifications);

  // Helper lambda to read data arrays from list<struct{name,values: list<T>}>
  auto readFloatDataArrays = [](const std::shared_ptr<arrow::Array>& array, int64_t row,
                                ProteinIdentification::ProteinGroup& group)
  {
    if (!array || array->IsNull(row)) return;
    auto list_arr = std::static_pointer_cast<arrow::ListArray>(array);
    auto struct_arr = std::static_pointer_cast<arrow::StructArray>(list_arr->value_slice(row));
    if (!struct_arr || struct_arr->length() == 0) return;

    auto name_arr = std::static_pointer_cast<arrow::StringArray>(struct_arr->field(0));
    auto values_list_arr = std::static_pointer_cast<arrow::ListArray>(struct_arr->field(1));

    for (int64_t i = 0; i < struct_arr->length(); ++i)
    {
      ProteinIdentification::ProteinGroup::FloatDataArray fda;
      fda.setName(name_arr->GetString(i));
      auto vals = std::static_pointer_cast<arrow::DoubleArray>(values_list_arr->value_slice(i));
      for (int64_t j = 0; j < vals->length(); ++j)
      {
        fda.push_back(static_cast<float>(vals->Value(j)));
      }
      group.getFloatDataArrays().push_back(std::move(fda));
    }
  };

  auto readStringDataArrays = [](const std::shared_ptr<arrow::Array>& array, int64_t row,
                                 ProteinIdentification::ProteinGroup& group)
  {
    if (!array || array->IsNull(row)) return;
    auto list_arr = std::static_pointer_cast<arrow::ListArray>(array);
    auto struct_arr = std::static_pointer_cast<arrow::StructArray>(list_arr->value_slice(row));
    if (!struct_arr || struct_arr->length() == 0) return;

    auto name_arr = std::static_pointer_cast<arrow::StringArray>(struct_arr->field(0));
    auto values_list_arr = std::static_pointer_cast<arrow::ListArray>(struct_arr->field(1));

    for (int64_t i = 0; i < struct_arr->length(); ++i)
    {
      ProteinIdentification::ProteinGroup::StringDataArray sda;
      sda.setName(name_arr->GetString(i));
      auto vals = std::static_pointer_cast<arrow::StringArray>(values_list_arr->value_slice(i));
      for (int64_t j = 0; j < vals->length(); ++j)
      {
        sda.push_back(vals->GetString(j));
      }
      group.getStringDataArrays().push_back(std::move(sda));
    }
  };

  auto readIntegerDataArrays = [](const std::shared_ptr<arrow::Array>& array, int64_t row,
                                  ProteinIdentification::ProteinGroup& group)
  {
    if (!array || array->IsNull(row)) return;
    auto list_arr = std::static_pointer_cast<arrow::ListArray>(array);
    auto struct_arr = std::static_pointer_cast<arrow::StructArray>(list_arr->value_slice(row));
    if (!struct_arr || struct_arr->length() == 0) return;

    auto name_arr = std::static_pointer_cast<arrow::StringArray>(struct_arr->field(0));
    auto values_list_arr = std::static_pointer_cast<arrow::ListArray>(struct_arr->field(1));

    for (int64_t i = 0; i < struct_arr->length(); ++i)
    {
      ProteinIdentification::ProteinGroup::IntegerDataArray ida;
      ida.setName(name_arr->GetString(i));
      auto vals = std::static_pointer_cast<arrow::Int64Array>(values_list_arr->value_slice(i));
      for (int64_t j = 0; j < vals->length(); ++j)
      {
        ida.push_back(static_cast<Int>(vals->Value(j)));
      }
      group.getIntegerDataArrays().push_back(std::move(ida));
    }
  };

  for (int64_t row = 0; row < table->num_rows(); ++row)
  {
    String run_id = getStringValue_(col_run_id, row);

    // Find or create matching ProteinIdentification
    auto it = run_id_map.find(run_id);
    size_t prot_id_idx;
    if (it != run_id_map.end())
    {
      prot_id_idx = it->second;
    }
    else
    {
      ProteinIdentification new_prot_id;
      new_prot_id.setIdentifier(run_id);
      prot_id_idx = protein_identifications.size();
      protein_identifications.push_back(std::move(new_prot_id));
      run_id_map[run_id] = prot_id_idx;
    }

    // Build ProteinGroup
    ProteinIdentification::ProteinGroup group;
    group.probability = getDoubleValue_(col_probability, row);
    auto accessions = readStringList_(col_accessions, row);
    group.accessions.assign(accessions.begin(), accessions.end());

    // Data arrays
    readFloatDataArrays(col_float_data, row, group);
    readStringDataArrays(col_string_data, row, group);
    readIntegerDataArrays(col_integer_data, row, group);

    // Insert into appropriate vector
    String group_type = getStringValue_(col_group_type, row);
    if (group_type == "indistinguishable")
    {
      protein_identifications[prot_id_idx].insertIndistinguishableProteins(std::move(group));
    }
    else
    {
      protein_identifications[prot_id_idx].insertProteinGroup(std::move(group));
    }
  }

  return true;
}


bool ProteinIdentificationArrowIO::importSearchParamsFromParquet(
  const String& filename,
  std::vector<ProteinIdentification>& protein_identifications)
{
  auto table = readParquetTable_(filename);
  if (!table) return false;
  return importSearchParamsFromArrow(table, protein_identifications);
}


bool ProteinIdentificationArrowIO::importProteinsFromParquet(
  const String& filename,
  std::vector<ProteinIdentification>& protein_identifications)
{
  auto table = readParquetTable_(filename);
  if (!table) return false;
  return importProteinsFromArrow(table, protein_identifications);
}


bool ProteinIdentificationArrowIO::importProteinGroupsFromParquet(
  const String& filename,
  std::vector<ProteinIdentification>& protein_identifications)
{
  auto table = readParquetTable_(filename);
  if (!table) return false;
  return importProteinGroupsFromArrow(table, protein_identifications);
}


bool ProteinIdentificationArrowIO::importFromParquet(
  const String& proteins_filename,
  const String& protein_groups_filename,
  const String& search_params_filename,
  std::vector<ProteinIdentification>& protein_identifications)
{
  protein_identifications.clear();

  // 1. Import search params first (creates ProteinIdentification shells)
  if (!importSearchParamsFromParquet(search_params_filename, protein_identifications))
  {
    OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: Failed to import search params" << std::endl;
    return false;
  }

  // 2. Import proteins second (adds ProteinHits to matching shells)
  if (!importProteinsFromParquet(proteins_filename, protein_identifications))
  {
    OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: Failed to import proteins" << std::endl;
    return false;
  }

  // 3. Import protein groups third (adds groups to matching shells)
  if (!importProteinGroupsFromParquet(protein_groups_filename, protein_identifications))
  {
    OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: Failed to import protein groups" << std::endl;
    return false;
  }

  return true;
}


// ==================== Identifier handling parity with XML lane ====================

std::map<String, String> ProteinIdentificationArrowIO::synthesizeRunIdentifiers(
  std::vector<ProteinIdentification>& protein_identifications)
{
  std::map<String, String> rename;
  for (auto& prot_id : protein_identifications)
  {
    // Mirror IdXMLFile.cpp:530 — `<search_engine>_<date>_<UniqueIdGenerator>`.
    const String se = prot_id.getSearchEngine().empty() ? String("unknown") : prot_id.getSearchEngine();
    const String dt = prot_id.getDateTime().isValid()
                          ? prot_id.getDateTime().toString()
                          : String("1900-01-01T00:00:00");
    const String stored = prot_id.getIdentifier();
    const String synthesized = se + "_" + dt + "_" + String(UniqueIdGenerator::getUniqueId());

    // Multiple ProtIDs sharing the stored identifier each get their own distinct
    // synthesized identifier; the rename map collapses to the last-seen entry
    // (matches XML-lane behavior should pre-fix files reach the loader). One
    // ProtID ends up orphaned of its pep_ids; warn once per collision so the
    // upstream data corruption is visible.
    auto existing = rename.find(stored);
    if (existing != rename.end())
    {
      OPENMS_LOG_WARN << "ProteinIdentificationArrowIO: multiple ProtIDs share stored identifier '"
                      << stored << "'; pep_id assignment ambiguity resolved by load order — "
                      << "regenerate input from a fixed-code build to eliminate." << std::endl;
      existing->second = synthesized;
    }
    else
    {
      rename[stored] = synthesized;
    }

    prot_id.setIdentifier(synthesized);
  }
  return rename;
}

void ProteinIdentificationArrowIO::applyRunIdentifierRename(
  const std::map<String, String>& rename,
  PeptideIdentificationList& pep_ids)
{
  for (auto& pid : pep_ids)
  {
    auto it = rename.find(pid.getIdentifier());
    if (it != rename.end())
    {
      pid.setIdentifier(it->second);
    }
    // Orphan pep_ids (no matching ProtID identifier) keep their stored identifier —
    // same as the XML lane, where they organically surface upstream-corruption signals.
  }
}

void ProteinIdentificationArrowIO::checkUniqueIdentifiers(
  const std::vector<ProteinIdentification>& protein_identifications)
{
  // Mirror XMLHandler::checkUniqueIdentifiers_ exactly — identical message text
  // so log-grepping finds both lanes.
  std::set<String> s;
  for (const auto& p : protein_identifications)
  {
    if (s.insert(p.getIdentifier()).second == false)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "ProteinIdentification run identifiers are not unique. This can lead to "
        "loss of unique PeptideIdentification assignment. Duplicated Protein-ID is:",
        p.getIdentifier());
    }
  }
}


} // namespace OpenMS
