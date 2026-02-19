// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ProteinIdentificationArrowExport.h>

#ifdef WITH_PARQUET

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>

#include <arrow/api.h>
#include <arrow/builder.h>
#include <arrow/io/file.h>
#include <parquet/arrow/writer.h>
#include <parquet/properties.h>

#include <cstdio>
#include <cstring>
#include <random>
#include <unordered_set>
#include <vector>

namespace OpenMS
{

namespace // anonymous
{
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
        case DataValue::DOUBLE_VALUE: (void)type_b->Append("float"); break;
        case DataValue::STRING_VALUE: (void)type_b->Append("str"); break;
        default: (void)type_b->Append("str"); break;
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
      OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: null table passed to writeArrowTableToParquet_" << std::endl;
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
      OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: Failed to open file: " << filename << std::endl;
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
      OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: Failed to write Parquet: "
                       << write_status.ToString() << std::endl;
      return false;
    }

    return true;
  }

} // anonymous namespace


std::shared_ptr<arrow::Table> ProteinIdentificationArrowExport::exportProteinsToArrow(
  const std::vector<ProteinIdentification>& protein_identifications)
{
  // -- Simple column builders --
  arrow::StringBuilder accession_builder, score_type_builder, sequence_builder;
  arrow::StringBuilder description_builder, run_identifier_builder;
  arrow::StringBuilder reference_file_builder, search_engine_builder;
  arrow::StringBuilder search_engine_version_builder;
  arrow::StringBuilder inference_engine_builder, inference_engine_version_builder;
  arrow::StringBuilder date_builder;
  arrow::DoubleBuilder score_builder, coverage_builder, significance_threshold_builder;
  arrow::BooleanBuilder higher_score_better_builder;
  arrow::Int32Builder rank_builder, is_decoy_builder;

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
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: accession_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = score_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: score_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = score_type_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: score_type_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = higher_score_better_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: higher_score_better_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = rank_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: rank_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = coverage_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: coverage_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = sequence_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: sequence_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = description_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: description_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = is_decoy_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: is_decoy_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = run_identifier_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: run_identifier_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = reference_file_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: reference_file_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = search_engine_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: search_engine_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = search_engine_version_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: search_engine_version_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = inference_engine_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: inference_engine_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = inference_engine_version_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: inference_engine_version_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = significance_threshold_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: significance_threshold_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = date_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: date_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }

  // Metavalue keys excluded from metavalues column (they have dedicated columns)
  static const std::unordered_set<std::string> excluded_hit_mvs = {
    "Description", "target_decoy"
  };

  for (const auto& prot_id : protein_identifications)
  {
    // Get shared per-run values
    const String& score_type = prot_id.getScoreType();
    bool higher_better = prot_id.isHigherScoreBetter();
    const String& run_id = prot_id.getIdentifier();
    const String& search_engine = prot_id.getSearchEngine();
    const String& search_engine_version = prot_id.getSearchEngineVersion();
    const String inference_engine = prot_id.getInferenceEngine();
    const String inference_engine_version = prot_id.getInferenceEngineVersion();
    double sig_threshold = prot_id.getSignificanceThreshold();
    String date_str = prot_id.getDateTime().toString("yyyy-MM-ddThh:mm:ss");

    // Get primary MS run path
    StringList ms_runs;
    prot_id.getPrimaryMSRunPath(ms_runs);
    String ref_file = ms_runs.empty() ? "" : ms_runs[0];

    for (const auto& hit : prot_id.getHits())
    {
      // === accession (not nullable) ===
      (void)accession_builder.Append(hit.getAccession());

      // === score (not nullable) ===
      (void)score_builder.Append(hit.getScore());

      // === score_type (not nullable) ===
      (void)score_type_builder.Append(score_type);

      // === higher_score_better (not nullable) ===
      (void)higher_score_better_builder.Append(higher_better);

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
          (void)is_decoy_builder.Append(td_type == ProteinHit::TargetDecoyType::DECOY ? 1 : 0);
        }
      }

      // === run_identifier (not nullable) ===
      (void)run_identifier_builder.Append(run_id);

      // === reference_file_name (nullable) ===
      if (ref_file.empty())
      {
        (void)reference_file_builder.AppendNull();
      }
      else
      {
        (void)reference_file_builder.Append(ref_file);
      }

      // === search_engine (not nullable) ===
      (void)search_engine_builder.Append(search_engine);

      // === search_engine_version (nullable, null if empty) ===
      if (search_engine_version.empty())
      {
        (void)search_engine_version_builder.AppendNull();
      }
      else
      {
        (void)search_engine_version_builder.Append(search_engine_version);
      }

      // === inference_engine (nullable, null if empty) ===
      if (inference_engine.empty())
      {
        (void)inference_engine_builder.AppendNull();
      }
      else
      {
        (void)inference_engine_builder.Append(inference_engine);
      }

      // === inference_engine_version (nullable, null if empty) ===
      if (inference_engine_version.empty())
      {
        (void)inference_engine_version_builder.AppendNull();
      }
      else
      {
        (void)inference_engine_version_builder.Append(inference_engine_version);
      }

      // === significance_threshold (nullable) ===
      (void)significance_threshold_builder.Append(sig_threshold);

      // === date (nullable, null if empty) ===
      if (date_str.empty())
      {
        (void)date_builder.AppendNull();
      }
      else
      {
        (void)date_builder.Append(date_str);
      }

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
  std::shared_ptr<arrow::Array> arr_accession, arr_score, arr_score_type;
  std::shared_ptr<arrow::Array> arr_higher_better, arr_rank, arr_coverage;
  std::shared_ptr<arrow::Array> arr_sequence, arr_description, arr_is_decoy;
  std::shared_ptr<arrow::Array> arr_run_id, arr_ref_file, arr_search_engine;
  std::shared_ptr<arrow::Array> arr_se_version, arr_inf_engine, arr_inf_version;
  std::shared_ptr<arrow::Array> arr_sig_threshold, arr_date;
  std::shared_ptr<arrow::Array> arr_modifications, arr_metavalues;

  status = accession_builder.Finish(&arr_accession);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: accession_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = score_builder.Finish(&arr_score);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: score_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = score_type_builder.Finish(&arr_score_type);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: score_type_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = higher_score_better_builder.Finish(&arr_higher_better);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: higher_score_better_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = rank_builder.Finish(&arr_rank);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: rank_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = coverage_builder.Finish(&arr_coverage);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: coverage_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = sequence_builder.Finish(&arr_sequence);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: sequence_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = description_builder.Finish(&arr_description);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: description_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = is_decoy_builder.Finish(&arr_is_decoy);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: is_decoy_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = run_identifier_builder.Finish(&arr_run_id);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: run_identifier_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = reference_file_builder.Finish(&arr_ref_file);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: reference_file_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = search_engine_builder.Finish(&arr_search_engine);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: search_engine_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = search_engine_version_builder.Finish(&arr_se_version);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: search_engine_version_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = inference_engine_builder.Finish(&arr_inf_engine);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: inference_engine_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = inference_engine_version_builder.Finish(&arr_inf_version);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: inference_engine_version_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = significance_threshold_builder.Finish(&arr_sig_threshold);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: significance_threshold_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = date_builder.Finish(&arr_date);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: date_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = modifications_builder.Finish(&arr_modifications);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: modifications_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = metavalues_builder.Finish(&arr_metavalues);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: metavalues_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }

  // Build schema (19 columns)
  auto schema = arrow::schema({
    arrow::field("accession", arrow::utf8(), /*nullable=*/false),
    arrow::field("score", arrow::float64(), /*nullable=*/false),
    arrow::field("score_type", arrow::utf8(), /*nullable=*/false),
    arrow::field("higher_score_better", arrow::boolean(), /*nullable=*/false),
    arrow::field("rank", arrow::int32(), /*nullable=*/true),
    arrow::field("coverage", arrow::float64(), /*nullable=*/true),
    arrow::field("sequence", arrow::utf8(), /*nullable=*/true),
    arrow::field("description", arrow::utf8(), /*nullable=*/true),
    arrow::field("is_decoy", arrow::int32(), /*nullable=*/true),
    arrow::field("run_identifier", arrow::utf8(), /*nullable=*/false),
    arrow::field("reference_file_name", arrow::utf8(), /*nullable=*/true),
    arrow::field("search_engine", arrow::utf8(), /*nullable=*/false),
    arrow::field("search_engine_version", arrow::utf8(), /*nullable=*/true),
    arrow::field("inference_engine", arrow::utf8(), /*nullable=*/true),
    arrow::field("inference_engine_version", arrow::utf8(), /*nullable=*/true),
    arrow::field("significance_threshold", arrow::float64(), /*nullable=*/true),
    arrow::field("date", arrow::utf8(), /*nullable=*/true),
    arrow::field("modifications", arrow::list(mod_struct_type), /*nullable=*/true),
    arrow::field("metavalues", metavalues_builder.type(), /*nullable=*/false),
  });

  auto table = arrow::Table::Make(schema, {
    arr_accession, arr_score, arr_score_type,
    arr_higher_better, arr_rank, arr_coverage,
    arr_sequence, arr_description, arr_is_decoy,
    arr_run_id, arr_ref_file, arr_search_engine,
    arr_se_version, arr_inf_engine, arr_inf_version,
    arr_sig_threshold, arr_date,
    arr_modifications, arr_metavalues
  });

  return table;
}


bool ProteinIdentificationArrowExport::exportProteinsToParquet(
  const std::vector<ProteinIdentification>& protein_identifications,
  const String& filename,
  const ParquetWriteConfig& config)
{
  auto table = exportProteinsToArrow(protein_identifications);
  if (!table)
  {
    OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: Failed to create Arrow table for proteins" << std::endl;
    return false;
  }
  return writeArrowTableToParquet_(table, filename, "proteins", config);
}


std::shared_ptr<arrow::Table> ProteinIdentificationArrowExport::exportProteinGroupsToArrow(
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
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: group_type_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = probability_builder.Finish(&arr_probability);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: probability_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = accessions_builder.Finish(&arr_accessions);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: accessions_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = run_identifier_builder.Finish(&arr_run_id);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: run_identifier_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = group_index_builder.Finish(&arr_group_index);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: group_index_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = float_data_builder.Finish(&arr_float_data);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: float_data_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = string_data_builder.Finish(&arr_string_data);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: string_data_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = integer_data_builder.Finish(&arr_integer_data);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: integer_data_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }

  // Build schema (8 columns)
  auto schema = arrow::schema({
    arrow::field("group_type", arrow::utf8(), /*nullable=*/false),
    arrow::field("probability", arrow::float64(), /*nullable=*/false),
    arrow::field("accessions", arrow::list(arrow::utf8()), /*nullable=*/false),
    arrow::field("run_identifier", arrow::utf8(), /*nullable=*/false),
    arrow::field("group_index", arrow::int32(), /*nullable=*/false),
    arrow::field("float_data", arrow::list(fd_struct_type), /*nullable=*/true),
    arrow::field("string_data", arrow::list(sd_struct_type), /*nullable=*/true),
    arrow::field("integer_data", arrow::list(id_struct_type), /*nullable=*/true),
  });

  auto table = arrow::Table::Make(schema, {
    arr_group_type, arr_probability, arr_accessions,
    arr_run_id, arr_group_index,
    arr_float_data, arr_string_data, arr_integer_data
  });

  return table;
}


bool ProteinIdentificationArrowExport::exportProteinGroupsToParquet(
  const std::vector<ProteinIdentification>& protein_identifications,
  const String& filename,
  const ParquetWriteConfig& config)
{
  auto table = exportProteinGroupsToArrow(protein_identifications);
  if (!table)
  {
    OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: Failed to create Arrow table for protein groups" << std::endl;
    return false;
  }
  return writeArrowTableToParquet_(table, filename, "protein_groups", config);
}


std::shared_ptr<arrow::Table> ProteinIdentificationArrowExport::exportSearchParamsToArrow(
  const std::vector<ProteinIdentification>& protein_identifications)
{
  // -- Simple column builders --
  arrow::StringBuilder run_identifier_builder, search_engine_builder;
  arrow::StringBuilder search_engine_version_builder;
  arrow::StringBuilder inference_engine_builder, inference_engine_version_builder;
  arrow::StringBuilder date_builder, score_type_builder;
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
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: run_identifier_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = search_engine_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: search_engine_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = search_engine_version_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: search_engine_version_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = inference_engine_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: inference_engine_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = inference_engine_version_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: inference_engine_version_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = date_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: date_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = score_type_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: score_type_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = higher_score_better_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: higher_score_better_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = significance_threshold_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: significance_threshold_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = db_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: db_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = db_version_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: db_version_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = taxonomy_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: taxonomy_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = charges_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: charges_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = mass_type_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: mass_type_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = precursor_mass_tolerance_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: precursor_mass_tolerance_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = precursor_mass_tolerance_ppm_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: precursor_mass_tolerance_ppm_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = fragment_mass_tolerance_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: fragment_mass_tolerance_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = fragment_mass_tolerance_ppm_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: fragment_mass_tolerance_ppm_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = digestion_enzyme_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: digestion_enzyme_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = enzyme_term_specificity_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: enzyme_term_specificity_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = missed_cleavages_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: missed_cleavages_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }

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

    // === date (nullable, null if empty) ===
    String date_str = prot_id.getDateTime().toString("yyyy-MM-ddThh:mm:ss");
    if (date_str.empty())
    {
      (void)date_builder.AppendNull();
    }
    else
    {
      (void)date_builder.Append(date_str);
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
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: run_identifier_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = search_engine_builder.Finish(&arr_search_engine);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: search_engine_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = search_engine_version_builder.Finish(&arr_se_version);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: search_engine_version_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = inference_engine_builder.Finish(&arr_inf_engine);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: inference_engine_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = inference_engine_version_builder.Finish(&arr_inf_version);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: inference_engine_version_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = date_builder.Finish(&arr_date);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: date_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = score_type_builder.Finish(&arr_score_type);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: score_type_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = higher_score_better_builder.Finish(&arr_higher_better);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: higher_score_better_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = significance_threshold_builder.Finish(&arr_sig_threshold);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: significance_threshold_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = db_builder.Finish(&arr_db);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: db_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = db_version_builder.Finish(&arr_db_version);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: db_version_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = taxonomy_builder.Finish(&arr_taxonomy);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: taxonomy_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = charges_builder.Finish(&arr_charges);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: charges_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = mass_type_builder.Finish(&arr_mass_type);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: mass_type_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = precursor_mass_tolerance_builder.Finish(&arr_precursor_tol);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: precursor_mass_tolerance_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = precursor_mass_tolerance_ppm_builder.Finish(&arr_precursor_ppm);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: precursor_mass_tolerance_ppm_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = fragment_mass_tolerance_builder.Finish(&arr_fragment_tol);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: fragment_mass_tolerance_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = fragment_mass_tolerance_ppm_builder.Finish(&arr_fragment_ppm);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: fragment_mass_tolerance_ppm_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = digestion_enzyme_builder.Finish(&arr_enzyme);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: digestion_enzyme_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = enzyme_term_specificity_builder.Finish(&arr_enzyme_spec);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: enzyme_term_specificity_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = missed_cleavages_builder.Finish(&arr_missed_cleavages);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: missed_cleavages_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = fixed_modifications_builder.Finish(&arr_fixed_mods);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: fixed_modifications_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = variable_modifications_builder.Finish(&arr_var_mods);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: variable_modifications_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = primary_ms_run_paths_builder.Finish(&arr_ms_run_paths);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: primary_ms_run_paths_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = metavalues_builder.Finish(&arr_metavalues);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: metavalues_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = sp_metavalues_builder.Finish(&arr_sp_metavalues);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: sp_metavalues_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }

  // Build schema (26 columns)
  auto schema = arrow::schema({
    arrow::field("run_identifier", arrow::utf8(), /*nullable=*/false),
    arrow::field("search_engine", arrow::utf8(), /*nullable=*/false),
    arrow::field("search_engine_version", arrow::utf8(), /*nullable=*/true),
    arrow::field("inference_engine", arrow::utf8(), /*nullable=*/true),
    arrow::field("inference_engine_version", arrow::utf8(), /*nullable=*/true),
    arrow::field("date", arrow::utf8(), /*nullable=*/true),
    arrow::field("score_type", arrow::utf8(), /*nullable=*/false),
    arrow::field("higher_score_better", arrow::boolean(), /*nullable=*/false),
    arrow::field("significance_threshold", arrow::float64(), /*nullable=*/true),
    arrow::field("db", arrow::utf8(), /*nullable=*/true),
    arrow::field("db_version", arrow::utf8(), /*nullable=*/true),
    arrow::field("taxonomy", arrow::utf8(), /*nullable=*/true),
    arrow::field("charges", arrow::utf8(), /*nullable=*/true),
    arrow::field("mass_type", arrow::utf8(), /*nullable=*/false),
    arrow::field("precursor_mass_tolerance", arrow::float64(), /*nullable=*/false),
    arrow::field("precursor_mass_tolerance_ppm", arrow::boolean(), /*nullable=*/false),
    arrow::field("fragment_mass_tolerance", arrow::float64(), /*nullable=*/false),
    arrow::field("fragment_mass_tolerance_ppm", arrow::boolean(), /*nullable=*/false),
    arrow::field("digestion_enzyme", arrow::utf8(), /*nullable=*/true),
    arrow::field("enzyme_term_specificity", arrow::utf8(), /*nullable=*/true),
    arrow::field("missed_cleavages", arrow::int32(), /*nullable=*/false),
    arrow::field("fixed_modifications", arrow::list(arrow::utf8()), /*nullable=*/false),
    arrow::field("variable_modifications", arrow::list(arrow::utf8()), /*nullable=*/false),
    arrow::field("primary_ms_run_paths", arrow::list(arrow::utf8()), /*nullable=*/false),
    arrow::field("metavalues", metavalues_builder.type(), /*nullable=*/false),
    arrow::field("sp_metavalues", sp_metavalues_builder.type(), /*nullable=*/false),
  });

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

  return table;
}


bool ProteinIdentificationArrowExport::exportSearchParamsToParquet(
  const std::vector<ProteinIdentification>& protein_identifications,
  const String& filename,
  const ParquetWriteConfig& config)
{
  auto table = exportSearchParamsToArrow(protein_identifications);
  if (!table)
  {
    OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport: Failed to create Arrow table for search params" << std::endl;
    return false;
  }
  return writeArrowTableToParquet_(table, filename, "search_params", config);
}


} // namespace OpenMS

#endif // WITH_PARQUET
