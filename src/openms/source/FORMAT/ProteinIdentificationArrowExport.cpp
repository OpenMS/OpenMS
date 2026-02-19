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
  const std::vector<ProteinIdentification>& /*protein_identifications*/)
{
  // TODO: implement in a future task
  OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport::exportProteinGroupsToArrow not yet implemented" << std::endl;
  return nullptr;
}


bool ProteinIdentificationArrowExport::exportProteinGroupsToParquet(
  const std::vector<ProteinIdentification>& /*protein_identifications*/,
  const String& /*filename*/,
  const ParquetWriteConfig& /*config*/)
{
  // TODO: implement in a future task
  OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport::exportProteinGroupsToParquet not yet implemented" << std::endl;
  return false;
}


std::shared_ptr<arrow::Table> ProteinIdentificationArrowExport::exportSearchParamsToArrow(
  const std::vector<ProteinIdentification>& /*protein_identifications*/)
{
  // TODO: implement in a future task
  OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport::exportSearchParamsToArrow not yet implemented" << std::endl;
  return nullptr;
}


bool ProteinIdentificationArrowExport::exportSearchParamsToParquet(
  const std::vector<ProteinIdentification>& /*protein_identifications*/,
  const String& /*filename*/,
  const ParquetWriteConfig& /*config*/)
{
  // TODO: implement in a future task
  OPENMS_LOG_ERROR << "ProteinIdentificationArrowExport::exportSearchParamsToParquet not yet implemented" << std::endl;
  return false;
}


} // namespace OpenMS

#endif // WITH_PARQUET
