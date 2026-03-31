// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ProteinGroupArrowExport.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>

#include <arrow/api.h>
#include <arrow/builder.h>
#include <arrow/io/file.h>
#include <parquet/arrow/writer.h>
#include <parquet/properties.h>

#include <cstdio>
#include <cstring>
#include <map>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace OpenMS
{

std::shared_ptr<arrow::Table> ProteinGroupArrowExport::exportToArrow(const ConsensusMap& cmap)
{
  if (cmap.getProteinIdentifications().empty())
  {
    OPENMS_LOG_WARN << "ProteinGroupArrowExport: No protein identifications found" << std::endl;
    return nullptr;
  }

  const auto& prot_id = cmap.getProteinIdentifications()[0];
  const auto& groups = prot_id.getIndistinguishableProteins();

  if (groups.empty())
  {
    OPENMS_LOG_WARN << "ProteinGroupArrowExport: No indistinguishable protein groups found" << std::endl;
    return nullptr;
  }

  // Build accession -> ProteinHit lookup
  std::unordered_map<std::string, const ProteinHit*> hit_lookup;
  for (const auto& hit : prot_id.getHits())
  {
    hit_lookup[hit.getAccession()] = &hit;
  }

  // Estimate number of rows: one per group per unique run file (or 1 if no quant data)
  Size estimated_rows = 0;
  for (const auto& group : groups)
  {
    if (group.getFloatDataArrays().size() >= 4 && !group.getStringDataArrays().empty())
    {
      // Count unique filenames
      const auto& filenames = group.getStringDataArrays()[0];
      std::unordered_set<std::string> unique_files;
      for (const auto& fn : filenames) { unique_files.insert(fn); }
      estimated_rows += std::max(unique_files.size(), (size_t)1);
    }
    else
    {
      estimated_rows += 1;
    }
  }

  // -- Simple column builders --
  arrow::StringBuilder anchor_protein_builder, run_file_name_builder;
  arrow::DoubleBuilder global_qvalue_builder, pg_qvalue_builder, gg_qvalue_builder;
  arrow::FloatBuilder sequence_coverage_builder, molecular_weight_builder;
  arrow::BooleanBuilder is_decoy_builder, contaminant_builder;

  // -- list<utf8> builders --
  auto pg_acc_vb = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder pg_accessions_builder(arrow::default_memory_pool(), pg_acc_vb);

  auto pg_names_vb = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder pg_names_builder(arrow::default_memory_pool(), pg_names_vb);

  auto gg_acc_vb = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder gg_accessions_builder(arrow::default_memory_pool(), gg_acc_vb);

  auto gg_names_vb = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder gg_names_builder(arrow::default_memory_pool(), gg_names_vb);

  // -- intensities: list<struct{label: utf8 not null, intensity: float32 not null}> --
  auto intensity_struct_type = arrow::struct_({
    arrow::field("label", arrow::utf8(), /*nullable=*/false),
    arrow::field("intensity", arrow::float32(), /*nullable=*/false)
  });
  auto int_label_b = std::make_shared<arrow::StringBuilder>();
  auto int_value_b = std::make_shared<arrow::FloatBuilder>();
  auto int_struct_b = std::make_shared<arrow::StructBuilder>(
    intensity_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{int_label_b, int_value_b});
  arrow::ListBuilder intensities_builder(arrow::default_memory_pool(), int_struct_b);

  // -- additional_intensities: list<struct{label: utf8 not null, intensities: list<struct{intensity_name: utf8 not null, intensity_value: float32 not null}> not null}> --
  auto int_pair_type = arrow::struct_({
    arrow::field("intensity_name", arrow::utf8(), /*nullable=*/false),
    arrow::field("intensity_value", arrow::float32(), /*nullable=*/false)
  });
  auto ip_name_b = std::make_shared<arrow::StringBuilder>();
  auto ip_value_b = std::make_shared<arrow::FloatBuilder>();
  auto ip_struct_b = std::make_shared<arrow::StructBuilder>(
    int_pair_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{ip_name_b, ip_value_b});
  auto ip_list_b = std::make_shared<arrow::ListBuilder>(arrow::default_memory_pool(), ip_struct_b);

  auto add_int_struct_type = arrow::struct_({
    arrow::field("label", arrow::utf8(), /*nullable=*/false),
    arrow::field("intensities", arrow::list(int_pair_type), /*nullable=*/false)
  });
  auto ai_label_b = std::make_shared<arrow::StringBuilder>();
  auto ai_struct_b = std::make_shared<arrow::StructBuilder>(
    add_int_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{ai_label_b, ip_list_b});
  arrow::ListBuilder additional_intensities_builder(arrow::default_memory_pool(), ai_struct_b);

  // -- peptides: list<struct{protein_name: utf8 not null, peptide_count: int32 not null}> --
  auto pep_struct_type = arrow::struct_({
    arrow::field("protein_name", arrow::utf8(), /*nullable=*/false),
    arrow::field("peptide_count", arrow::int32(), /*nullable=*/false)
  });
  auto pep_name_b = std::make_shared<arrow::StringBuilder>();
  auto pep_count_b = std::make_shared<arrow::Int32Builder>();
  auto pep_struct_b = std::make_shared<arrow::StructBuilder>(
    pep_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{pep_name_b, pep_count_b});
  arrow::ListBuilder peptides_builder(arrow::default_memory_pool(), pep_struct_b);

  // -- peptide_counts: struct{unique_sequences: int32 not null, total_sequences: int32 not null} --
  auto pep_counts_type = arrow::struct_({
    arrow::field("unique_sequences", arrow::int32(), /*nullable=*/false),
    arrow::field("total_sequences", arrow::int32(), /*nullable=*/false)
  });
  auto pc_unique_b = std::make_shared<arrow::Int32Builder>();
  auto pc_total_b = std::make_shared<arrow::Int32Builder>();
  arrow::StructBuilder peptide_counts_builder(
    pep_counts_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{pc_unique_b, pc_total_b});

  // -- feature_counts: struct{unique_features: int32 not null, total_features: int32 not null} --
  auto feat_counts_type = arrow::struct_({
    arrow::field("unique_features", arrow::int32(), /*nullable=*/false),
    arrow::field("total_features", arrow::int32(), /*nullable=*/false)
  });
  auto fc_unique_b = std::make_shared<arrow::Int32Builder>();
  auto fc_total_b = std::make_shared<arrow::Int32Builder>();
  arrow::StructBuilder feature_counts_builder(
    feat_counts_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{fc_unique_b, fc_total_b});

  // -- additional_scores: list<struct{score_name: utf8 not null, score_value: float64 not null, higher_better: bool}> --
  auto as_struct_type = arrow::struct_({
    arrow::field("score_name", arrow::utf8(), /*nullable=*/false),
    arrow::field("score_value", arrow::float64(), /*nullable=*/false),
    arrow::field("higher_better", arrow::boolean())
  });
  auto as_name_b = std::make_shared<arrow::StringBuilder>();
  auto as_value_b = std::make_shared<arrow::DoubleBuilder>();
  auto as_hb_b = std::make_shared<arrow::BooleanBuilder>();
  auto as_struct_b = std::make_shared<arrow::StructBuilder>(
    as_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{as_name_b, as_value_b, as_hb_b});
  arrow::ListBuilder additional_scores_builder(arrow::default_memory_pool(), as_struct_b);

  // -- cv_params: list<struct{cv_name: utf8 not null, cv_value: utf8 not null}> --
  auto cv_struct_type = arrow::struct_({
    arrow::field("cv_name", arrow::utf8(), /*nullable=*/false),
    arrow::field("cv_value", arrow::utf8(), /*nullable=*/false)
  });
  auto cv_name_b = std::make_shared<arrow::StringBuilder>();
  auto cv_value_b = std::make_shared<arrow::StringBuilder>();
  auto cv_struct_b = std::make_shared<arrow::StructBuilder>(
    cv_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{cv_name_b, cv_value_b});
  arrow::ListBuilder cv_params_builder(arrow::default_memory_pool(), cv_struct_b);

  // Reserve capacity
  arrow::Status status;
  status = anchor_protein_builder.Reserve(estimated_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: anchor_protein Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = run_file_name_builder.Reserve(estimated_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: run_file_name Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = global_qvalue_builder.Reserve(estimated_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: global_qvalue Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = pg_qvalue_builder.Reserve(estimated_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: pg_qvalue Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = gg_qvalue_builder.Reserve(estimated_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: gg_qvalue Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = sequence_coverage_builder.Reserve(estimated_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: sequence_coverage Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = molecular_weight_builder.Reserve(estimated_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: molecular_weight Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = is_decoy_builder.Reserve(estimated_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: is_decoy Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = contaminant_builder.Reserve(estimated_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: contaminant Reserve failed: " << status.ToString() << std::endl; return nullptr; }

  // Helper lambda to emit one row for a protein group + run file
  auto emitRow = [&](const ProteinIdentification::ProteinGroup& group,
                     const std::string& run_file,
                     const std::vector<std::pair<std::string, float>>& channel_intensities,
                     int distinct_peptides_for_run)
  {
    const std::string& anchor = group.accessions.empty() ? "" : group.accessions[0];

    // anchor_protein
    (void)anchor_protein_builder.Append(anchor);

    // run_file_name
    (void)run_file_name_builder.Append(run_file);

    // pg_accessions
    (void)pg_accessions_builder.Append();
    for (const auto& acc : group.accessions)
    {
      (void)pg_acc_vb->Append(acc);
    }

    // pg_names
    (void)pg_names_builder.Append();
    for (const auto& acc : group.accessions)
    {
      auto it = hit_lookup.find(acc);
      if (it != hit_lookup.end())
      {
        (void)pg_names_vb->Append(it->second->getDescription());
      }
      else
      {
        (void)pg_names_vb->Append("");
      }
    }

    // gg_accessions, gg_names (keep lists symmetric so they can be zipped)
    (void)gg_accessions_builder.Append();
    (void)gg_names_builder.Append();
    for (const auto& acc : group.accessions)
    {
      auto it = hit_lookup.find(acc);
      if (it != hit_lookup.end())
      {
        bool has_ga = it->second->metaValueExists("gene_accession");
        bool has_gn = it->second->metaValueExists("gene_name");
        if (has_ga || has_gn)
        {
          (void)gg_acc_vb->Append(has_ga ? it->second->getMetaValue("gene_accession").toString() : "");
          (void)gg_names_vb->Append(has_gn ? it->second->getMetaValue("gene_name").toString() : "");
        }
      }
    }

    // gg_qvalue - not available
    (void)gg_qvalue_builder.AppendNull();

    // global_qvalue
    (void)global_qvalue_builder.Append(group.probability);

    // pg_qvalue - not per-run in OpenMS
    (void)pg_qvalue_builder.AppendNull();

    // intensities
    (void)intensities_builder.Append();
    for (const auto& [label, intensity] : channel_intensities)
    {
      (void)int_struct_b->Append();
      (void)int_label_b->Append(label);
      (void)int_value_b->Append(intensity);
    }

    // additional_intensities - empty for now
    (void)additional_intensities_builder.Append();

    // is_decoy
    auto anchor_it = hit_lookup.find(anchor);
    if (anchor_it != hit_lookup.end())
    {
      (void)is_decoy_builder.Append(anchor_it->second->isDecoy());
    }
    else
    {
      (void)is_decoy_builder.Append(false);
    }

    // contaminant - not tracked
    (void)contaminant_builder.AppendNull();

    // peptides - empty for now
    (void)peptides_builder.Append();

    // peptide_counts
    if (distinct_peptides_for_run >= 0)
    {
      (void)peptide_counts_builder.Append();
      (void)pc_unique_b->Append(distinct_peptides_for_run);
      (void)pc_total_b->Append(0);
    }
    else
    {
      (void)peptide_counts_builder.AppendNull();
    }

    // feature_counts - not available
    (void)feature_counts_builder.AppendNull();

    // sequence_coverage
    if (anchor_it != hit_lookup.end() && anchor_it->second->getCoverage() != ProteinHit::COVERAGE_UNKNOWN)
    {
      (void)sequence_coverage_builder.Append(static_cast<float>(anchor_it->second->getCoverage()));
    }
    else
    {
      (void)sequence_coverage_builder.AppendNull();
    }

    // molecular_weight - not available
    (void)molecular_weight_builder.AppendNull();

    // additional_scores - anchor protein score
    (void)additional_scores_builder.Append();
    if (anchor_it != hit_lookup.end())
    {
      (void)as_struct_b->Append();
      (void)as_name_b->Append(prot_id.getScoreType());
      (void)as_value_b->Append(anchor_it->second->getScore());
      (void)as_hb_b->Append(prot_id.isHigherScoreBetter());
    }

    // cv_params - empty for now
    (void)cv_params_builder.Append();
  };

  // Iterate protein groups and emit rows
  for (const auto& group : groups)
  {
    bool has_quant = group.getFloatDataArrays().size() >= 4
                     && !group.getStringDataArrays().empty()
                     && !group.getStringDataArrays()[0].empty();

    if (has_quant)
    {
      const auto& abundances = group.getFloatDataArrays()[3]; // file_channel_level_abundance
      const auto& filenames = group.getStringDataArrays()[0]; // file_channel_level_filename
      const auto& channels = group.getIntegerDataArrays()[0]; // file_channel_level_channel

      // Group by filename: filename -> [(channel_label, intensity)]
      std::map<std::string, std::vector<std::pair<std::string, float>>> file_intensities;
      for (Size i = 0; i < filenames.size() && i < abundances.size() && i < channels.size(); ++i)
      {
        file_intensities[filenames[i]].emplace_back(std::to_string(channels[i]), abundances[i]);
      }

      // Emit one row per unique filename
      for (const auto& [filename, channel_int] : file_intensities)
      {
        // We don't have per-run distinct_peptides directly - use -1 to indicate unknown
        emitRow(group, filename, channel_int, -1);
      }
    }
    else
    {
      // Unquantified group: single row with empty run_file_name
      emitRow(group, "", {}, -1);
    }
  }

  // Finalize all arrays
  std::shared_ptr<arrow::Array> arr_pg_acc, arr_pg_names, arr_gg_acc, arr_gg_names;
  std::shared_ptr<arrow::Array> arr_gg_qval, arr_anchor, arr_run_file;
  std::shared_ptr<arrow::Array> arr_global_qval, arr_pg_qval;
  std::shared_ptr<arrow::Array> arr_intensities, arr_add_intensities;
  std::shared_ptr<arrow::Array> arr_is_decoy, arr_contaminant;
  std::shared_ptr<arrow::Array> arr_peptides, arr_pep_counts, arr_feat_counts;
  std::shared_ptr<arrow::Array> arr_seq_cov, arr_mol_weight;
  std::shared_ptr<arrow::Array> arr_add_scores, arr_cv_params;

  status = pg_accessions_builder.Finish(&arr_pg_acc);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: pg_accessions Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = pg_names_builder.Finish(&arr_pg_names);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: pg_names Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = gg_accessions_builder.Finish(&arr_gg_acc);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: gg_accessions Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = gg_names_builder.Finish(&arr_gg_names);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: gg_names Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = gg_qvalue_builder.Finish(&arr_gg_qval);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: gg_qvalue Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = anchor_protein_builder.Finish(&arr_anchor);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: anchor_protein Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = run_file_name_builder.Finish(&arr_run_file);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: run_file_name Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = global_qvalue_builder.Finish(&arr_global_qval);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: global_qvalue Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = pg_qvalue_builder.Finish(&arr_pg_qval);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: pg_qvalue Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = intensities_builder.Finish(&arr_intensities);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: intensities Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = additional_intensities_builder.Finish(&arr_add_intensities);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: additional_intensities Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = is_decoy_builder.Finish(&arr_is_decoy);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: is_decoy Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = contaminant_builder.Finish(&arr_contaminant);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: contaminant Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = peptides_builder.Finish(&arr_peptides);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: peptides Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = peptide_counts_builder.Finish(&arr_pep_counts);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: peptide_counts Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = feature_counts_builder.Finish(&arr_feat_counts);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: feature_counts Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = sequence_coverage_builder.Finish(&arr_seq_cov);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: sequence_coverage Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = molecular_weight_builder.Finish(&arr_mol_weight);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: molecular_weight Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = additional_scores_builder.Finish(&arr_add_scores);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: additional_scores Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = cv_params_builder.Finish(&arr_cv_params);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport: cv_params Finish failed: " << status.ToString() << std::endl; return nullptr; }

  // Build schema matching quantms pg spec field order
  auto schema = arrow::schema({
    arrow::field("pg_accessions", arrow::list(arrow::utf8()), /*nullable=*/false),
    arrow::field("pg_names", arrow::list(arrow::utf8())),
    arrow::field("gg_accessions", arrow::list(arrow::utf8())),
    arrow::field("gg_names", arrow::list(arrow::utf8())),
    arrow::field("gg_qvalue", arrow::float64()),
    arrow::field("anchor_protein", arrow::utf8(), /*nullable=*/false),
    arrow::field("run_file_name", arrow::utf8(), /*nullable=*/false),
    arrow::field("global_qvalue", arrow::float64()),
    arrow::field("pg_qvalue", arrow::float64()),
    arrow::field("intensities", intensities_builder.type()),
    arrow::field("additional_intensities", additional_intensities_builder.type()),
    arrow::field("is_decoy", arrow::boolean(), /*nullable=*/false),
    arrow::field("contaminant", arrow::boolean()),
    arrow::field("peptides", peptides_builder.type(), /*nullable=*/false),
    arrow::field("peptide_counts", pep_counts_type),
    arrow::field("feature_counts", feat_counts_type),
    arrow::field("sequence_coverage", arrow::float32()),
    arrow::field("molecular_weight", arrow::float32()),
    arrow::field("additional_scores", additional_scores_builder.type()),
    arrow::field("cv_params", cv_params_builder.type()),
  });

  auto table = arrow::Table::Make(schema, {
    arr_pg_acc, arr_pg_names, arr_gg_acc, arr_gg_names,
    arr_gg_qval, arr_anchor, arr_run_file,
    arr_global_qval, arr_pg_qval,
    arr_intensities, arr_add_intensities,
    arr_is_decoy, arr_contaminant,
    arr_peptides, arr_pep_counts, arr_feat_counts,
    arr_seq_cov, arr_mol_weight,
    arr_add_scores, arr_cv_params,
  });

  return table;
}


bool ProteinGroupArrowExport::exportToParquet(
  const ConsensusMap& cmap,
  const String& filename,
  const ParquetWriteConfig& config)
{
  auto table = exportToArrow(cmap);
  if (!table)
  {
    OPENMS_LOG_ERROR << "ProteinGroupArrowExport: Failed to create Arrow table" << std::endl;
    return false;
  }

  // Add QPX metadata to schema
  {
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

    auto metadata = arrow::key_value_metadata({
      {"qpx_version", "1.0"},
      {"creator", "OpenMS"},
      {"file_type", "pg"},
      {"creation_date", DateTime::nowUTC().toString("yyyy-MM-ddThh:mm:ssZ")},
      {"uuid", std::string(buf)},
      {"software_provider", "OpenMS"}
    });
    table = table->ReplaceSchemaMetadata(metadata);
  }

  // Open output file
  auto result = arrow::io::FileOutputStream::Open(filename);
  if (!result.ok())
  {
    OPENMS_LOG_ERROR << "ProteinGroupArrowExport: Failed to open file: " << filename << std::endl;
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

  auto write_status = parquet::arrow::WriteTable(
    *table, arrow::default_memory_pool(), outfile,
    config.row_group_size, writer_props, arrow_props);

  if (!write_status.ok())
  {
    OPENMS_LOG_ERROR << "ProteinGroupArrowExport: Failed to write Parquet: "
                     << write_status.ToString() << std::endl;
    return false;
  }

  return true;
}

} // namespace OpenMS
