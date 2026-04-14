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
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/FORMAT/ArrowIOHelpers.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>

#include <arrow/api.h>
#include <arrow/builder.h>

#include <map>
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


namespace
{
  /// Attach the canonical QPX "pg" file metadata to a table's schema.
  std::shared_ptr<arrow::Table> attachQPXPgMetadata(const std::shared_ptr<arrow::Table>& table)
  {
    auto metadata = arrow::key_value_metadata({
      {"qpx_version", "1.0"},
      {"creator", "OpenMS"},
      {"file_type", "pg"},
      {"creation_date", DateTime::nowUTC().toString("yyyy-MM-ddThh:mm:ssZ")},
      {"uuid", std::string(ArrowIOHelpers::generateUuidV4())},
      {"software_provider", "OpenMS"}
    });
    return table->ReplaceSchemaMetadata(metadata);
  }
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
  return exportToParquet(table, filename, config);
}

bool ProteinGroupArrowExport::exportToParquet(
  const std::shared_ptr<arrow::Table>& table,
  const String& filename,
  const ParquetWriteConfig& config)
{
  if (!table)
  {
    OPENMS_LOG_ERROR << "ProteinGroupArrowExport: null table passed to exportToParquet (" << filename << ")" << std::endl;
    return false;
  }

  // Guard: this overload stamps file_type="pg", so the caller must pass a
  // QPXPgSchema table.
  auto validation = ArrowSchemaValidation::validate(table, QPXPgSchema::schema(), ArrowSchemaValidation::Mode::Strict);
  if (!validation.valid)
  {
    OPENMS_LOG_ERROR << "ProteinGroupArrowExport: table schema does not match QPXPgSchema ("
                     << filename << "): " << validation.toString() << std::endl;
    return false;
  }

  return ArrowIOHelpers::writeTableToParquet(attachQPXPgMetadata(table), filename, config);
}


std::shared_ptr<arrow::Table> ProteinGroupArrowExport::exportToArrow(
  const std::vector<ProteinIdentification>& protein_identifications,
  const PeptideIdentificationList& peptide_identifications)
{
  auto target_schema = QPXPgSchema::schema();

  if (protein_identifications.empty())
  {
    OPENMS_LOG_WARN << "ProteinGroupArrowExport (id-only): No protein identifications found" << std::endl;
    return arrow::Table::MakeEmpty(target_schema).ValueOrDie();
  }

  // Compute total estimated rows across all protein identifications
  Size estimated_rows = 0;
  for (const auto& pid : protein_identifications)
  {
    estimated_rows += pid.getIndistinguishableProteins().size();
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

  // -- intensities: null (not available for id-only) --
  // We still need a properly typed builder to match schema; append nulls per row
  auto intensities_value_type = std::dynamic_pointer_cast<arrow::ListType>(QPXPgSchema::intensitiesType())->value_type();
  auto int_label_b = std::make_shared<arrow::StringBuilder>();
  auto int_value_b = std::make_shared<arrow::FloatBuilder>();
  auto int_struct_b = std::make_shared<arrow::StructBuilder>(
    intensities_value_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{int_label_b, int_value_b});
  arrow::ListBuilder intensities_builder(arrow::default_memory_pool(), int_struct_b, QPXPgSchema::intensitiesType());

  // -- additional_intensities: null (not available for id-only) --
  auto add_int_value_type = std::dynamic_pointer_cast<arrow::ListType>(QPXPgSchema::additionalIntensitiesType())->value_type();
  // inner list type for the "intensities" field of the struct
  auto int_pair_type = arrow::struct_({
    arrow::field("intensity_name", arrow::utf8(), /*nullable=*/false),
    arrow::field("intensity_value", arrow::float32(), /*nullable=*/false)
  });
  auto ip_name_b = std::make_shared<arrow::StringBuilder>();
  auto ip_value_b = std::make_shared<arrow::FloatBuilder>();
  auto ip_struct_b = std::make_shared<arrow::StructBuilder>(
    int_pair_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{ip_name_b, ip_value_b});
  auto ip_list_b = std::make_shared<arrow::ListBuilder>(arrow::default_memory_pool(), ip_struct_b, arrow::list(int_pair_type));
  auto ai_label_b = std::make_shared<arrow::StringBuilder>();
  auto ai_struct_b = std::make_shared<arrow::StructBuilder>(
    add_int_value_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{ai_label_b, ip_list_b});
  arrow::ListBuilder additional_intensities_builder(arrow::default_memory_pool(), ai_struct_b, QPXPgSchema::additionalIntensitiesType());

  // -- peptides: list<struct{protein_name: utf8, peptide_count: int32}> --
  auto pep_value_type = std::dynamic_pointer_cast<arrow::ListType>(QPXPgSchema::peptidesType())->value_type();
  auto pep_name_b = std::make_shared<arrow::StringBuilder>();
  auto pep_count_b = std::make_shared<arrow::Int32Builder>();
  auto pep_struct_b = std::make_shared<arrow::StructBuilder>(
    pep_value_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{pep_name_b, pep_count_b});
  arrow::ListBuilder peptides_builder(arrow::default_memory_pool(), pep_struct_b, QPXPgSchema::peptidesType());

  // -- peptide_counts: struct{unique_sequences, total_sequences} --
  auto pep_counts_type = QPXPgSchema::peptideCountsType();
  auto pc_unique_b = std::make_shared<arrow::Int32Builder>();
  auto pc_total_b = std::make_shared<arrow::Int32Builder>();
  arrow::StructBuilder peptide_counts_builder(
    pep_counts_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{pc_unique_b, pc_total_b});

  // -- feature_counts: null (not available for id-only) --
  auto feat_counts_type = QPXPgSchema::featureCountsType();
  auto fc_unique_b = std::make_shared<arrow::Int32Builder>();
  auto fc_total_b = std::make_shared<arrow::Int32Builder>();
  arrow::StructBuilder feature_counts_builder(
    feat_counts_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{fc_unique_b, fc_total_b});

  // -- additional_scores: list<struct{score_name, score_value, higher_better}> --
  auto as_type = QPXPgSchema::additionalScoresType();
  auto as_value_type = std::dynamic_pointer_cast<arrow::ListType>(as_type)->value_type();
  auto as_name_b = std::make_shared<arrow::StringBuilder>();
  auto as_value_b = std::make_shared<arrow::DoubleBuilder>();
  auto as_hb_b = std::make_shared<arrow::BooleanBuilder>();
  auto as_struct_b = std::make_shared<arrow::StructBuilder>(
    as_value_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{as_name_b, as_value_b, as_hb_b});
  arrow::ListBuilder additional_scores_builder(arrow::default_memory_pool(), as_struct_b, as_type);

  // -- cv_params: list<struct{cv_name, cv_value}> --
  auto cv_type = QPXPgSchema::cvParamsType();
  auto cv_value_type = std::dynamic_pointer_cast<arrow::ListType>(cv_type)->value_type();
  auto cv_name_b = std::make_shared<arrow::StringBuilder>();
  auto cv_value_b = std::make_shared<arrow::StringBuilder>();
  auto cv_struct_b = std::make_shared<arrow::StructBuilder>(
    cv_value_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{cv_name_b, cv_value_b});
  arrow::ListBuilder cv_params_builder(arrow::default_memory_pool(), cv_struct_b, cv_type);

  // Reserve capacity
  arrow::Status status;
  status = anchor_protein_builder.Reserve(estimated_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): anchor_protein Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = run_file_name_builder.Reserve(estimated_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): run_file_name Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = global_qvalue_builder.Reserve(estimated_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): global_qvalue Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = is_decoy_builder.Reserve(estimated_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): is_decoy Reserve failed: " << status.ToString() << std::endl; return nullptr; }

  // Iterate over all protein identifications; for each, build scoped lookup
  // structures and append one row per indistinguishable protein group.
  Size total_groups_found = 0;
  for (const auto& prot_id : protein_identifications)
  {
    const auto& groups = prot_id.getIndistinguishableProteins();
    if (groups.empty()) continue;
    total_groups_found += groups.size();

    // Derive run file name from primary MS run path for this prot_id
    StringList run_paths;
    prot_id.getPrimaryMSRunPath(run_paths);
    std::string run_file = run_paths.empty() ? "" : std::string(run_paths[0]);

    // Build accession -> ProteinHit lookup scoped to this prot_id
    std::unordered_map<std::string, const ProteinHit*> hit_lookup;
    for (const auto& hit : prot_id.getHits())
    {
      hit_lookup[hit.getAccession()] = &hit;
    }

    // Build per-protein peptide counts from peptide_identifications whose
    // identifier matches this prot_id's identifier
    std::unordered_map<std::string, std::unordered_set<std::string>> acc_to_peptides;
    for (const auto& pep_id : peptide_identifications)
    {
      if (pep_id.getIdentifier() != prot_id.getIdentifier()) continue;
      if (pep_id.getHits().empty()) continue;
      // Use the best hit (first after assuming sorted by score)
      const auto& best_hit = pep_id.getHits()[0];
      std::string seq = best_hit.getSequence().toString();
      for (const auto& ev : best_hit.getPeptideEvidences())
      {
        const std::string& acc = ev.getProteinAccession();
        if (!acc.empty())
        {
          acc_to_peptides[acc].insert(seq);
        }
      }
    }

    // Iterate protein groups for this prot_id and emit rows
    for (const auto& group : groups)
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

      // gg_accessions, gg_names (gene groups — one entry per group member, keep
      // parallel to pg_accessions so consumers can zip positionally; empty string
      // when the meta value is missing).
      (void)gg_accessions_builder.Append();
      (void)gg_names_builder.Append();
      for (const auto& acc : group.accessions)
      {
        auto it = hit_lookup.find(acc);
        bool has_ga = (it != hit_lookup.end()) && it->second->metaValueExists("gene_accession");
        bool has_gn = (it != hit_lookup.end()) && it->second->metaValueExists("gene_name");
        (void)gg_acc_vb->Append(has_ga ? it->second->getMetaValue("gene_accession").toString() : "");
        (void)gg_names_vb->Append(has_gn ? it->second->getMetaValue("gene_name").toString() : "");
      }

      // gg_qvalue - not available
      (void)gg_qvalue_builder.AppendNull();

      // global_qvalue — use group probability
      (void)global_qvalue_builder.Append(group.probability);

      // pg_qvalue - not available for id-only
      (void)pg_qvalue_builder.AppendNull();

      // intensities — null for id-only
      (void)intensities_builder.AppendNull();

      // additional_intensities — null for id-only
      (void)additional_intensities_builder.AppendNull();

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

      // peptides — one {protein_name, peptide_count} per group member
      (void)peptides_builder.Append();
      int total_sequences = 0;
      std::unordered_set<std::string> group_unique_peptides;
      for (const auto& acc : group.accessions)
      {
        auto pit = acc_to_peptides.find(acc);
        int count = (pit != acc_to_peptides.end()) ? static_cast<int>(pit->second.size()) : 0;
        (void)pep_struct_b->Append();
        (void)pep_name_b->Append(acc);
        (void)pep_count_b->Append(count);
        total_sequences += count;
        if (pit != acc_to_peptides.end())
        {
          group_unique_peptides.insert(pit->second.begin(), pit->second.end());
        }
      }

      // peptide_counts — unique sequences across the union of group members,
      // total_sequences is the sum of per-protein counts (double-counts shared peptides).
      {
        (void)peptide_counts_builder.Append();
        (void)pc_unique_b->Append(static_cast<int>(group_unique_peptides.size()));
        (void)pc_total_b->Append(total_sequences);
      }

      // feature_counts - not available for id-only
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

      // additional_scores — anchor protein score
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
    }
  }

  if (total_groups_found == 0)
  {
    OPENMS_LOG_WARN << "ProteinGroupArrowExport (id-only): No indistinguishable protein groups found" << std::endl;
    return arrow::Table::MakeEmpty(target_schema).ValueOrDie();
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
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): pg_accessions Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = pg_names_builder.Finish(&arr_pg_names);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): pg_names Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = gg_accessions_builder.Finish(&arr_gg_acc);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): gg_accessions Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = gg_names_builder.Finish(&arr_gg_names);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): gg_names Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = gg_qvalue_builder.Finish(&arr_gg_qval);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): gg_qvalue Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = anchor_protein_builder.Finish(&arr_anchor);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): anchor_protein Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = run_file_name_builder.Finish(&arr_run_file);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): run_file_name Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = global_qvalue_builder.Finish(&arr_global_qval);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): global_qvalue Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = pg_qvalue_builder.Finish(&arr_pg_qval);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): pg_qvalue Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = intensities_builder.Finish(&arr_intensities);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): intensities Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = additional_intensities_builder.Finish(&arr_add_intensities);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): additional_intensities Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = is_decoy_builder.Finish(&arr_is_decoy);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): is_decoy Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = contaminant_builder.Finish(&arr_contaminant);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): contaminant Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = peptides_builder.Finish(&arr_peptides);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): peptides Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = peptide_counts_builder.Finish(&arr_pep_counts);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): peptide_counts Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = feature_counts_builder.Finish(&arr_feat_counts);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): feature_counts Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = sequence_coverage_builder.Finish(&arr_seq_cov);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): sequence_coverage Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = molecular_weight_builder.Finish(&arr_mol_weight);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): molecular_weight Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = additional_scores_builder.Finish(&arr_add_scores);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): additional_scores Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = cv_params_builder.Finish(&arr_cv_params);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): cv_params Finish failed: " << status.ToString() << std::endl; return nullptr; }

  auto table = arrow::Table::Make(target_schema, {
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
  const std::vector<ProteinIdentification>& protein_identifications,
  const PeptideIdentificationList& peptide_identifications,
  const String& filename,
  const ParquetWriteConfig& config)
{
  auto table = exportToArrow(protein_identifications, peptide_identifications);
  if (!table)
  {
    OPENMS_LOG_ERROR << "ProteinGroupArrowExport (id-only): Failed to create Arrow table" << std::endl;
    return false;
  }
  return exportToParquet(table, filename, config);
}

} // namespace OpenMS
