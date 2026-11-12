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
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FORMAT/ArrowIOHelpers.h>
#include "ArrowFileIOHelpers.h"
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/FORMAT/ModificationDefinitionIO.h>

#include <arrow/api.h>
#include <arrow/builder.h>
#include <arrow/io/file.h>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>
#include <parquet/properties.h>

#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace OpenMS
{

namespace // anonymous
{
  /// Read a list<utf8> column at a given row into a vector of Strings.
  std::vector<std::string> readStringList_(const std::shared_ptr<arrow::Array>& array, int64_t row)
  {
    std::vector<std::string> result;
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
    const std::string& run_id = prot_id.getIdentifier();

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
      const std::string& seq = hit.getSequence();
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
        std::string desc = hit.getDescription();
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
      Internal::ArrowFileIO::appendMetaValues(hit, mv_name_b, mv_value_b, mv_type_b, mv_struct_b, excluded_hit_mvs);

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
  const std::string& filename,
  const ParquetWriteConfig& config)
{
  auto table = exportProteinsToArrow(protein_identifications);
  if (!table)
  {
    OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: Failed to create Arrow table for proteins" << std::endl;
    return false;
  }
  return Internal::ArrowFileIO::writeParquetTable(table, filename, "proteins", config, "ProteinIdentificationArrowIO");
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
                          const std::string& run_id,
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
    const std::string& run_id = prot_id.getIdentifier();

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
  const std::string& filename,
  const ParquetWriteConfig& config)
{
  auto table = exportProteinGroupsToArrow(protein_identifications);
  if (!table)
  {
    OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: Failed to create Arrow table for protein groups" << std::endl;
    return false;
  }
  return Internal::ArrowFileIO::writeParquetTable(table, filename, "protein_groups", config, "ProteinIdentificationArrowIO");
}


std::shared_ptr<arrow::Table> ProteinIdentificationArrowIO::exportSearchParamsToArrow(
  const std::vector<ProteinIdentification>& protein_identifications,
  const std::map<std::string, std::string>& definitions_by_run)
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
    const std::string& se_version = prot_id.getSearchEngineVersion();
    if (se_version.empty())
    {
      (void)search_engine_version_builder.AppendNull();
    }
    else
    {
      (void)search_engine_version_builder.Append(se_version);
    }

    // === inference_engine (nullable, null if empty) ===
    const std::string inference_engine = prot_id.getInferenceEngine();
    if (inference_engine.empty())
    {
      (void)inference_engine_builder.AppendNull();
    }
    else
    {
      (void)inference_engine_builder.Append(inference_engine);
    }

    // === inference_engine_version (nullable, null if empty) ===
    const std::string inference_engine_version = prot_id.getInferenceEngineVersion();
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
      std::string date_str = dt.toString("yyyy-MM-ddThh:mm:ss");
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
    const std::string& enzyme_name = sp.digestion_enzyme.getName();
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
    Internal::ArrowFileIO::appendMetaValues(prot_id, mv_name_b, mv_value_b, mv_type_b, mv_struct_b, no_exclusions);

    // === sp_metavalues (not nullable) - from SearchParameters ===
    (void)sp_metavalues_builder.Append();
    const auto defs_it = definitions_by_run.find(prot_id.getIdentifier());
    if (defs_it == definitions_by_run.end())
    {
      Internal::ArrowFileIO::appendMetaValues(sp, sp_mv_name_b, sp_mv_value_b, sp_mv_type_b, sp_mv_struct_b, no_exclusions);
    }
    else
    { // the freshly collected modification definitions replace an inbound value
      static const std::unordered_set<std::string> exclude_definitions{Constants::UserParam::MODIFICATION_DEFINITIONS};
      Internal::ArrowFileIO::appendMetaValues(sp, sp_mv_name_b, sp_mv_value_b, sp_mv_type_b, sp_mv_struct_b, exclude_definitions);
      (void)sp_mv_struct_b->Append();
      (void)sp_mv_name_b->Append(Constants::UserParam::MODIFICATION_DEFINITIONS);
      (void)sp_mv_value_b->Append(defs_it->second);
      (void)sp_mv_type_b->Append("string");
    }

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
  const std::string& filename,
  const ParquetWriteConfig& config,
  const std::map<std::string, std::string>& definitions_by_run)
{
  auto table = exportSearchParamsToArrow(protein_identifications, definitions_by_run);
  if (!table)
  {
    OPENMS_LOG_ERROR << "ProteinIdentificationArrowIO: Failed to create Arrow table for search params" << std::endl;
    return false;
  }
  return Internal::ArrowFileIO::writeParquetTable(table, filename, "search_params", config, "ProteinIdentificationArrowIO");
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
  auto col_run_id = ArrowIOHelpers::getColumn(table, SearchParamsSchema::RUN_IDENTIFIER, true, "ProteinIdentificationArrowIO");
  auto col_search_engine = ArrowIOHelpers::getColumn(table, SearchParamsSchema::SEARCH_ENGINE, true, "ProteinIdentificationArrowIO");
  auto col_se_version = ArrowIOHelpers::getColumn(table, SearchParamsSchema::SEARCH_ENGINE_VERSION, false, "ProteinIdentificationArrowIO");
  auto col_inf_engine = ArrowIOHelpers::getColumn(table, SearchParamsSchema::INFERENCE_ENGINE, false, "ProteinIdentificationArrowIO");
  auto col_inf_version = ArrowIOHelpers::getColumn(table, SearchParamsSchema::INFERENCE_ENGINE_VERSION, false, "ProteinIdentificationArrowIO");
  auto col_date = ArrowIOHelpers::getColumn(table, SearchParamsSchema::DATE, false, "ProteinIdentificationArrowIO");
  auto col_score_type = ArrowIOHelpers::getColumn(table, SearchParamsSchema::SCORE_TYPE, true, "ProteinIdentificationArrowIO");
  auto col_higher_better = ArrowIOHelpers::getColumn(table, SearchParamsSchema::HIGHER_SCORE_BETTER, true, "ProteinIdentificationArrowIO");
  auto col_sig_threshold = ArrowIOHelpers::getColumn(table, SearchParamsSchema::SIGNIFICANCE_THRESHOLD, false, "ProteinIdentificationArrowIO");
  auto col_db = ArrowIOHelpers::getColumn(table, SearchParamsSchema::DB, false, "ProteinIdentificationArrowIO");
  auto col_db_version = ArrowIOHelpers::getColumn(table, SearchParamsSchema::DB_VERSION, false, "ProteinIdentificationArrowIO");
  auto col_taxonomy = ArrowIOHelpers::getColumn(table, SearchParamsSchema::TAXONOMY, false, "ProteinIdentificationArrowIO");
  auto col_charges = ArrowIOHelpers::getColumn(table, SearchParamsSchema::CHARGES, false, "ProteinIdentificationArrowIO");
  auto col_mass_type = ArrowIOHelpers::getColumn(table, SearchParamsSchema::MASS_TYPE, true, "ProteinIdentificationArrowIO");
  auto col_precursor_tol = ArrowIOHelpers::getColumn(table, SearchParamsSchema::PRECURSOR_MASS_TOLERANCE, true, "ProteinIdentificationArrowIO");
  auto col_precursor_ppm = ArrowIOHelpers::getColumn(table, SearchParamsSchema::PRECURSOR_MASS_TOLERANCE_PPM, true, "ProteinIdentificationArrowIO");
  auto col_fragment_tol = ArrowIOHelpers::getColumn(table, SearchParamsSchema::FRAGMENT_MASS_TOLERANCE, true, "ProteinIdentificationArrowIO");
  auto col_fragment_ppm = ArrowIOHelpers::getColumn(table, SearchParamsSchema::FRAGMENT_MASS_TOLERANCE_PPM, true, "ProteinIdentificationArrowIO");
  auto col_enzyme = ArrowIOHelpers::getColumn(table, SearchParamsSchema::DIGESTION_ENZYME, false, "ProteinIdentificationArrowIO");
  auto col_enzyme_spec = ArrowIOHelpers::getColumn(table, SearchParamsSchema::ENZYME_TERM_SPECIFICITY, false, "ProteinIdentificationArrowIO");
  auto col_missed_cleavages = ArrowIOHelpers::getColumn(table, SearchParamsSchema::MISSED_CLEAVAGES, true, "ProteinIdentificationArrowIO");
  auto col_fixed_mods = ArrowIOHelpers::getColumn(table, SearchParamsSchema::FIXED_MODIFICATIONS, false, "ProteinIdentificationArrowIO");
  auto col_var_mods = ArrowIOHelpers::getColumn(table, SearchParamsSchema::VARIABLE_MODIFICATIONS, false, "ProteinIdentificationArrowIO");
  auto col_ms_run_paths = ArrowIOHelpers::getColumn(table, SearchParamsSchema::PRIMARY_MS_RUN_PATHS, false, "ProteinIdentificationArrowIO");
  auto col_metavalues = ArrowIOHelpers::getColumn(table, SearchParamsSchema::METAVALUES, false, "ProteinIdentificationArrowIO");
  auto col_sp_metavalues = ArrowIOHelpers::getColumn(table, SearchParamsSchema::SP_METAVALUES, false, "ProteinIdentificationArrowIO");

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
    prot_id.setIdentifier(ArrowIOHelpers::getStringValue(col_run_id, row));
    prot_id.setSearchEngine(ArrowIOHelpers::getStringValue(col_search_engine, row));
    prot_id.setSearchEngineVersion(ArrowIOHelpers::getStringValue(col_se_version, row));
    prot_id.setScoreType(ArrowIOHelpers::getStringValue(col_score_type, row));
    prot_id.setHigherScoreBetter(ArrowIOHelpers::getBoolValue(col_higher_better, row));
    prot_id.setSignificanceThreshold(ArrowIOHelpers::getDoubleValue(col_sig_threshold, row));

    // Date (timestamp; normalize to seconds based on the column's actual time unit).
    // The schema declares timestamp(SECOND), but Parquet may persist/restore the column
    // as nanoseconds, so we have to inspect the runtime unit rather than assume seconds.
    if (!ArrowIOHelpers::isNull(col_date, row))
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

    // Primary MS run paths (an empty list stores none, see ProteinIdentification::setPrimaryMSRunPath)
    if (col_ms_run_paths)
    {
      auto ms_run_paths = readStringList_(col_ms_run_paths, row);
      StringList sl(ms_run_paths.begin(), ms_run_paths.end());
      prot_id.setPrimaryMSRunPath(sl);
    }

    // SearchParameters
    ProteinIdentification::SearchParameters sp;
    sp.db = ArrowIOHelpers::getStringValue(col_db, row);
    sp.db_version = ArrowIOHelpers::getStringValue(col_db_version, row);
    sp.taxonomy = ArrowIOHelpers::getStringValue(col_taxonomy, row);
    sp.charges = ArrowIOHelpers::getStringValue(col_charges, row);

    std::string mass_type_str = ArrowIOHelpers::getStringValue(col_mass_type, row);
    sp.mass_type = (mass_type_str == "AVERAGE") ?
      ProteinIdentification::PeakMassType::AVERAGE :
      ProteinIdentification::PeakMassType::MONOISOTOPIC;

    sp.precursor_mass_tolerance = ArrowIOHelpers::getDoubleValue(col_precursor_tol, row);
    sp.precursor_mass_tolerance_ppm = ArrowIOHelpers::getBoolValue(col_precursor_ppm, row);
    sp.fragment_mass_tolerance = ArrowIOHelpers::getDoubleValue(col_fragment_tol, row);
    sp.fragment_mass_tolerance_ppm = ArrowIOHelpers::getBoolValue(col_fragment_ppm, row);
    sp.missed_cleavages = static_cast<UInt>(ArrowIOHelpers::getInt32Value(col_missed_cleavages, row));

    // Enzyme
    std::string enzyme_name = ArrowIOHelpers::getStringValue(col_enzyme, row);
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
    std::string spec_str = ArrowIOHelpers::getStringValue(col_enzyme_spec, row);
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
    ModificationDefinitionIO::registerFrom(sp); // before any peptidoform is parsed

    prot_id.setSearchParameters(sp);

    // Inference engine (must be set after setSearchParameters)
    std::string inf_engine = ArrowIOHelpers::getStringValue(col_inf_engine, row);
    if (!inf_engine.empty())
    {
      prot_id.setInferenceEngine(inf_engine);
      std::string inf_version = ArrowIOHelpers::getStringValue(col_inf_version, row);
      if (!inf_version.empty())
      {
        prot_id.setInferenceEngineVersion(inf_version);
      }
    }

    // ProteinIdentification metavalues (exclude keys handled via dedicated columns)
    // spectra_data has its own column; spectra_data_raw only travels as a meta value.
    static const std::unordered_set<std::string> excluded_prot_id_mvs = {
      "InferenceEngine", "InferenceEngineVersion",
      "spectra_data"
    };
    ArrowIOHelpers::readMetaValues(col_metavalues, row, prot_id, excluded_prot_id_mvs);
    prot_id.removeEmptyPrimaryMSRunPaths();

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
  auto col_accession = ArrowIOHelpers::getColumn(table, ProteinSchema::ACCESSION, true, "ProteinIdentificationArrowIO");
  auto col_score = ArrowIOHelpers::getColumn(table, ProteinSchema::SCORE, true, "ProteinIdentificationArrowIO");
  auto col_rank = ArrowIOHelpers::getColumn(table, ProteinSchema::RANK, false, "ProteinIdentificationArrowIO");
  auto col_coverage = ArrowIOHelpers::getColumn(table, ProteinSchema::COVERAGE, false, "ProteinIdentificationArrowIO");
  auto col_sequence = ArrowIOHelpers::getColumn(table, ProteinSchema::SEQUENCE, false, "ProteinIdentificationArrowIO");
  auto col_description = ArrowIOHelpers::getColumn(table, ProteinSchema::DESCRIPTION, false, "ProteinIdentificationArrowIO");
  auto col_is_decoy = ArrowIOHelpers::getColumn(table, ProteinSchema::IS_DECOY, false, "ProteinIdentificationArrowIO");
  auto col_run_id = ArrowIOHelpers::getColumn(table, ProteinSchema::RUN_IDENTIFIER, true, "ProteinIdentificationArrowIO");
  auto col_modifications = ArrowIOHelpers::getColumn(table, ProteinSchema::MODIFICATIONS, false, "ProteinIdentificationArrowIO");
  auto col_metavalues = ArrowIOHelpers::getColumn(table, ProteinSchema::METAVALUES, false, "ProteinIdentificationArrowIO");

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
    std::string run_id = ArrowIOHelpers::getStringValue(col_run_id, row);

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
    hit.setAccession(ArrowIOHelpers::getStringValue(col_accession, row));
    hit.setScore(ArrowIOHelpers::getDoubleValue(col_score, row));

    if (!ArrowIOHelpers::isNull(col_rank, row))
    {
      hit.setRank(static_cast<UInt>(ArrowIOHelpers::getInt32Value(col_rank, row)));
    }

    if (!ArrowIOHelpers::isNull(col_coverage, row))
    {
      hit.setCoverage(ArrowIOHelpers::getDoubleValue(col_coverage, row));
    }

    if (!ArrowIOHelpers::isNull(col_sequence, row))
    {
      hit.setSequence(ArrowIOHelpers::getStringValue(col_sequence, row));
    }

    if (!ArrowIOHelpers::isNull(col_description, row))
    {
      hit.setDescription(ArrowIOHelpers::getStringValue(col_description, row));
    }

    // is_decoy -> target_decoy metavalue
    if (!ArrowIOHelpers::isNull(col_is_decoy, row))
    {
      bool is_decoy = ArrowIOHelpers::getBoolValue(col_is_decoy, row);
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
          std::string mod_id = mod_name_arr->GetString(i);
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
  auto col_group_type = ArrowIOHelpers::getColumn(table, ProteinGroupSchema::GROUP_TYPE, true, "ProteinIdentificationArrowIO");
  auto col_probability = ArrowIOHelpers::getColumn(table, ProteinGroupSchema::PROBABILITY, true, "ProteinIdentificationArrowIO");
  auto col_accessions = ArrowIOHelpers::getColumn(table, ProteinGroupSchema::ACCESSIONS, true, "ProteinIdentificationArrowIO");
  auto col_run_id = ArrowIOHelpers::getColumn(table, ProteinGroupSchema::RUN_IDENTIFIER, true, "ProteinIdentificationArrowIO");
  auto col_group_index = ArrowIOHelpers::getColumn(table, ProteinGroupSchema::GROUP_INDEX, false, "ProteinIdentificationArrowIO"); // informational, not needed for reconstruction
  auto col_float_data = ArrowIOHelpers::getColumn(table, ProteinGroupSchema::FLOAT_DATA, false, "ProteinIdentificationArrowIO");
  auto col_string_data = ArrowIOHelpers::getColumn(table, ProteinGroupSchema::STRING_DATA, false, "ProteinIdentificationArrowIO");
  auto col_integer_data = ArrowIOHelpers::getColumn(table, ProteinGroupSchema::INTEGER_DATA, false, "ProteinIdentificationArrowIO");

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
    std::string run_id = ArrowIOHelpers::getStringValue(col_run_id, row);

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
    group.probability = ArrowIOHelpers::getDoubleValue(col_probability, row);
    auto accessions = readStringList_(col_accessions, row);
    group.accessions.assign(accessions.begin(), accessions.end());

    // Data arrays
    readFloatDataArrays(col_float_data, row, group);
    readStringDataArrays(col_string_data, row, group);
    readIntegerDataArrays(col_integer_data, row, group);

    // Insert into appropriate vector
    std::string group_type = ArrowIOHelpers::getStringValue(col_group_type, row);
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
  const std::string& filename,
  std::vector<ProteinIdentification>& protein_identifications)
{
  auto table = Internal::ArrowFileIO::readParquetTable(filename, "ProteinIdentificationArrowIO");
  if (!table) return false;
  return importSearchParamsFromArrow(table, protein_identifications);
}


bool ProteinIdentificationArrowIO::importProteinsFromParquet(
  const std::string& filename,
  std::vector<ProteinIdentification>& protein_identifications)
{
  auto table = Internal::ArrowFileIO::readParquetTable(filename, "ProteinIdentificationArrowIO");
  if (!table) return false;
  return importProteinsFromArrow(table, protein_identifications);
}


bool ProteinIdentificationArrowIO::importProteinGroupsFromParquet(
  const std::string& filename,
  std::vector<ProteinIdentification>& protein_identifications)
{
  auto table = Internal::ArrowFileIO::readParquetTable(filename, "ProteinIdentificationArrowIO");
  if (!table) return false;
  return importProteinGroupsFromArrow(table, protein_identifications);
}


bool ProteinIdentificationArrowIO::importFromParquet(
  const std::string& proteins_filename,
  const std::string& protein_groups_filename,
  const std::string& search_params_filename,
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

std::map<std::string, std::string> ProteinIdentificationArrowIO::synthesizeRunIdentifiers(
  std::vector<ProteinIdentification>& protein_identifications)
{
  std::map<std::string, std::string> rename;
  for (auto& prot_id : protein_identifications)
  {
    // Mirror IdXMLFile.cpp:530 — `<search_engine>_<date>_<UniqueIdGenerator>`.
    const std::string se = prot_id.getSearchEngine().empty() ? std::string("unknown") : prot_id.getSearchEngine();
    const std::string dt = prot_id.getDateTime().isValid()
                          ? prot_id.getDateTime().toString()
                          : std::string("1900-01-01T00:00:00");
    const std::string stored = prot_id.getIdentifier();
    const std::string synthesized = se + "_" + dt + "_" + StringUtils::toStr(UniqueIdGenerator::getUniqueId());

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
  const std::map<std::string, std::string>& rename,
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
  std::set<std::string> s;
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
