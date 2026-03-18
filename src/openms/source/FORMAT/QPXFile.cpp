// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/QPXFile.h>

#ifdef WITH_PARQUET

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/Scores.h>
#include <OpenMS/CHEMISTRY/ProForma.h>
#include <OpenMS/METADATA/SpectrumNativeIDParser.h>

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
#include <map>
#include <OpenMS/DATASTRUCTURES/DateTime.h>

namespace OpenMS
{

namespace // anonymous
{
  /// Try to extract a scan number from a native ID string.
  /// Uses SpectrumNativeIDParser to auto-detect the native ID format.
  /// Returns -1 on failure.
  Int extractScan(const String& native_id)
  {
    std::string regex_str = SpectrumNativeIDParser::getRegExFromNativeID(native_id);
    if (regex_str.empty())
    {
      return -1;
    }
    boost::regex scan_regex(regex_str);
    return SpectrumNativeIDParser::extractScanNumber(native_id, scan_regex, true);
  }
} // anonymous namespace


std::shared_ptr<arrow::Table> QPXFile::exportToArrow(
  const std::vector<ProteinIdentification>& protein_identifications,
  const PeptideIdentificationList& peptide_identifications,
  bool export_all_psms)
{
  // -- Simple column builders --
  arrow::StringBuilder sequence_builder, peptidoform_builder;
  arrow::StringBuilder reference_file_builder, score_type_builder;
  arrow::StringBuilder spectrum_reference_builder, cv_params_builder;
  arrow::StringBuilder run_identifier_builder;
  arrow::Int32Builder precursor_charge_builder, rank_builder, p_id_builder;
  arrow::BooleanBuilder is_decoy_builder, higher_score_better_builder;
  arrow::Int32Builder scan_builder;
  arrow::DoubleBuilder pep_builder, calculated_mz_builder, observed_mz_builder;
  arrow::DoubleBuilder rt_builder, ion_mobility_builder, predicted_rt_builder, score_builder;

  // -- protein_accessions list<utf8> --
  auto pa_vb = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder protein_accessions_builder(arrow::default_memory_pool(), pa_vb);

  // -- Modifications list<struct{name, accession, positions: list<struct{position, scores}>}> --
  auto pos_position_b = std::make_shared<arrow::StringBuilder>();
  auto pos_scores_b = std::make_shared<arrow::DoubleBuilder>();
  auto position_struct_type = arrow::struct_({
    arrow::field("position", arrow::utf8()),
    arrow::field("scores", arrow::float64())
  });
  auto position_struct_b = std::make_shared<arrow::StructBuilder>(
    position_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{pos_position_b, pos_scores_b});
  auto positions_list_b = std::make_shared<arrow::ListBuilder>(arrow::default_memory_pool(), position_struct_b);

  auto mod_name_b = std::make_shared<arrow::StringBuilder>();
  auto mod_acc_b = std::make_shared<arrow::StringBuilder>();
  auto modification_struct_type = arrow::struct_({
    arrow::field("name", arrow::utf8()),
    arrow::field("accession", arrow::utf8()),
    arrow::field("positions", arrow::list(position_struct_type))
  });
  auto modification_struct_b = std::make_shared<arrow::StructBuilder>(
    modification_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{mod_name_b, mod_acc_b, positions_list_b});
  arrow::ListBuilder modifications_builder(arrow::default_memory_pool(), modification_struct_b);

  // -- additional_scores list<struct{score_name, score_value, higher_better}> --
  auto as_name_b = std::make_shared<arrow::StringBuilder>();
  auto as_value_b = std::make_shared<arrow::DoubleBuilder>();
  auto as_hb_b = std::make_shared<arrow::BooleanBuilder>();
  auto as_struct_type = arrow::struct_({
    arrow::field("score_name", arrow::utf8()),
    arrow::field("score_value", arrow::float64()),
    arrow::field("higher_better", arrow::boolean())
  });
  auto as_struct_b = std::make_shared<arrow::StructBuilder>(
    as_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{as_name_b, as_value_b, as_hb_b});
  arrow::ListBuilder additional_scores_builder(arrow::default_memory_pool(), as_struct_b);

  // -- psm_metavalues list<struct{name, value, value_type}> --
  auto pmv_name_b = std::make_shared<arrow::StringBuilder>();
  auto pmv_value_b = std::make_shared<arrow::StringBuilder>();
  auto pmv_type_b = std::make_shared<arrow::StringBuilder>();
  auto pmv_struct_type = arrow::struct_({
    arrow::field("name", arrow::utf8()),
    arrow::field("value", arrow::utf8()),
    arrow::field("value_type", arrow::utf8())
  });
  auto pmv_struct_b = std::make_shared<arrow::StructBuilder>(
    pmv_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{pmv_name_b, pmv_value_b, pmv_type_b});
  arrow::ListBuilder psm_metavalues_builder(arrow::default_memory_pool(), pmv_struct_b);

  // -- spectrum_metavalues list<struct{name, value, value_type}> --
  auto smv_name_b = std::make_shared<arrow::StringBuilder>();
  auto smv_value_b = std::make_shared<arrow::StringBuilder>();
  auto smv_type_b = std::make_shared<arrow::StringBuilder>();
  auto smv_struct_type = arrow::struct_({
    arrow::field("name", arrow::utf8()),
    arrow::field("value", arrow::utf8()),
    arrow::field("value_type", arrow::utf8())
  });
  auto smv_struct_b = std::make_shared<arrow::StructBuilder>(
    smv_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{smv_name_b, smv_value_b, smv_type_b});
  arrow::ListBuilder spectrum_metavalues_builder(arrow::default_memory_pool(), smv_struct_b);

  // Estimate total rows for capacity reservation
  Size num_rows = 0;
  for (const auto& pep_id : peptide_identifications)
  {
    if (pep_id.getHits().empty()) continue;
    num_rows += export_all_psms ? pep_id.getHits().size() : 1;
  }

  // Reserve capacity for simple builders
  arrow::Status status;

  status = sequence_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: sequence_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = peptidoform_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: peptidoform_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = precursor_charge_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: precursor_charge_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = pep_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: pep_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = is_decoy_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: is_decoy_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = higher_score_better_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: higher_score_better_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = calculated_mz_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: calculated_mz_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = observed_mz_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: observed_mz_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = predicted_rt_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: predicted_rt_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = reference_file_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: reference_file_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = cv_params_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: cv_params_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = scan_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: scan_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = rt_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: rt_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = ion_mobility_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: ion_mobility_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = spectrum_reference_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: spectrum_reference_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = score_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: score_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = score_type_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: score_type_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = rank_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: rank_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = p_id_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: p_id_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = run_identifier_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: run_identifier_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }

  // Build file name lookup from ProteinIdentification primary MS run paths
  std::map<String, String> id_to_filename;
  for (const auto& prot_id : protein_identifications)
  {
    StringList ms_runs;
    prot_id.getPrimaryMSRunPath(ms_runs);
    if (!ms_runs.empty())
    {
      id_to_filename[prot_id.getIdentifier()] = ms_runs[0];
    }
  }

  // Metavalue keys excluded from psm_metavalues (they have dedicated columns)
  static const std::unordered_set<std::string> excluded_hit_mvs = {
    "target_decoy", "predicted_RT", "predicted_rt", "ion_mobility", "IM",
    "scan", "reference_file_name"
  };

  IDScoreSwitcherAlgorithm idsa;

  Int p_id_index = 0;
  for (const auto& pep_id : peptide_identifications)
  {
    const auto& hits = pep_id.getHits();
    if (hits.empty())
    {
      ++p_id_index;
      continue;
    }

    Size num_hits = export_all_psms ? hits.size() : 1;

    // Get shared per-spectrum values
    std::string spec_ref = pep_id.getSpectrumReference();

    // Detect PEP score type for this identification
    IDScoreSwitcherAlgorithm::ScoreSearchResult pep_result;
    try
    {
      pep_result = idsa.findScoreType(pep_id, IDScoreSwitcherAlgorithm::ScoreType::PEP);
    }
    catch (...) {} // No PEP score available

    // Lookup reference file name
    auto fn_it = id_to_filename.find(pep_id.getIdentifier());

    for (Size hit_idx = 0; hit_idx < num_hits; ++hit_idx)
    {
      const PeptideHit& hit = hits[hit_idx];
      const auto& seq = hit.getSequence();

      // === sequence ===
      (void)sequence_builder.Append(seq.toUnmodifiedString());

      // === peptidoform (ProForma canonical) ===
      auto pf = ProForma::fromAASequence(seq);
      (void)peptidoform_builder.Append(ProForma::toString(pf, ProForma::WriteMode::CANONICAL));

      // === modifications (structured) ===
      (void)modifications_builder.Append();
      if (seq.isModified())
      {
        // Collect modifications keyed by name to group positions
        std::map<std::string, std::pair<std::string, std::vector<std::string>>> mod_map;

        // N-terminal
        if (seq.hasNTerminalModification())
        {
          const auto* mod = seq.getNTerminalModification();
          std::string name = seq.getNTerminalModificationName();
          std::string acc_str;
          if (mod && mod->getUniModRecordId() > 0)
          {
            acc_str = "UNIMOD:" + std::to_string(mod->getUniModRecordId());
          }
          mod_map[name] = {acc_str, {}};
          mod_map[name].second.push_back("N-term.0");
        }

        // Residue modifications
        for (Size pos = 0; pos < seq.size(); ++pos)
        {
          const auto& residue = seq[pos];
          if (residue.isModified())
          {
            const auto* mod = residue.getModification();
            std::string name = residue.getModificationName();
            std::string acc_str;
            if (mod && mod->getUniModRecordId() > 0)
            {
              acc_str = "UNIMOD:" + std::to_string(mod->getUniModRecordId());
            }
            std::string pos_str = std::string(residue.getOneLetterCode()) + "." + std::to_string(pos + 1);
            if (mod_map.find(name) == mod_map.end())
            {
              mod_map[name] = {acc_str, {}};
            }
            mod_map[name].second.push_back(pos_str);
          }
        }

        // C-terminal
        if (seq.hasCTerminalModification())
        {
          const auto* mod = seq.getCTerminalModification();
          std::string name = seq.getCTerminalModificationName();
          std::string acc_str;
          if (mod && mod->getUniModRecordId() > 0)
          {
            acc_str = "UNIMOD:" + std::to_string(mod->getUniModRecordId());
          }
          std::string pos_str = "C-term." + std::to_string(seq.size() + 1);
          if (mod_map.find(name) == mod_map.end())
          {
            mod_map[name] = {acc_str, {}};
          }
          mod_map[name].second.push_back(pos_str);
        }

        // Write modification structs
        for (const auto& [name, acc_positions] : mod_map)
        {
          (void)modification_struct_b->Append();
          (void)mod_name_b->Append(name);
          if (acc_positions.first.empty())
          {
            (void)mod_acc_b->AppendNull();
          }
          else
          {
            (void)mod_acc_b->Append(acc_positions.first);
          }
          (void)positions_list_b->Append();
          for (const auto& pos_str : acc_positions.second)
          {
            (void)position_struct_b->Append();
            (void)pos_position_b->Append(pos_str);
            (void)pos_scores_b->AppendNull();
          }
        }
      }
      // else: empty modifications list (Append() already called)

      // === precursor_charge ===
      int charge = hit.getCharge();
      (void)precursor_charge_builder.Append(charge);

      // === posterior_error_probability ===
      if (pep_result.is_main_score_type)
      {
        (void)pep_builder.Append(hit.getScore());
      }
      else if (!pep_result.score_name.empty() && hit.metaValueExists(pep_result.score_name))
      {
        (void)pep_builder.Append(static_cast<double>(hit.getMetaValue(pep_result.score_name)));
      }
      else
      {
        (void)pep_builder.AppendNull();
      }

      // === is_decoy (nullable) ===
      if (hit.metaValueExists("target_decoy"))
      {
        std::string td = hit.getMetaValue("target_decoy").toString();
        (void)is_decoy_builder.Append(td.substr(0, 5) == "decoy");
      }
      else
      {
        (void)is_decoy_builder.AppendNull();
      }

      // === calculated_mz ===
      if (charge > 0)
      {
        (void)calculated_mz_builder.Append(seq.getMZ(charge));
      }
      else
      {
        (void)calculated_mz_builder.AppendNull();
      }

      // === observed_mz ===
      {
        double obs_mz = pep_id.getMZ();
        if (obs_mz == obs_mz) // not NaN
        {
          (void)observed_mz_builder.Append(obs_mz);
        }
        else
        {
          (void)observed_mz_builder.AppendNull();
        }
      }

      // === additional_scores from hit metavalues ===
      (void)additional_scores_builder.Append();
      {
        std::vector<String> keys;
        hit.getKeys(keys);
        for (const auto& key : keys)
        {
          if (excluded_hit_mvs.count(key)) continue;
          const DataValue& val = hit.getMetaValue(key);
          if ((val.valueType() == DataValue::INT_VALUE || val.valueType() == DataValue::DOUBLE_VALUE)
              && Scores::isKnownScoreType(key))
          {
            (void)as_struct_b->Append();
            (void)as_name_b->Append(key);
            (void)as_value_b->Append(static_cast<double>(val));
            try
            {
              auto st = IDScoreSwitcherAlgorithm::toScoreTypeEnum(key);
              (void)as_hb_b->Append(idsa.isScoreTypeHigherBetter(st));
            }
            catch (...)
            {
              (void)as_hb_b->AppendNull();
            }
          }
        }
      }

      // === protein_accessions ===
      (void)protein_accessions_builder.Append();
      for (const auto& ev : hit.getPeptideEvidences())
      {
        (void)pa_vb->Append(ev.getProteinAccession());
      }

      // === predicted_rt ===
      if (hit.metaValueExists("predicted_RT"))
      {
        (void)predicted_rt_builder.Append(static_cast<double>(hit.getMetaValue("predicted_RT")));
      }
      else if (hit.metaValueExists("predicted_rt"))
      {
        (void)predicted_rt_builder.Append(static_cast<double>(hit.getMetaValue("predicted_rt")));
      }
      else
      {
        (void)predicted_rt_builder.AppendNull();
      }

      // === reference_file_name ===
      if (fn_it != id_to_filename.end())
      {
        (void)reference_file_builder.Append(fn_it->second);
      }
      else
      {
        (void)reference_file_builder.AppendNull();
      }

      // === cv_params (reserved) ===
      (void)cv_params_builder.AppendNull();

      // === scan / spectrum_reference ===
      if (spec_ref.empty())
      {
        (void)spectrum_reference_builder.AppendNull();
        (void)scan_builder.AppendNull();
      }
      else
      {
        (void)spectrum_reference_builder.Append(spec_ref);
        Int scan_num = extractScan(spec_ref);
        if (scan_num >= 0)
        {
          (void)scan_builder.Append(scan_num);
        }
        else
        {
          (void)scan_builder.AppendNull();
        }
      }

      // === rt ===
      double rt = pep_id.getRT();
      if (rt == rt) // not NaN
      {
        (void)rt_builder.Append(rt);
      }
      else
      {
        (void)rt_builder.AppendNull();
      }

      // === ion_mobility ===
      if (hit.metaValueExists("ion_mobility"))
      {
        (void)ion_mobility_builder.Append(static_cast<double>(hit.getMetaValue("ion_mobility")));
      }
      else if (hit.metaValueExists("IM"))
      {
        (void)ion_mobility_builder.Append(static_cast<double>(hit.getMetaValue("IM")));
      }
      else if (pep_id.metaValueExists("ion_mobility"))
      {
        (void)ion_mobility_builder.Append(static_cast<double>(pep_id.getMetaValue("ion_mobility")));
      }
      else if (pep_id.metaValueExists("IM"))
      {
        (void)ion_mobility_builder.Append(static_cast<double>(pep_id.getMetaValue("IM")));
      }
      else
      {
        (void)ion_mobility_builder.AppendNull();
      }

      // === score / score_type / higher_score_better ===
      (void)score_builder.Append(hit.getScore());
      (void)score_type_builder.Append(pep_id.getScoreType());
      (void)higher_score_better_builder.Append(pep_id.isHigherScoreBetter());

      // === rank (0-based) ===
      (void)rank_builder.Append(static_cast<int32_t>(hit_idx));

      // === P_ID (parent spectrum index) ===
      (void)p_id_builder.Append(p_id_index);

      // === run_identifier ===
      (void)run_identifier_builder.Append(pep_id.getIdentifier());

      // === psm_metavalues ===
      (void)psm_metavalues_builder.Append();
      {
        std::vector<String> keys;
        hit.getKeys(keys);
        for (const auto& key : keys)
        {
          if (excluded_hit_mvs.count(key)) continue;
          const DataValue& val = hit.getMetaValue(key);
          // Skip known scores (already in additional_scores)
          if ((val.valueType() == DataValue::INT_VALUE || val.valueType() == DataValue::DOUBLE_VALUE)
              && Scores::isKnownScoreType(key))
          {
            continue;
          }
          (void)pmv_struct_b->Append();
          (void)pmv_name_b->Append(key);
          (void)pmv_value_b->Append(val.toString());

          switch (val.valueType())
          {
            case DataValue::INT_VALUE: (void)pmv_type_b->Append("int"); break;
            case DataValue::DOUBLE_VALUE: (void)pmv_type_b->Append("double"); break;
            case DataValue::STRING_VALUE: (void)pmv_type_b->Append("string"); break;
            case DataValue::INT_LIST: (void)pmv_type_b->Append("int_list"); break;
            case DataValue::DOUBLE_LIST: (void)pmv_type_b->Append("double_list"); break;
            case DataValue::STRING_LIST: (void)pmv_type_b->Append("string_list"); break;
            default: (void)pmv_type_b->Append("string"); break;
          }
        }
      }

      // === spectrum_metavalues ===
      (void)spectrum_metavalues_builder.Append();
      {
        std::vector<String> keys;
        pep_id.getKeys(keys);
        for (const auto& key : keys)
        {
          // Skip keys that have dedicated columns
          if (key == "spectrum_reference" || key == "ion_mobility" || key == "IM") continue;
          const DataValue& val = pep_id.getMetaValue(key);
          (void)smv_struct_b->Append();
          (void)smv_name_b->Append(key);
          (void)smv_value_b->Append(val.toString());

          switch (val.valueType())
          {
            case DataValue::INT_VALUE: (void)smv_type_b->Append("int"); break;
            case DataValue::DOUBLE_VALUE: (void)smv_type_b->Append("double"); break;
            case DataValue::STRING_VALUE: (void)smv_type_b->Append("string"); break;
            case DataValue::INT_LIST: (void)smv_type_b->Append("int_list"); break;
            case DataValue::DOUBLE_LIST: (void)smv_type_b->Append("double_list"); break;
            case DataValue::STRING_LIST: (void)smv_type_b->Append("string_list"); break;
            default: (void)smv_type_b->Append("string"); break;
          }
        }
      }
    } // end hit loop

    ++p_id_index;
  } // end peptide identification loop

  // Finalize all arrays
  std::shared_ptr<arrow::Array> arr_sequence, arr_peptidoform, arr_modifications;
  std::shared_ptr<arrow::Array> arr_charge, arr_pep, arr_is_decoy;
  std::shared_ptr<arrow::Array> arr_calc_mz, arr_obs_mz, arr_additional_scores;
  std::shared_ptr<arrow::Array> arr_protein_acc, arr_predicted_rt, arr_ref_file;
  std::shared_ptr<arrow::Array> arr_cv_params, arr_scan, arr_rt, arr_ion_mobility;
  std::shared_ptr<arrow::Array> arr_spectrum_ref, arr_score, arr_score_type;
  std::shared_ptr<arrow::Array> arr_higher_score_better;
  std::shared_ptr<arrow::Array> arr_rank, arr_p_id;
  std::shared_ptr<arrow::Array> arr_psm_mvs, arr_spectrum_mvs;

  status = sequence_builder.Finish(&arr_sequence);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: sequence_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = peptidoform_builder.Finish(&arr_peptidoform);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: peptidoform_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = modifications_builder.Finish(&arr_modifications);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: modifications_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = precursor_charge_builder.Finish(&arr_charge);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: precursor_charge_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = pep_builder.Finish(&arr_pep);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: pep_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = is_decoy_builder.Finish(&arr_is_decoy);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: is_decoy_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = calculated_mz_builder.Finish(&arr_calc_mz);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: calculated_mz_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = observed_mz_builder.Finish(&arr_obs_mz);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: observed_mz_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = additional_scores_builder.Finish(&arr_additional_scores);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: additional_scores_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = protein_accessions_builder.Finish(&arr_protein_acc);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: protein_accessions_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = predicted_rt_builder.Finish(&arr_predicted_rt);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: predicted_rt_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = reference_file_builder.Finish(&arr_ref_file);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: reference_file_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = cv_params_builder.Finish(&arr_cv_params);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: cv_params_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = scan_builder.Finish(&arr_scan);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: scan_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = rt_builder.Finish(&arr_rt);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: rt_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = ion_mobility_builder.Finish(&arr_ion_mobility);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: ion_mobility_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = spectrum_reference_builder.Finish(&arr_spectrum_ref);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: spectrum_reference_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = score_builder.Finish(&arr_score);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: score_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = score_type_builder.Finish(&arr_score_type);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: score_type_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = higher_score_better_builder.Finish(&arr_higher_score_better);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: higher_score_better_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = rank_builder.Finish(&arr_rank);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: rank_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = p_id_builder.Finish(&arr_p_id);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: p_id_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = psm_metavalues_builder.Finish(&arr_psm_mvs);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: psm_metavalues_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = spectrum_metavalues_builder.Finish(&arr_spectrum_mvs);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: spectrum_metavalues_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }

  std::shared_ptr<arrow::Array> arr_run_identifier;
  status = run_identifier_builder.Finish(&arr_run_identifier);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: run_identifier_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }

  // Build schema matching the Python to_psm_arrow() column order
  auto schema = arrow::schema({
    arrow::field("sequence", arrow::utf8()),
    arrow::field("peptidoform", arrow::utf8()),
    arrow::field("modifications", modifications_builder.type()),
    arrow::field("precursor_charge", arrow::int32()),
    arrow::field("posterior_error_probability", arrow::float64()),
    arrow::field("is_decoy", arrow::boolean()),
    arrow::field("calculated_mz", arrow::float64()),
    arrow::field("observed_mz", arrow::float64()),
    arrow::field("additional_scores", additional_scores_builder.type()),
    arrow::field("protein_accessions", arrow::list(arrow::utf8())),
    arrow::field("predicted_rt", arrow::float64()),
    arrow::field("reference_file_name", arrow::utf8()),
    arrow::field("cv_params", arrow::utf8()),
    arrow::field("scan", arrow::int32()),
    arrow::field("rt", arrow::float64()),
    arrow::field("ion_mobility", arrow::float64()),
    // OpenMS-specific columns
    arrow::field("spectrum_reference", arrow::utf8()),
    arrow::field("score", arrow::float64()),
    arrow::field("score_type", arrow::utf8()),
    arrow::field("higher_score_better", arrow::boolean()),
    arrow::field("rank", arrow::int32()),
    arrow::field("peptide_identification_index", arrow::int32()),
    arrow::field("psm_metavalues", psm_metavalues_builder.type()),
    arrow::field("spectrum_metavalues", spectrum_metavalues_builder.type()),
    arrow::field("run_identifier", arrow::utf8()),
  });

  auto table = arrow::Table::Make(schema, {
    arr_sequence, arr_peptidoform, arr_modifications,
    arr_charge, arr_pep, arr_is_decoy,
    arr_calc_mz, arr_obs_mz, arr_additional_scores,
    arr_protein_acc, arr_predicted_rt, arr_ref_file,
    arr_cv_params, arr_scan, arr_rt, arr_ion_mobility,
    arr_spectrum_ref, arr_score, arr_score_type, arr_higher_score_better,
    arr_rank, arr_p_id, arr_psm_mvs, arr_spectrum_mvs,
    arr_run_identifier
  });

  return table;
}


bool QPXFile::exportToParquet(
  const std::vector<ProteinIdentification>& protein_identifications,
  const PeptideIdentificationList& peptide_identifications,
  const String& filename,
  bool export_all_psms,
  const ParquetWriteConfig& config)
{
  auto table = exportToArrow(protein_identifications, peptide_identifications, export_all_psms);
  if (!table)
  {
    OPENMS_LOG_ERROR << "QPXFile: Failed to create Arrow table" << std::endl;
    return false;
  }

  // Add QPX file metadata to the table schema (matches Python to_psm_qpx())
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
    std::string uuid_str(buf);
    
    auto metadata = arrow::key_value_metadata({
      {"qpx_version", "1.0"},
      {"creator", "OpenMS"},
      {"file_type", "psm"},
      {"creation_date", DateTime::nowUTC().toString("yyyy-MM-ddThh:mm:ssZ")},
      {"uuid", uuid_str},
      {"scan_format", "scan"},
      {"software_provider", "OpenMS"}
    });
    table = table->ReplaceSchemaMetadata(metadata);
  }

  // Open output file
  auto result = arrow::io::FileOutputStream::Open(filename);
  if (!result.ok())
  {
    OPENMS_LOG_ERROR << "QPXFile: Failed to open file: " << filename << std::endl;
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
    OPENMS_LOG_ERROR << "QPXFile: Failed to write Parquet: "
                     << write_status.ToString() << std::endl;
    return false;
  }

  return true;
}

} // namespace OpenMS

#endif // WITH_PARQUET
