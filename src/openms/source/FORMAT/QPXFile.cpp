// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/QPXFile.h>

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/FORMAT/ArrowIOHelpers.h>
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

  // -- fragment annotation arrays: mz_array, intensity_array, charge_array, ion_type_array --
  auto psm_mz_arr_vb = std::make_shared<arrow::FloatBuilder>();
  arrow::ListBuilder psm_mz_array_builder(arrow::default_memory_pool(), psm_mz_arr_vb);
  auto psm_int_arr_vb = std::make_shared<arrow::FloatBuilder>();
  arrow::ListBuilder psm_intensity_array_builder(arrow::default_memory_pool(), psm_int_arr_vb);
  auto psm_chg_arr_vb = std::make_shared<arrow::Int32Builder>();
  arrow::ListBuilder psm_charge_array_builder(arrow::default_memory_pool(), psm_chg_arr_vb);
  auto psm_ion_arr_vb = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder psm_ion_type_array_builder(arrow::default_memory_pool(), psm_ion_arr_vb);

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
  static const std::unordered_set<std::string> excluded_hit_mvs_psm = {
    "target_decoy", "predicted_RT", "predicted_rt", "ion_mobility", "IM",
    "scan", "reference_file_name",
    Constants::UserParam::FRAGMENT_ANNOTATION_USERPARAM  // dedicated mz/intensity/charge/ion_type arrays
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
          const std::string& name = seq.getNTerminalModificationName();
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
            const std::string& name = residue.getModificationName();
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
          const std::string& name = seq.getCTerminalModificationName();
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
          if (excluded_hit_mvs_psm.count(key)) continue;
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
          if (excluded_hit_mvs_psm.count(key)) continue;
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

      // === fragment annotation arrays (from PeakAnnotations) ===
      const auto& psm_peak_annotations = hit.getPeakAnnotations();
      if (!psm_peak_annotations.empty())
      {
        (void)psm_mz_array_builder.Append();
        (void)psm_intensity_array_builder.Append();
        (void)psm_charge_array_builder.Append();
        (void)psm_ion_type_array_builder.Append();
        for (const auto& pa : psm_peak_annotations)
        {
          (void)psm_mz_arr_vb->Append(static_cast<float>(pa.mz));
          (void)psm_int_arr_vb->Append(static_cast<float>(pa.intensity));
          (void)psm_chg_arr_vb->Append(pa.charge);
          (void)psm_ion_arr_vb->Append(pa.annotation);
        }
      }
      else
      {
        (void)psm_mz_array_builder.AppendNull();
        (void)psm_intensity_array_builder.AppendNull();
        (void)psm_charge_array_builder.AppendNull();
        (void)psm_ion_type_array_builder.AppendNull();
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
  std::shared_ptr<arrow::Array> arr_psm_mz_array, arr_psm_intensity_array, arr_psm_charge_array, arr_psm_ion_type_array;
  status = psm_mz_array_builder.Finish(&arr_psm_mz_array);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: psm_mz_array_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = psm_intensity_array_builder.Finish(&arr_psm_intensity_array);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: psm_intensity_array_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = psm_charge_array_builder.Finish(&arr_psm_charge_array);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: psm_charge_array_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = psm_ion_type_array_builder.Finish(&arr_psm_ion_type_array);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: psm_ion_type_array_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }

  // Build schema from registry
  auto schema = PSMSchema::schema();

  auto table = arrow::Table::Make(schema, {
    arr_sequence, arr_peptidoform, arr_modifications,
    arr_charge, arr_pep, arr_is_decoy,
    arr_calc_mz, arr_obs_mz, arr_additional_scores,
    arr_protein_acc, arr_predicted_rt, arr_ref_file,
    arr_cv_params, arr_scan, arr_rt, arr_ion_mobility,
    arr_spectrum_ref, arr_score, arr_score_type, arr_higher_score_better,
    arr_rank, arr_p_id, arr_psm_mvs, arr_spectrum_mvs,
    arr_run_identifier,
    arr_psm_mz_array, arr_psm_intensity_array, arr_psm_charge_array, arr_psm_ion_type_array
  });

  // Validate table against registry schema (strict — write path must match exactly)
  auto validation = ArrowSchemaValidation::validate(table, PSMSchema::schema());
  if (!validation.valid)
  {
    OPENMS_LOG_ERROR << "QPXFile: Schema validation failed: " << validation.toString() << "\n";
    return nullptr;
  }

  return table;
}

std::shared_ptr<arrow::Table> QPXFile::exportPSMsToQPXArrow(
  const std::vector<ProteinIdentification>& protein_identifications,
  const PeptideIdentificationList& peptide_identifications,
  bool export_all_psms)
{
  // -- Simple column builders --
  arrow::StringBuilder sequence_builder, peptidoform_builder;
  arrow::StringBuilder run_file_name_builder;
  arrow::Int16Builder charge_builder, missed_cleavages_builder;
  arrow::BooleanBuilder is_decoy_builder;
  arrow::FloatBuilder calculated_mz_builder, observed_mz_builder;
  arrow::FloatBuilder mass_error_ppm_builder;
  arrow::FloatBuilder rt_builder, ion_mobility_builder, predicted_rt_builder;
  arrow::DoubleBuilder pep_builder;

  // -- protein_accessions list<utf8> --
  auto pa_vb = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder protein_accessions_builder(arrow::default_memory_pool(), pa_vb);

  // -- Modifications list<struct{name, accession, positions: list<struct{position(int32), amino_acid, scores: list<struct>}>}> --
  auto score_name_b = std::make_shared<arrow::StringBuilder>();
  auto score_value_b = std::make_shared<arrow::DoubleBuilder>();
  auto score_hb_b = std::make_shared<arrow::BooleanBuilder>();
  auto scores_struct_type = arrow::struct_({
    arrow::field("score_name", arrow::utf8(), /*nullable=*/false),
    arrow::field("score_value", arrow::float64(), /*nullable=*/false),
    arrow::field("higher_better", arrow::boolean())
  });
  auto scores_struct_b = std::make_shared<arrow::StructBuilder>(
    scores_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{score_name_b, score_value_b, score_hb_b});
  auto scores_list_b = std::make_shared<arrow::ListBuilder>(arrow::default_memory_pool(), scores_struct_b);

  auto pos_position_b = std::make_shared<arrow::Int32Builder>();
  auto pos_amino_acid_b = std::make_shared<arrow::StringBuilder>();
  auto position_struct_type = arrow::struct_({
    arrow::field("position", arrow::int32(), /*nullable=*/false),
    arrow::field("amino_acid", arrow::utf8()),
    arrow::field("scores", arrow::list(scores_struct_type))
  });
  auto position_struct_b = std::make_shared<arrow::StructBuilder>(
    position_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{pos_position_b, pos_amino_acid_b, scores_list_b});
  auto positions_list_b = std::make_shared<arrow::ListBuilder>(arrow::default_memory_pool(), position_struct_b);

  auto mod_name_b = std::make_shared<arrow::StringBuilder>();
  auto mod_acc_b = std::make_shared<arrow::StringBuilder>();
  auto modification_struct_type = arrow::struct_({
    arrow::field("name", arrow::utf8(), /*nullable=*/false),
    arrow::field("accession", arrow::utf8()),
    arrow::field("positions", arrow::list(position_struct_type), /*nullable=*/false)
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
    arrow::field("score_name", arrow::utf8(), /*nullable=*/false),
    arrow::field("score_value", arrow::float64(), /*nullable=*/false),
    arrow::field("higher_better", arrow::boolean())
  });
  auto as_struct_b = std::make_shared<arrow::StructBuilder>(
    as_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{as_name_b, as_value_b, as_hb_b});
  arrow::ListBuilder additional_scores_builder(arrow::default_memory_pool(), as_struct_b);

  // -- cv_params list<struct{cv_name, cv_value}> --
  auto cv_name_b = std::make_shared<arrow::StringBuilder>();
  auto cv_value_b = std::make_shared<arrow::StringBuilder>();
  auto cv_struct_type = arrow::struct_({
    arrow::field("cv_name", arrow::utf8(), /*nullable=*/false),
    arrow::field("cv_value", arrow::utf8(), /*nullable=*/false)
  });
  auto cv_struct_b = std::make_shared<arrow::StructBuilder>(
    cv_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{cv_name_b, cv_value_b});
  arrow::ListBuilder cv_params_builder(arrow::default_memory_pool(), cv_struct_b);

  // -- scan list<int32> --
  auto scan_value_b = std::make_shared<arrow::Int32Builder>();
  arrow::ListBuilder scan_builder(arrow::default_memory_pool(), scan_value_b);

  // -- cross_links (null for now) --
  auto cross_links_type = QPXPSMSchema::crossLinksType();
  // We build a null array later

  // -- mz_array list<float32>, intensity_array list<float32>, charge_array list<int32>, ion_type_array list<utf8> --
  // Populated from PeptideHit::getPeakAnnotations() (fragment ion annotations).
  auto mz_arr_vb = std::make_shared<arrow::FloatBuilder>();
  arrow::ListBuilder mz_array_builder(arrow::default_memory_pool(), mz_arr_vb);
  auto int_arr_vb = std::make_shared<arrow::FloatBuilder>();
  arrow::ListBuilder intensity_array_builder(arrow::default_memory_pool(), int_arr_vb);
  auto chg_arr_vb = std::make_shared<arrow::Int32Builder>();
  arrow::ListBuilder charge_array_builder(arrow::default_memory_pool(), chg_arr_vb);
  auto ion_arr_vb = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder ion_type_array_builder(arrow::default_memory_pool(), ion_arr_vb);
  // -- ion_mobility_array (null for now — no per-fragment IM data) --

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
  status = charge_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: charge_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = pep_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: pep_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = is_decoy_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: is_decoy_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = calculated_mz_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: calculated_mz_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = observed_mz_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: observed_mz_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = mass_error_ppm_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: mass_error_ppm_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = predicted_rt_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: predicted_rt_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = run_file_name_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: run_file_name_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = rt_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: rt_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = ion_mobility_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: ion_mobility_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = missed_cleavages_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: missed_cleavages_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }

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

  // Metavalue keys excluded from additional_scores (they have dedicated columns)
  static const std::unordered_set<std::string> excluded_hit_mvs = {
    "target_decoy", "predicted_RT", "predicted_rt", "ion_mobility", "IM",
    "scan", "reference_file_name",
    Constants::UserParam::FRAGMENT_ANNOTATION_USERPARAM  // dedicated mz/intensity/charge/ion_type arrays
  };

  IDScoreSwitcherAlgorithm idsa;

  for (const auto& pep_id : peptide_identifications)
  {
    const auto& hits = pep_id.getHits();
    if (hits.empty())
    {
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

      // === sequence (non-nullable) ===
      (void)sequence_builder.Append(seq.toUnmodifiedString());

      // === peptidoform (ProForma canonical, non-nullable) ===
      auto pf = ProForma::fromAASequence(seq);
      (void)peptidoform_builder.Append(ProForma::toString(pf, ProForma::WriteMode::CANONICAL));

      // === modifications (structured with int32 position and amino_acid) ===
      (void)modifications_builder.Append();
      if (seq.isModified())
      {
        // Collect modifications keyed by name to group positions
        // value: {accession, vector<{position_int, amino_acid_char}>}
        struct PosInfo { int32_t position; std::string amino_acid; };
        std::map<std::string, std::pair<std::string, std::vector<PosInfo>>> mod_map;

        // N-terminal
        if (seq.hasNTerminalModification())
        {
          const auto* mod = seq.getNTerminalModification();
          const std::string& name = seq.getNTerminalModificationName();
          std::string acc_str;
          if (mod && mod->getUniModRecordId() > 0)
          {
            acc_str = "UNIMOD:" + std::to_string(mod->getUniModRecordId());
          }
          mod_map[name] = {acc_str, {}};
          mod_map[name].second.push_back({0, ""}); // position 0 for N-term
        }

        // Residue modifications
        for (Size pos = 0; pos < seq.size(); ++pos)
        {
          const auto& residue = seq[pos];
          if (residue.isModified())
          {
            const auto* mod = residue.getModification();
            const std::string& name = residue.getModificationName();
            std::string acc_str;
            if (mod && mod->getUniModRecordId() > 0)
            {
              acc_str = "UNIMOD:" + std::to_string(mod->getUniModRecordId());
            }
            if (mod_map.find(name) == mod_map.end())
            {
              mod_map[name] = {acc_str, {}};
            }
            mod_map[name].second.push_back({static_cast<int32_t>(pos + 1), std::string(residue.getOneLetterCode())});
          }
        }

        // C-terminal
        if (seq.hasCTerminalModification())
        {
          const auto* mod = seq.getCTerminalModification();
          const std::string& name = seq.getCTerminalModificationName();
          std::string acc_str;
          if (mod && mod->getUniModRecordId() > 0)
          {
            acc_str = "UNIMOD:" + std::to_string(mod->getUniModRecordId());
          }
          if (mod_map.find(name) == mod_map.end())
          {
            mod_map[name] = {acc_str, {}};
          }
          mod_map[name].second.push_back({static_cast<int32_t>(seq.size() + 1), ""}); // C-term
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
          for (const auto& pi : acc_positions.second)
          {
            (void)position_struct_b->Append();
            (void)pos_position_b->Append(pi.position);
            if (pi.amino_acid.empty())
            {
              (void)pos_amino_acid_b->AppendNull();
            }
            else
            {
              (void)pos_amino_acid_b->Append(pi.amino_acid);
            }
            (void)scores_list_b->Append(); // empty scores list per position
          }
        }
      }
      // else: empty modifications list (Append() already called)

      // === charge (int16, non-nullable) ===
      int charge = hit.getCharge();
      (void)charge_builder.Append(static_cast<int16_t>(charge));

      // === posterior_error_probability (float64, nullable) ===
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

      // === is_decoy (bool, non-nullable) ===
      if (hit.metaValueExists("target_decoy"))
      {
        std::string td = hit.getMetaValue("target_decoy").toString();
        (void)is_decoy_builder.Append(td.substr(0, 5) == "decoy");
      }
      else
      {
        (void)is_decoy_builder.Append(false); // default to false (target) when unknown
      }

      // === calculated_mz (float32, non-nullable) ===
      float calc_mz = 0.0f;
      if (charge > 0)
      {
        calc_mz = static_cast<float>(seq.getMZ(charge));
      }
      (void)calculated_mz_builder.Append(calc_mz);

      // === observed_mz (float32, non-nullable) ===
      float obs_mz = static_cast<float>(pep_id.getMZ());
      if (obs_mz != obs_mz) obs_mz = 0.0f; // NaN -> 0
      (void)observed_mz_builder.Append(obs_mz);

      // === mass_error_ppm (float32, nullable) ===
      if (charge > 0 && obs_mz > 0.0f && calc_mz > 0.0f)
      {
        float error_ppm = static_cast<float>((obs_mz - calc_mz) / calc_mz * 1e6);
        (void)mass_error_ppm_builder.Append(error_ppm);
      }
      else
      {
        (void)mass_error_ppm_builder.AppendNull();
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

      // === predicted_rt (float32, nullable) ===
      if (hit.metaValueExists("predicted_RT"))
      {
        (void)predicted_rt_builder.Append(static_cast<float>(static_cast<double>(hit.getMetaValue("predicted_RT"))));
      }
      else if (hit.metaValueExists("predicted_rt"))
      {
        (void)predicted_rt_builder.Append(static_cast<float>(static_cast<double>(hit.getMetaValue("predicted_rt"))));
      }
      else
      {
        (void)predicted_rt_builder.AppendNull();
      }

      // === run_file_name (non-nullable) ===
      // Prefer per-hit or per-PeptideIdentification file reference, fall back to identifier-level path
      {
        std::string run_file;
        if (hit.metaValueExists("reference_file_name"))
        {
          run_file = hit.getMetaValue("reference_file_name").toString();
        }
        else if (hit.metaValueExists("run_file_name"))
        {
          run_file = hit.getMetaValue("run_file_name").toString();
        }
        else if (pep_id.metaValueExists("reference_file_name"))
        {
          run_file = pep_id.getMetaValue("reference_file_name").toString();
        }
        if (run_file.empty() && fn_it != id_to_filename.end())
        {
          run_file = fn_it->second;
        }
        (void)run_file_name_builder.Append(run_file); // non-nullable, default to empty
      }

      // === cv_params (list<struct>, nullable - null for now) ===
      (void)cv_params_builder.AppendNull();

      // === scan (list<int32>, non-nullable) ===
      (void)scan_builder.Append();
      if (!spec_ref.empty())
      {
        Int scan_num = extractScan(spec_ref);
        if (scan_num >= 0)
        {
          (void)scan_value_b->Append(scan_num);
        }
      }
      // else: empty list (Append() already called for the list start)

      // === rt (float32, nullable) ===
      double rt = pep_id.getRT();
      if (rt == rt) // not NaN
      {
        (void)rt_builder.Append(static_cast<float>(rt));
      }
      else
      {
        (void)rt_builder.AppendNull();
      }

      // === ion_mobility (float32, nullable) ===
      if (hit.metaValueExists("ion_mobility"))
      {
        (void)ion_mobility_builder.Append(static_cast<float>(static_cast<double>(hit.getMetaValue("ion_mobility"))));
      }
      else if (hit.metaValueExists("IM"))
      {
        (void)ion_mobility_builder.Append(static_cast<float>(static_cast<double>(hit.getMetaValue("IM"))));
      }
      else if (pep_id.metaValueExists("ion_mobility"))
      {
        (void)ion_mobility_builder.Append(static_cast<float>(static_cast<double>(pep_id.getMetaValue("ion_mobility"))));
      }
      else if (pep_id.metaValueExists("IM"))
      {
        (void)ion_mobility_builder.Append(static_cast<float>(static_cast<double>(pep_id.getMetaValue("IM"))));
      }
      else
      {
        (void)ion_mobility_builder.AppendNull();
      }

      // === missed_cleavages (int16, nullable) ===
      if (hit.metaValueExists("missed_cleavages"))
      {
        (void)missed_cleavages_builder.Append(static_cast<int16_t>(static_cast<int>(hit.getMetaValue("missed_cleavages"))));
      }
      else
      {
        (void)missed_cleavages_builder.AppendNull();
      }

      // === protein_accessions (list<utf8>, nullable) ===
      (void)protein_accessions_builder.Append();
      for (const auto& ev : hit.getPeptideEvidences())
      {
        (void)pa_vb->Append(ev.getProteinAccession());
      }

      // === fragment annotation arrays (from PeakAnnotations) ===
      const auto& peak_annotations = hit.getPeakAnnotations();
      if (!peak_annotations.empty())
      {
        (void)mz_array_builder.Append();
        (void)intensity_array_builder.Append();
        (void)charge_array_builder.Append();
        (void)ion_type_array_builder.Append();
        for (const auto& pa : peak_annotations)
        {
          (void)mz_arr_vb->Append(static_cast<float>(pa.mz));
          (void)int_arr_vb->Append(static_cast<float>(pa.intensity));
          (void)chg_arr_vb->Append(pa.charge);
          (void)ion_arr_vb->Append(pa.annotation);
        }
      }
      else
      {
        (void)mz_array_builder.AppendNull();
        (void)intensity_array_builder.AppendNull();
        (void)charge_array_builder.AppendNull();
        (void)ion_type_array_builder.AppendNull();
      }

    } // end hit loop

  } // end peptide identification loop

  // Finalize all arrays
  std::shared_ptr<arrow::Array> arr_sequence, arr_peptidoform, arr_modifications;
  std::shared_ptr<arrow::Array> arr_charge, arr_pep, arr_is_decoy;
  std::shared_ptr<arrow::Array> arr_calc_mz, arr_obs_mz, arr_mass_error_ppm;
  std::shared_ptr<arrow::Array> arr_additional_scores;
  std::shared_ptr<arrow::Array> arr_predicted_rt, arr_run_file_name;
  std::shared_ptr<arrow::Array> arr_cv_params, arr_scan, arr_rt, arr_ion_mobility;
  std::shared_ptr<arrow::Array> arr_missed_cleavages;
  std::shared_ptr<arrow::Array> arr_protein_acc;

  status = sequence_builder.Finish(&arr_sequence);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: sequence_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = peptidoform_builder.Finish(&arr_peptidoform);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: peptidoform_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = modifications_builder.Finish(&arr_modifications);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: modifications_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = charge_builder.Finish(&arr_charge);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: charge_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = pep_builder.Finish(&arr_pep);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: pep_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = is_decoy_builder.Finish(&arr_is_decoy);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: is_decoy_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = calculated_mz_builder.Finish(&arr_calc_mz);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: calculated_mz_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = observed_mz_builder.Finish(&arr_obs_mz);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: observed_mz_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = mass_error_ppm_builder.Finish(&arr_mass_error_ppm);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: mass_error_ppm_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = additional_scores_builder.Finish(&arr_additional_scores);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: additional_scores_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = predicted_rt_builder.Finish(&arr_predicted_rt);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: predicted_rt_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = run_file_name_builder.Finish(&arr_run_file_name);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: run_file_name_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = cv_params_builder.Finish(&arr_cv_params);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: cv_params_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = scan_builder.Finish(&arr_scan);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: scan_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = rt_builder.Finish(&arr_rt);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: rt_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = ion_mobility_builder.Finish(&arr_ion_mobility);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: ion_mobility_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = missed_cleavages_builder.Finish(&arr_missed_cleavages);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: missed_cleavages_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = protein_accessions_builder.Finish(&arr_protein_acc);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: protein_accessions_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  std::shared_ptr<arrow::Array> arr_mz_array, arr_intensity_array, arr_charge_array, arr_ion_type_array;
  status = mz_array_builder.Finish(&arr_mz_array);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: mz_array_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = intensity_array_builder.Finish(&arr_intensity_array);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: intensity_array_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = charge_array_builder.Finish(&arr_charge_array);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: charge_array_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = ion_type_array_builder.Finish(&arr_ion_type_array);
  if (!status.ok()) { OPENMS_LOG_ERROR << "QPXFile: ion_type_array_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }

  // Build null arrays for columns not yet populated.
  // Use actual row count from a finalized array to avoid coupling with the pre-loop estimate.
  const int64_t actual_rows = arr_sequence->length();
  auto make_null_array = [&](const std::shared_ptr<arrow::DataType>& type) -> std::shared_ptr<arrow::Array>
  {
    auto result = arrow::MakeArrayOfNull(type, actual_rows);
    if (!result.ok())
    {
      OPENMS_LOG_ERROR << "QPXFile: MakeArrayOfNull failed: " << result.status().ToString() << std::endl;
      return nullptr;
    }
    return *result;
  };

  auto arr_cross_links = make_null_array(QPXPSMSchema::crossLinksType());
  auto arr_ion_mobility_array = make_null_array(arrow::list(arrow::float32()));

  if (!arr_cross_links || !arr_mz_array || !arr_intensity_array ||
      !arr_charge_array || !arr_ion_type_array || !arr_ion_mobility_array)
  {
    return nullptr;
  }

  // Build schema from registry
  auto schema = QPXPSMSchema::schema();

  auto table = arrow::Table::Make(schema, {
    arr_sequence, arr_peptidoform, arr_modifications,
    arr_charge, arr_pep, arr_is_decoy,
    arr_calc_mz, arr_obs_mz, arr_mass_error_ppm,
    arr_additional_scores, arr_predicted_rt, arr_run_file_name,
    arr_cv_params, arr_scan, arr_rt, arr_ion_mobility,
    arr_missed_cleavages, arr_protein_acc,
    arr_cross_links,
    arr_mz_array, arr_intensity_array, arr_charge_array,
    arr_ion_type_array, arr_ion_mobility_array
  });

  // Validate table against registry schema (strict -- write path must match exactly)
  auto validation = ArrowSchemaValidation::validate(table, QPXPSMSchema::schema());
  if (!validation.valid)
  {
    OPENMS_LOG_ERROR << "QPXFile: Schema validation failed: " << validation.toString() << "\n";
    return nullptr;
  }

  return table;
}


namespace
{
  /// Attach the canonical QPX "psm" file metadata (qpx_version, file_type="psm",
  /// UUID, creation date, scan_format, creator) to the schema of @p table.
  std::shared_ptr<arrow::Table> attachQPXPsmMetadata(const std::shared_ptr<arrow::Table>& table)
  {
    auto metadata = arrow::key_value_metadata({
      {"qpx_version", "1.0"},
      {"creator", "OpenMS"},
      {"file_type", "psm"},
      {"creation_date", DateTime::nowUTC().toString("yyyy-MM-ddThh:mm:ssZ")},
      {"uuid", std::string(ArrowIOHelpers::generateUuidV4())},
      {"scan_format", "scan"},
      {"software_provider", "OpenMS"}
    });
    return table->ReplaceSchemaMetadata(metadata);
  }
}

bool QPXFile::exportToParquet(
  const std::vector<ProteinIdentification>& protein_identifications,
  const PeptideIdentificationList& peptide_identifications,
  const String& filename,
  bool export_all_psms,
  const ParquetWriteConfig& config)
{
  auto table = exportPSMsToQPXArrow(protein_identifications, peptide_identifications, export_all_psms);
  if (!table)
  {
    OPENMS_LOG_ERROR << "QPXFile: Failed to create Arrow table" << std::endl;
    return false;
  }
  return exportToParquet(table, filename, config);
}

bool QPXFile::exportToParquet(
  const std::shared_ptr<arrow::Table>& table,
  const String& filename,
  const ParquetWriteConfig& config)
{
  if (!table)
  {
    OPENMS_LOG_ERROR << "QPXFile: null table passed to exportToParquet (" << filename << ")" << std::endl;
    return false;
  }

  // Guard: the table-taking overload attaches file_type="psm" metadata, so the
  // caller must actually pass a QPXPSMSchema table (not, e.g., the internal
  // PSMSchema produced by exportToArrow).
  auto validation = ArrowSchemaValidation::validate(table, QPXPSMSchema::schema(), ArrowSchemaValidation::Mode::Strict);
  if (!validation.valid)
  {
    OPENMS_LOG_ERROR << "QPXFile: table schema does not match QPXPSMSchema ("
                     << filename << "): " << validation.toString() << std::endl;
    return false;
  }

  return ArrowIOHelpers::writeTableToParquet(attachQPXPsmMetadata(table), filename, config);
}

} // namespace OpenMS
