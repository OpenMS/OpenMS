// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/PSMParquetIO.h>

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/FORMAT/ArrowIOHelpers.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/Scores.h>
#include <OpenMS/CHEMISTRY/ProForma.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/METADATA/SpectrumNativeIDParser.h>

#include <arrow/api.h>
#include <arrow/builder.h>

#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace OpenMS
{

namespace
{
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


std::shared_ptr<arrow::Table> PSMParquetIO::exportToArrow(
  const std::vector<ProteinIdentification>& protein_identifications,
  const PeptideIdentificationList& peptide_identifications,
  bool export_all_psms)
{
  // -- Simple column builders --
  arrow::StringBuilder sequence_builder, peptidoform_builder;
  arrow::StringBuilder reference_file_builder, score_type_builder;
  arrow::StringBuilder spectrum_reference_builder, cv_params_builder;
  arrow::StringBuilder run_identifier_builder;
  arrow::Int32Builder precursor_charge_builder, hit_index_builder, p_id_builder;
  arrow::BooleanBuilder is_decoy_builder, higher_score_better_builder;
  arrow::Int32Builder scan_builder;
  arrow::DoubleBuilder pep_builder, calculated_mz_builder, observed_mz_builder;
  arrow::DoubleBuilder rt_builder, ion_mobility_builder, predicted_rt_builder, score_builder;

  // -- protein_accessions list<struct{accession, aa_before, aa_after, start, end}> --
  auto pa_struct_type = std::static_pointer_cast<arrow::ListType>(PSMSchema::proteinAccessionsType())->value_type();
  auto pa_acc_b      = std::make_shared<arrow::StringBuilder>();
  auto pa_before_b   = std::make_shared<arrow::StringBuilder>();
  auto pa_after_b    = std::make_shared<arrow::StringBuilder>();
  auto pa_start_b    = std::make_shared<arrow::Int32Builder>();
  auto pa_end_b      = std::make_shared<arrow::Int32Builder>();
  auto pa_struct_b   = std::make_shared<arrow::StructBuilder>(
    pa_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{pa_acc_b, pa_before_b, pa_after_b, pa_start_b, pa_end_b});
  arrow::ListBuilder protein_accessions_builder(arrow::default_memory_pool(), pa_struct_b);

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

  // -- fragment annotation arrays --
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
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: sequence_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = peptidoform_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: peptidoform_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = precursor_charge_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: precursor_charge_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = pep_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: pep_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = is_decoy_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: is_decoy_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = higher_score_better_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: higher_score_better_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = calculated_mz_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: calculated_mz_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = observed_mz_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: observed_mz_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = predicted_rt_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: predicted_rt_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = reference_file_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: reference_file_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = cv_params_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: cv_params_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = scan_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: scan_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = rt_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: rt_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = ion_mobility_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: ion_mobility_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = spectrum_reference_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: spectrum_reference_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = score_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: score_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = score_type_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: score_type_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = hit_index_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: hit_index_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = p_id_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: p_id_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = run_identifier_builder.Reserve(num_rows);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: run_identifier_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }

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
    Constants::UserParam::FRAGMENT_ANNOTATION_USERPARAM
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

    std::string spec_ref = pep_id.getSpectrumReference();

    IDScoreSwitcherAlgorithm::ScoreSearchResult pep_result;
    try
    {
      pep_result = idsa.findScoreType(pep_id, IDScoreSwitcherAlgorithm::ScoreType::PEP);
    }
    catch (...) {}

    auto fn_it = id_to_filename.find(pep_id.getIdentifier());

    for (Size hit_idx = 0; hit_idx < num_hits; ++hit_idx)
    {
      const PeptideHit& hit = hits[hit_idx];
      const auto& seq = hit.getSequence();

      // === sequence (full modified OpenMS format for fast round-trip) ===
      (void)sequence_builder.Append(seq.toString());

      // === peptidoform (ProForma canonical for external interop) ===
      auto pf = ProForma::fromAASequence(seq);
      (void)peptidoform_builder.Append(ProForma::toString(pf, ProForma::WriteMode::CANONICAL));

      // === modifications (structured) ===
      (void)modifications_builder.Append();
      if (seq.isModified())
      {
        std::map<std::string, std::pair<std::string, std::vector<std::string>>> mod_map;

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

      // === is_decoy ===
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
        if (obs_mz == obs_mz)
        {
          (void)observed_mz_builder.Append(obs_mz);
        }
        else
        {
          (void)observed_mz_builder.AppendNull();
        }
      }

      // === additional_scores ===
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
        (void)pa_struct_b->Append();
        (void)pa_acc_b->Append(ev.getProteinAccession());

        const char b = ev.getAABefore();
        if (b == PeptideEvidence::UNKNOWN_AA) { (void)pa_before_b->AppendNull(); }
        else { (void)pa_before_b->Append(std::string(1, b)); }

        const char a = ev.getAAAfter();
        if (a == PeptideEvidence::UNKNOWN_AA) { (void)pa_after_b->AppendNull(); }
        else { (void)pa_after_b->Append(std::string(1, a)); }

        const int s = ev.getStart();
        if (s == PeptideEvidence::UNKNOWN_POSITION) { (void)pa_start_b->AppendNull(); }
        else { (void)pa_start_b->Append(s); }

        const int e = ev.getEnd();
        if (e == PeptideEvidence::UNKNOWN_POSITION) { (void)pa_end_b->AppendNull(); }
        else { (void)pa_end_b->Append(e); }
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
      if (rt == rt)
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

      // === hit_index ===
      (void)hit_index_builder.Append(static_cast<int32_t>(hit_idx));

      // === P_ID ===
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

      // === fragment annotation arrays ===
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
  std::shared_ptr<arrow::Array> arr_hit_index, arr_p_id;
  std::shared_ptr<arrow::Array> arr_psm_mvs, arr_spectrum_mvs;

  status = sequence_builder.Finish(&arr_sequence);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: sequence_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = peptidoform_builder.Finish(&arr_peptidoform);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: peptidoform_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = modifications_builder.Finish(&arr_modifications);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: modifications_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = precursor_charge_builder.Finish(&arr_charge);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: precursor_charge_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = pep_builder.Finish(&arr_pep);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: pep_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = is_decoy_builder.Finish(&arr_is_decoy);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: is_decoy_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = calculated_mz_builder.Finish(&arr_calc_mz);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: calculated_mz_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = observed_mz_builder.Finish(&arr_obs_mz);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: observed_mz_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = additional_scores_builder.Finish(&arr_additional_scores);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: additional_scores_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = protein_accessions_builder.Finish(&arr_protein_acc);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: protein_accessions_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = predicted_rt_builder.Finish(&arr_predicted_rt);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: predicted_rt_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = reference_file_builder.Finish(&arr_ref_file);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: reference_file_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = cv_params_builder.Finish(&arr_cv_params);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: cv_params_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = scan_builder.Finish(&arr_scan);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: scan_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = rt_builder.Finish(&arr_rt);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: rt_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = ion_mobility_builder.Finish(&arr_ion_mobility);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: ion_mobility_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = spectrum_reference_builder.Finish(&arr_spectrum_ref);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: spectrum_reference_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = score_builder.Finish(&arr_score);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: score_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = score_type_builder.Finish(&arr_score_type);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: score_type_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = higher_score_better_builder.Finish(&arr_higher_score_better);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: higher_score_better_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = hit_index_builder.Finish(&arr_hit_index);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: hit_index_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = p_id_builder.Finish(&arr_p_id);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: p_id_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = psm_metavalues_builder.Finish(&arr_psm_mvs);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: psm_metavalues_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = spectrum_metavalues_builder.Finish(&arr_spectrum_mvs);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: spectrum_metavalues_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }

  std::shared_ptr<arrow::Array> arr_run_identifier;
  status = run_identifier_builder.Finish(&arr_run_identifier);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: run_identifier_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  std::shared_ptr<arrow::Array> arr_psm_mz_array, arr_psm_intensity_array, arr_psm_charge_array, arr_psm_ion_type_array;
  status = psm_mz_array_builder.Finish(&arr_psm_mz_array);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: psm_mz_array_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = psm_intensity_array_builder.Finish(&arr_psm_intensity_array);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: psm_intensity_array_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = psm_charge_array_builder.Finish(&arr_psm_charge_array);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: psm_charge_array_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = psm_ion_type_array_builder.Finish(&arr_psm_ion_type_array);
  if (!status.ok()) { OPENMS_LOG_ERROR << "PSMParquetIO: psm_ion_type_array_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }

  auto schema = PSMSchema::schema();

  auto table = arrow::Table::Make(schema, {
    arr_sequence, arr_peptidoform, arr_modifications,
    arr_charge, arr_pep, arr_is_decoy,
    arr_calc_mz, arr_obs_mz, arr_additional_scores,
    arr_protein_acc, arr_predicted_rt, arr_ref_file,
    arr_cv_params, arr_scan, arr_rt, arr_ion_mobility,
    arr_spectrum_ref, arr_score, arr_score_type, arr_higher_score_better,
    arr_hit_index, arr_p_id, arr_psm_mvs, arr_spectrum_mvs,
    arr_run_identifier,
    arr_psm_mz_array, arr_psm_intensity_array, arr_psm_charge_array, arr_psm_ion_type_array
  });

  auto validation = ArrowSchemaValidation::validate(table, PSMSchema::schema());
  if (!validation.valid)
  {
    OPENMS_LOG_ERROR << "PSMParquetIO: Schema validation failed: " << validation.toString() << "\n";
    return nullptr;
  }

  return table;
}


bool PSMParquetIO::importFromArrow(
  const std::shared_ptr<arrow::Table>& table,
  std::vector<ProteinIdentification>& protein_identifications,
  PeptideIdentificationList& peptide_identifications)
{
  if (!table)
  {
    OPENMS_LOG_ERROR << "PSMParquetIO::importFromArrow: null table" << std::endl;
    return false;
  }

  auto combined_result = table->CombineChunks(arrow::default_memory_pool());
  if (!combined_result.ok())
  {
    OPENMS_LOG_ERROR << "PSMParquetIO::importFromArrow: Failed to combine chunks" << std::endl;
    return false;
  }
  const auto& tbl = *combined_result;
  int64_t num_rows = tbl->num_rows();

  auto psm_validation = ArrowSchemaValidation::validate(
    tbl, PSMSchema::schema(), ArrowSchemaValidation::Mode::Subset);
  if (!psm_validation.valid)
  {
    OPENMS_LOG_ERROR << "PSMParquetIO::importFromArrow: Incompatible PSM schema: "
                     << psm_validation.toString() << std::endl;
    return false;
  }
  if (num_rows == 0) { return true; }

  auto col_p_id = ArrowIOHelpers::getColumn(tbl, PSMSchema::PEPTIDE_IDENTIFICATION_INDEX);
  auto col_peptidoform = ArrowIOHelpers::getColumn(tbl, PSMSchema::PEPTIDOFORM, /*required=*/false);
  auto col_sequence = ArrowIOHelpers::getColumn(tbl, PSMSchema::SEQUENCE, /*required=*/false);
  auto col_charge = ArrowIOHelpers::getColumn(tbl, PSMSchema::PRECURSOR_CHARGE);
  auto col_score = ArrowIOHelpers::getColumn(tbl, PSMSchema::SCORE);
  auto col_score_type = ArrowIOHelpers::getColumn(tbl, PSMSchema::SCORE_TYPE);
  auto col_rt = ArrowIOHelpers::getColumn(tbl, PSMSchema::RT, /*required=*/false);
  auto col_mz = ArrowIOHelpers::getColumn(tbl, PSMSchema::OBSERVED_MZ, /*required=*/false);
  auto col_spec_ref = ArrowIOHelpers::getColumn(tbl, PSMSchema::SPECTRUM_REFERENCE, /*required=*/false);
  auto col_run_id = ArrowIOHelpers::getColumn(tbl, PSMSchema::RUN_IDENTIFIER, /*required=*/false);
  auto col_is_decoy = ArrowIOHelpers::getColumn(tbl, PSMSchema::IS_DECOY, /*required=*/false);
  auto col_protein_accs = ArrowIOHelpers::getColumn(tbl, PSMSchema::PROTEIN_ACCESSIONS, /*required=*/false);
  auto col_additional_scores = ArrowIOHelpers::getColumn(tbl, PSMSchema::ADDITIONAL_SCORES, /*required=*/false);
  auto col_psm_metavalues = ArrowIOHelpers::getColumn(tbl, PSMSchema::PSM_METAVALUES, /*required=*/false);
  auto col_spectrum_metavalues = ArrowIOHelpers::getColumn(tbl, PSMSchema::SPECTRUM_METAVALUES, /*required=*/false);
  auto col_predicted_rt = ArrowIOHelpers::getColumn(tbl, PSMSchema::PREDICTED_RT, /*required=*/false);
  auto col_ion_mobility = ArrowIOHelpers::getColumn(tbl, PSMSchema::ION_MOBILITY, /*required=*/false);
  auto col_hsb = ArrowIOHelpers::getColumn(tbl, PSMSchema::HIGHER_SCORE_BETTER, /*required=*/false);
  auto col_scan = ArrowIOHelpers::getColumn(tbl, PSMSchema::SCAN, /*required=*/false);
  auto col_ref_file = ArrowIOHelpers::getColumn(tbl, PSMSchema::REFERENCE_FILE_NAME, /*required=*/false);

  if (!col_p_id || !col_charge || !col_score || !col_score_type)
  {
    OPENMS_LOG_ERROR << "PSMParquetIO::importFromArrow: Missing required PSM columns" << std::endl;
    return false;
  }

  std::unordered_map<std::string, bool> higher_score_better_lookup;
  for (const auto& prot_id : protein_identifications)
  {
    higher_score_better_lookup[prot_id.getIdentifier()] = prot_id.isHigherScoreBetter();
  }

  std::unordered_map<int32_t, size_t> p_id_to_idx;
  const size_t pep_ids_start_size = peptide_identifications.size();

  std::unordered_map<std::string, std::string> run_id_to_ref_file;

  for (int64_t row = 0; row < num_rows; ++row)
  {
    int32_t p_id = ArrowIOHelpers::getInt32Value(col_p_id, row, -1);

    auto [it, inserted] = p_id_to_idx.try_emplace(p_id, peptide_identifications.size());
    if (inserted)
    {
      peptide_identifications.emplace_back();
      PeptideIdentification& pid = peptide_identifications.back();
      pid.setScoreType(ArrowIOHelpers::getStringValue(col_score_type, row));

      if (col_run_id && !ArrowIOHelpers::isNull(col_run_id, row))
      {
        pid.setIdentifier(ArrowIOHelpers::getStringValue(col_run_id, row));
      }

      if (col_hsb && !ArrowIOHelpers::isNull(col_hsb, row))
      {
        pid.setHigherScoreBetter(ArrowIOHelpers::getBoolValue(col_hsb, row, true));
      }
      else if (col_run_id && !ArrowIOHelpers::isNull(col_run_id, row))
      {
        auto hsb_it = higher_score_better_lookup.find(pid.getIdentifier());
        pid.setHigherScoreBetter(
          hsb_it != higher_score_better_lookup.end() ? hsb_it->second : true);
      }
      else
      {
        pid.setHigherScoreBetter(true);
      }

      if (col_rt && !ArrowIOHelpers::isNull(col_rt, row)) { pid.setRT(ArrowIOHelpers::getDoubleValue(col_rt, row)); }
      if (col_mz && !ArrowIOHelpers::isNull(col_mz, row)) { pid.setMZ(ArrowIOHelpers::getDoubleValue(col_mz, row)); }
      if (col_spec_ref && !ArrowIOHelpers::isNull(col_spec_ref, row))
      {
        pid.setSpectrumReference(ArrowIOHelpers::getStringValue(col_spec_ref, row));
      }
      if (col_spectrum_metavalues)
      {
        ArrowIOHelpers::readMetaValues(col_spectrum_metavalues, row, pid);
      }
    }

    PeptideHit hit;

    // Prefer the sequence column (native OpenMS modified format) for fast reconstruction.
    // Falls back to peptidoform (ProForma) only if sequence is unavailable.
    bool sequence_set = false;
    if (col_sequence && !ArrowIOHelpers::isNull(col_sequence, row))
    {
      const String seq_str = ArrowIOHelpers::getStringValue(col_sequence, row);
      if (!seq_str.empty())
      {
        hit.setSequence(AASequence::fromString(seq_str));
        sequence_set = true;
      }
    }
    if (!sequence_set && col_peptidoform && !ArrowIOHelpers::isNull(col_peptidoform, row))
    {
      const String peptidoform_str = ArrowIOHelpers::getStringValue(col_peptidoform, row);
      if (!peptidoform_str.empty())
      {
        try
        {
          auto pf = ProForma::parse(peptidoform_str);
          hit.setSequence(ProForma::toAASequence(pf, ProForma::ConversionPolicy::BEST_EFFORT));
          sequence_set = true;
        }
        catch (...)
        {
        }
      }
    }

    hit.setCharge(static_cast<Int>(ArrowIOHelpers::getInt32Value(col_charge, row, 0)));
    hit.setScore(ArrowIOHelpers::getDoubleValue(col_score, row, 0.0));

    if (col_is_decoy && !ArrowIOHelpers::isNull(col_is_decoy, row))
    {
      bool is_decoy = ArrowIOHelpers::getBoolValue(col_is_decoy, row, false);
      hit.setMetaValue("target_decoy", is_decoy ? "decoy" : "target");
    }

    if (col_protein_accs && !ArrowIOHelpers::isNull(col_protein_accs, row))
    {
      auto list_arr = std::static_pointer_cast<arrow::ListArray>(col_protein_accs);
      auto struct_arr = std::static_pointer_cast<arrow::StructArray>(list_arr->values());
      auto acc_arr    = std::static_pointer_cast<arrow::StringArray>(struct_arr->GetFieldByName("accession"));
      auto before_arr = std::static_pointer_cast<arrow::StringArray>(struct_arr->GetFieldByName("aa_before"));
      auto after_arr  = std::static_pointer_cast<arrow::StringArray>(struct_arr->GetFieldByName("aa_after"));
      auto start_arr  = std::static_pointer_cast<arrow::Int32Array>(struct_arr->GetFieldByName("start"));
      auto end_arr    = std::static_pointer_cast<arrow::Int32Array>(struct_arr->GetFieldByName("end"));
      int64_t lstart = list_arr->value_offset(row);
      int64_t lend = lstart + list_arr->value_length(row);
      for (int64_t k = lstart; k < lend; ++k)
      {
        PeptideEvidence ev;
        ev.setProteinAccession(acc_arr->GetString(k));
        const std::string before_s = before_arr->IsNull(k) ? std::string{} : before_arr->GetString(k);
        ev.setAABefore(before_s.empty() ? PeptideEvidence::UNKNOWN_AA : before_s[0]);
        const std::string after_s = after_arr->IsNull(k) ? std::string{} : after_arr->GetString(k);
        ev.setAAAfter(after_s.empty() ? PeptideEvidence::UNKNOWN_AA : after_s[0]);
        ev.setStart(start_arr->IsNull(k) ? PeptideEvidence::UNKNOWN_POSITION : start_arr->Value(k));
        ev.setEnd  (end_arr  ->IsNull(k) ? PeptideEvidence::UNKNOWN_POSITION : end_arr  ->Value(k));
        hit.addPeptideEvidence(ev);
      }
    }

    if (col_additional_scores && !ArrowIOHelpers::isNull(col_additional_scores, row))
    {
      auto list_arr = std::static_pointer_cast<arrow::ListArray>(col_additional_scores);
      auto struct_arr = std::static_pointer_cast<arrow::StructArray>(list_arr->values());
      auto names_arr = std::static_pointer_cast<arrow::StringArray>(struct_arr->GetFieldByName("score_name"));
      auto values_arr = std::static_pointer_cast<arrow::DoubleArray>(struct_arr->GetFieldByName("score_value"));
      int64_t start = list_arr->value_offset(row);
      int64_t end = start + list_arr->value_length(row);
      for (int64_t k = start; k < end; ++k)
      {
        hit.setMetaValue(names_arr->GetString(k), values_arr->Value(k));
      }
    }

    if (col_predicted_rt && !ArrowIOHelpers::isNull(col_predicted_rt, row))
    {
      hit.setMetaValue("predicted_RT", ArrowIOHelpers::getDoubleValue(col_predicted_rt, row));
    }
    if (col_ion_mobility && !ArrowIOHelpers::isNull(col_ion_mobility, row))
    {
      hit.setMetaValue("ion_mobility", ArrowIOHelpers::getDoubleValue(col_ion_mobility, row));
    }
    if (col_scan && !ArrowIOHelpers::isNull(col_scan, row))
    {
      hit.setMetaValue("scan", static_cast<int>(ArrowIOHelpers::getInt32Value(col_scan, row)));
    }
    if (col_ref_file && !ArrowIOHelpers::isNull(col_ref_file, row))
    {
      const String ref_file = ArrowIOHelpers::getStringValue(col_ref_file, row);
      hit.setMetaValue("reference_file_name", ref_file);
      if (!ref_file.empty())
      {
        const std::string& run_id = peptide_identifications[it->second].getIdentifier();
        run_id_to_ref_file.try_emplace(run_id, ref_file);
      }
    }

    if (col_psm_metavalues)
    {
      static const std::unordered_set<std::string> psm_excluded_mvs =
        {"target_decoy", "predicted_RT", "predicted_rt", "ion_mobility", "IM",
         "scan", "reference_file_name"};
      ArrowIOHelpers::readMetaValues(col_psm_metavalues, row, hit, psm_excluded_mvs);
    }

    peptide_identifications[it->second].getHits().push_back(std::move(hit));
  }

  // Append a shell ProteinIdentification for any run_identifier without one.
  std::unordered_set<std::string> known;
  for (const auto& p : protein_identifications) { known.insert(p.getIdentifier()); }
  for (size_t i = pep_ids_start_size; i < peptide_identifications.size(); ++i)
  {
    const auto& pid = peptide_identifications[i];
    const std::string& id = pid.getIdentifier();
    if (!id.empty() && known.insert(id).second)
    {
      ProteinIdentification shell;
      shell.setIdentifier(id);
      shell.setScoreType(pid.getScoreType());
      shell.setHigherScoreBetter(pid.isHigherScoreBetter());
      auto rf_it = run_id_to_ref_file.find(id);
      if (rf_it != run_id_to_ref_file.end())
      {
        shell.setPrimaryMSRunPath({rf_it->second});
      }
      protein_identifications.push_back(std::move(shell));
    }
  }

  return true;
}

} // namespace OpenMS
