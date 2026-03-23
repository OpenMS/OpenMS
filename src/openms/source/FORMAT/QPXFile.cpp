// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/QPXFile.h>

#ifdef WITH_PARQUET

#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
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

  // -- mz_array, intensity_array, charge_array, ion_type_array, ion_mobility_array (null for now) --
  // We build null arrays later

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
    "scan", "reference_file_name"
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
      if (fn_it != id_to_filename.end())
      {
        (void)run_file_name_builder.Append(fn_it->second);
      }
      else
      {
        (void)run_file_name_builder.Append(""); // non-nullable, default to empty
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

  // Build null arrays for columns not yet populated
  auto make_null_array = [&](const std::shared_ptr<arrow::DataType>& type) -> std::shared_ptr<arrow::Array>
  {
    auto result = arrow::MakeArrayOfNull(type, num_rows);
    if (!result.ok())
    {
      OPENMS_LOG_ERROR << "QPXFile: MakeArrayOfNull failed: " << result.status().ToString() << std::endl;
      return nullptr;
    }
    return *result;
  };

  auto arr_cross_links = make_null_array(QPXPSMSchema::crossLinksType());
  auto arr_mz_array = make_null_array(arrow::list(arrow::float32()));
  auto arr_intensity_array = make_null_array(arrow::list(arrow::float32()));
  auto arr_charge_array = make_null_array(arrow::list(arrow::int32()));
  auto arr_ion_type_array = make_null_array(arrow::list(arrow::utf8()));
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
    OPENMS_LOG_ERROR << "QPXFile: Failed to write Parquet: "
                     << write_status.ToString() << std::endl;
    return false;
  }

  return true;
}

} // namespace OpenMS

#endif // WITH_PARQUET
