// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ConsensusMapArrowExport.h>

#ifdef WITH_PARQUET

#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/Scores.h>
#include <OpenMS/CHEMISTRY/ProForma.h>
#include <OpenMS/METADATA/SpectrumLookup.h>
#include <OpenMS/METADATA/SpectrumNativeIDParser.h>

#include <arrow/api.h>
#include <arrow/builder.h>
#include <arrow/io/file.h>
#include <parquet/arrow/writer.h>
#include <parquet/properties.h>

#include <limits>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>

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


std::shared_ptr<arrow::Table> ConsensusMapArrowExport::exportToArrow(const ConsensusMap& cmap)
{
  // -- Simple column builders --
  arrow::StringBuilder sequence_builder, peptidoform_builder, anchor_protein_builder;
  arrow::StringBuilder reference_file_builder, scan_builder, score_type_builder;
  arrow::StringBuilder spectrum_reference_builder, scan_ref_file_builder;
  arrow::StringBuilder cv_params_builder, additional_intensities_builder;
  arrow::Int32Builder precursor_charge_builder, is_decoy_builder, unique_builder;
  arrow::FloatBuilder calculated_mz_builder, observed_mz_builder, rt_builder;
  arrow::FloatBuilder quality_builder, rt_start_builder, rt_stop_builder;
  arrow::FloatBuilder ion_mobility_builder, predicted_rt_builder;
  arrow::FloatBuilder start_ion_mobility_builder, stop_ion_mobility_builder;
  arrow::DoubleBuilder pep_builder, score_builder, pg_qvalue_builder;

  // -- List<utf8> builders for pg_accessions, gg_accessions, gg_names --
  auto pg_acc_vb = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder pg_accessions_builder(arrow::default_memory_pool(), pg_acc_vb);

  auto gg_acc_vb = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder gg_accessions_builder(arrow::default_memory_pool(), gg_acc_vb);

  auto gg_names_vb = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder gg_names_builder(arrow::default_memory_pool(), gg_names_vb);

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

  // -- Intensity list<struct{sample_accession, channel, intensity}> --
  auto intensity_struct_type = arrow::struct_({
    arrow::field("sample_accession", arrow::utf8()),
    arrow::field("channel", arrow::utf8()),
    arrow::field("intensity", arrow::float32())
  });
  auto intensity_sa_b = std::make_shared<arrow::StringBuilder>();
  auto intensity_ch_b = std::make_shared<arrow::StringBuilder>();
  auto intensity_val_b = std::make_shared<arrow::FloatBuilder>();
  auto intensity_struct_b = std::make_shared<arrow::StructBuilder>(
    intensity_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{intensity_sa_b, intensity_ch_b, intensity_val_b});
  arrow::ListBuilder intensities_builder(arrow::default_memory_pool(), intensity_struct_b);

  // -- feature_metavalues list<struct{name, value, value_type}> --
  auto fmv_name_b = std::make_shared<arrow::StringBuilder>();
  auto fmv_value_b = std::make_shared<arrow::StringBuilder>();
  auto fmv_type_b = std::make_shared<arrow::StringBuilder>();
  auto fmv_struct_type = arrow::struct_({
    arrow::field("name", arrow::utf8()),
    arrow::field("value", arrow::utf8()),
    arrow::field("value_type", arrow::utf8())
  });
  auto fmv_struct_b = std::make_shared<arrow::StructBuilder>(
    fmv_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{fmv_name_b, fmv_value_b, fmv_type_b});
  arrow::ListBuilder feature_metavalues_builder(arrow::default_memory_pool(), fmv_struct_b);

  // Reserve capacity for all builders to avoid repeated allocations
  // This is where most memory allocation occurs - failures here indicate OOM
  arrow::Status status;
  const Size num_features = cmap.size();

  // Reserve for simple column builders (one value per feature)
  status = sequence_builder.Reserve(num_features);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: sequence_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = peptidoform_builder.Reserve(num_features);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: peptidoform_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = anchor_protein_builder.Reserve(num_features);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: anchor_protein_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = reference_file_builder.Reserve(num_features);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: reference_file_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = scan_builder.Reserve(num_features);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: scan_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = score_type_builder.Reserve(num_features);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: score_type_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = spectrum_reference_builder.Reserve(num_features);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: spectrum_reference_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = scan_ref_file_builder.Reserve(num_features);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: scan_ref_file_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = cv_params_builder.Reserve(num_features);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: cv_params_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = additional_intensities_builder.Reserve(num_features);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: additional_intensities_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }

  status = precursor_charge_builder.Reserve(num_features);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: precursor_charge_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = is_decoy_builder.Reserve(num_features);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: is_decoy_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = unique_builder.Reserve(num_features);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: unique_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }

  status = calculated_mz_builder.Reserve(num_features);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: calculated_mz_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = observed_mz_builder.Reserve(num_features);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: observed_mz_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = rt_builder.Reserve(num_features);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: rt_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = quality_builder.Reserve(num_features);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: quality_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = rt_start_builder.Reserve(num_features);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: rt_start_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = rt_stop_builder.Reserve(num_features);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: rt_stop_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = ion_mobility_builder.Reserve(num_features);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: ion_mobility_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = predicted_rt_builder.Reserve(num_features);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: predicted_rt_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = start_ion_mobility_builder.Reserve(num_features);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: start_ion_mobility_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = stop_ion_mobility_builder.Reserve(num_features);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: stop_ion_mobility_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }

  status = pep_builder.Reserve(num_features);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: pep_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = score_builder.Reserve(num_features);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: score_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }
  status = pg_qvalue_builder.Reserve(num_features);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: pg_qvalue_builder Reserve failed: " << status.ToString() << std::endl; return nullptr; }

  // Build column header lookup
  const auto& column_headers = cmap.getColumnHeaders();

  // Build protein group q-value lookup and gene info lookup
  std::unordered_map<std::string, double> pg_qvalue_lookup;
  std::unordered_map<std::string, std::vector<std::string>> gg_accessions_lookup;
  std::unordered_map<std::string, std::vector<std::string>> gg_names_lookup;

  for (const auto& prot_id : cmap.getProteinIdentifications())
  {
    // Gene info from ProteinHit metavalues
    for (const auto& ph : prot_id.getHits())
    {
      const std::string& acc = ph.getAccession();
      if (ph.metaValueExists("gene_accession"))
      {
        gg_accessions_lookup[acc].push_back(ph.getMetaValue("gene_accession").toString());
      }
      if (ph.metaValueExists("gene_name"))
      {
        gg_names_lookup[acc].push_back(ph.getMetaValue("gene_name").toString());
      }
    }

    for (const auto& pg : prot_id.getProteinGroups())
    {
      for (const auto& acc : pg.accessions)
      {
        pg_qvalue_lookup[acc] = pg.probability;
      }
    }
    for (const auto& ig : prot_id.getIndistinguishableProteins())
    {
      for (const auto& acc : ig.accessions)
      {
        if (pg_qvalue_lookup.find(acc) == pg_qvalue_lookup.end())
        {
          pg_qvalue_lookup[acc] = ig.probability;
        }
      }
    }
  }

  // Metavalue keys excluded from feature_metavalues (they have dedicated columns)
  static const std::unordered_set<std::string> excluded_feature_mvs = {
    "rt_start", "rt_stop", "start_ion_mobility", "stop_ion_mobility"
  };
  // Metavalue keys excluded from additional_scores on PeptideHit
  static const std::unordered_set<std::string> excluded_hit_mvs = {
    "target_decoy", "predicted_RT", "predicted_rt"
  };

  IDScoreSwitcherAlgorithm idsa;

  for (const auto& cf : cmap)
  {
    const auto& pep_ids = cf.getPeptideIdentifications();
    bool has_id = !pep_ids.empty() && !pep_ids[0].getHits().empty();
    const PeptideHit* best_hit = has_id ? &pep_ids[0].getHits()[0] : nullptr;

    // === sequence / peptidoform ===
    if (best_hit)
    {
      const auto& seq = best_hit->getSequence();
      (void)sequence_builder.Append(seq.toUnmodifiedString());
      auto pf = ProForma::fromAASequence(seq);
      (void)peptidoform_builder.Append(ProForma::toString(pf, ProForma::WriteMode::CANONICAL));

      // calculated_mz
      int charge = cf.getCharge();
      if (charge > 0)
      {
        (void)calculated_mz_builder.Append(static_cast<float>(seq.getMZ(charge)));
      }
      else
      {
        (void)calculated_mz_builder.AppendNull();
      }

      // === modifications ===
      (void)modifications_builder.Append();
      if (seq.isModified())
      {
        // Collect modifications keyed by name to group positions
        std::map<std::string, std::pair<std::string, std::vector<std::string>>> mod_map; // name -> (accession, [position_strings])

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
    }
    else
    {
      (void)sequence_builder.AppendNull();
      (void)peptidoform_builder.AppendNull();
      (void)calculated_mz_builder.AppendNull();
      // empty modifications list
      (void)modifications_builder.Append();
    }

    (void)precursor_charge_builder.Append(cf.getCharge());
    (void)observed_mz_builder.Append(static_cast<float>(cf.getMZ()));
    (void)rt_builder.Append(static_cast<float>(cf.getRT()));
    (void)quality_builder.Append(static_cast<float>(cf.getQuality()));

    // === PEP ===
    if (has_id)
    {
      auto result = idsa.findScoreType(pep_ids[0], IDScoreSwitcherAlgorithm::ScoreType::PEP);
      if (!result.score_name.empty() && best_hit)
      {
        if (result.is_main_score_type)
        {
          (void)pep_builder.Append(best_hit->getScore());
        }
        else if (best_hit->metaValueExists(result.score_name))
        {
          (void)pep_builder.Append(static_cast<double>(best_hit->getMetaValue(result.score_name)));
        }
        else
        {
          (void)pep_builder.AppendNull();
        }
      }
      else
      {
        (void)pep_builder.AppendNull();
      }
    }
    else
    {
      (void)pep_builder.AppendNull();
    }

    // === is_decoy ===
    if (best_hit && best_hit->metaValueExists("target_decoy"))
    {
      std::string td = best_hit->getMetaValue("target_decoy").toString();
      (void)is_decoy_builder.Append(td.substr(0, 5) == "decoy" ? 1 : 0);
    }
    else
    {
      (void)is_decoy_builder.AppendNull();
    }

    // === additional_scores from best hit ===
    (void)additional_scores_builder.Append();
    if (best_hit)
    {
      std::vector<String> keys;
      best_hit->getKeys(keys);
      for (const auto& key : keys)
      {
        if (excluded_hit_mvs.count(key)) continue;
        const DataValue& val = best_hit->getMetaValue(key);
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

    // === predicted_rt ===
    if (best_hit && best_hit->metaValueExists("predicted_RT"))
    {
      (void)predicted_rt_builder.Append(static_cast<float>(static_cast<double>(best_hit->getMetaValue("predicted_RT"))));
    }
    else if (best_hit && best_hit->metaValueExists("predicted_rt"))
    {
      (void)predicted_rt_builder.Append(static_cast<float>(static_cast<double>(best_hit->getMetaValue("predicted_rt"))));
    }
    else
    {
      (void)predicted_rt_builder.AppendNull();
    }

    // === scan / spectrum_reference ===
    std::string spec_ref;
    if (!pep_ids.empty() && pep_ids[0].metaValueExists("spectrum_reference"))
    {
      spec_ref = pep_ids[0].getMetaValue("spectrum_reference").toString();
    }
    // Use null for missing values (not empty string) to align with Python and semantic correctness
    if (spec_ref.empty())
    {
      (void)spectrum_reference_builder.AppendNull();
    }
    else
    {
      (void)spectrum_reference_builder.Append(spec_ref);
    }

    // Extract scan number from native ID
    if (!spec_ref.empty())
    {
      Int scan_num = extractScan(spec_ref);
      if (scan_num >= 0)
      {
        (void)scan_builder.Append(std::to_string(scan_num));
      }
      else
      {
        (void)scan_builder.AppendNull();
      }
    }
    else
    {
      (void)scan_builder.AppendNull();
    }

    // === ion_mobility, start/stop ===
    if (!pep_ids.empty() && pep_ids[0].metaValueExists("ion_mobility"))
    {
      (void)ion_mobility_builder.Append(static_cast<float>(static_cast<double>(pep_ids[0].getMetaValue("ion_mobility"))));
    }
    else if (!pep_ids.empty() && pep_ids[0].metaValueExists("IM"))
    {
      (void)ion_mobility_builder.Append(static_cast<float>(static_cast<double>(pep_ids[0].getMetaValue("IM"))));
    }
    else
    {
      (void)ion_mobility_builder.AppendNull();
    }

    if (cf.metaValueExists("start_ion_mobility"))
    {
      (void)start_ion_mobility_builder.Append(static_cast<float>(static_cast<double>(cf.getMetaValue("start_ion_mobility"))));
    }
    else
    {
      (void)start_ion_mobility_builder.AppendNull();
    }

    if (cf.metaValueExists("stop_ion_mobility"))
    {
      (void)stop_ion_mobility_builder.Append(static_cast<float>(static_cast<double>(cf.getMetaValue("stop_ion_mobility"))));
    }
    else
    {
      (void)stop_ion_mobility_builder.AppendNull();
    }

    // === reference_file_name, cv_params, additional_intensities (reserved) ===
    // Use null for columns that don't have meaningful values in C++ export
    // (reference_file_name is a parameter in Python API but not in C++)
    (void)reference_file_builder.AppendNull();
    (void)scan_ref_file_builder.AppendNull();
    (void)cv_params_builder.AppendNull();
    (void)additional_intensities_builder.AppendNull();

    // === Score / score_type ===
    if (best_hit)
    {
      (void)score_builder.Append(best_hit->getScore());
      (void)score_type_builder.Append(pep_ids[0].getScoreType());
    }
    else
    {
      (void)score_builder.AppendNull();
      (void)score_type_builder.AppendNull();
    }

    // === Protein accessions ===
    // Collect unique protein accessions, preserving order (first occurrence for anchor)
    std::vector<std::string> protein_accs;
    std::unordered_set<std::string> seen_accs;
    if (best_hit)
    {
      for (const auto& ev : best_hit->getPeptideEvidences())
      {
        const std::string& acc = ev.getProteinAccession();
        if (seen_accs.insert(acc).second) // only add if not already seen
        {
          protein_accs.push_back(acc);
        }
      }
    }

    (void)pg_accessions_builder.Append();
    for (const auto& acc : protein_accs)
    {
      (void)pg_acc_vb->Append(acc);
    }

    // anchor_protein: null if no proteins
    if (protein_accs.empty())
    {
      (void)anchor_protein_builder.AppendNull();
    }
    else
    {
      (void)anchor_protein_builder.Append(protein_accs[0]);
    }

    // unique: 1 if exactly one unique protein accession, 0 otherwise
    (void)unique_builder.Append(protein_accs.size() == 1 ? 1 : 0);

    // === pg_global_qvalue ===
    // Use minimum q-value across all protein accessions (not first match which is order-dependent)
    double min_qval = std::numeric_limits<double>::max();
    bool found_qval = false;
    for (const auto& acc : protein_accs)
    {
      auto it = pg_qvalue_lookup.find(acc);
      if (it != pg_qvalue_lookup.end())
      {
        if (it->second < min_qval)
        {
          min_qval = it->second;
        }
        found_qval = true;
      }
    }
    if (found_qval)
    {
      (void)pg_qvalue_builder.Append(min_qval);
    }
    else
    {
      (void)pg_qvalue_builder.AppendNull();
    }

    // === gg_accessions / gg_names ===
    std::unordered_set<std::string> seen_gg_acc, seen_gg_nam;
    (void)gg_accessions_builder.Append();
    (void)gg_names_builder.Append();
    for (const auto& acc : protein_accs)
    {
      auto it_a = gg_accessions_lookup.find(acc);
      if (it_a != gg_accessions_lookup.end())
      {
        for (const auto& ga : it_a->second)
        {
          if (seen_gg_acc.insert(ga).second)
          {
            (void)gg_acc_vb->Append(ga);
          }
        }
      }
      auto it_n = gg_names_lookup.find(acc);
      if (it_n != gg_names_lookup.end())
      {
        for (const auto& gn : it_n->second)
        {
          if (seen_gg_nam.insert(gn).second)
          {
            (void)gg_names_vb->Append(gn);
          }
        }
      }
    }

    // === rt_start / rt_stop ===
    if (cf.metaValueExists("rt_start"))
    {
      (void)rt_start_builder.Append(static_cast<float>(static_cast<double>(cf.getMetaValue("rt_start"))));
    }
    else if (cf.getWidth() > 0)
    {
      (void)rt_start_builder.Append(static_cast<float>(cf.getRT() - cf.getWidth() / 2.0));
    }
    else
    {
      (void)rt_start_builder.AppendNull();
    }

    if (cf.metaValueExists("rt_stop"))
    {
      (void)rt_stop_builder.Append(static_cast<float>(static_cast<double>(cf.getMetaValue("rt_stop"))));
    }
    else if (cf.getWidth() > 0)
    {
      (void)rt_stop_builder.Append(static_cast<float>(cf.getRT() + cf.getWidth() / 2.0));
    }
    else
    {
      (void)rt_stop_builder.AppendNull();
    }

    // === intensities ===
    (void)intensities_builder.Append();
    for (const auto& fh : cf.getFeatures())
    {
      auto it = column_headers.find(fh.getMapIndex());
      if (it != column_headers.end())
      {
        (void)intensity_struct_b->Append();
        (void)intensity_sa_b->Append(it->second.filename);
        (void)intensity_ch_b->Append(it->second.label);
        (void)intensity_val_b->Append(static_cast<float>(fh.getIntensity()));
      }
    }

    // === feature_metavalues ===
    (void)feature_metavalues_builder.Append();
    {
      std::vector<String> cf_keys;
      cf.getKeys(cf_keys);
      for (const auto& key : cf_keys)
      {
        if (excluded_feature_mvs.count(key)) continue;
        const DataValue& val = cf.getMetaValue(key);
        (void)fmv_struct_b->Append();
        (void)fmv_name_b->Append(key);
        (void)fmv_value_b->Append(val.toString());

        switch (val.valueType())
        {
          case DataValue::INT_VALUE: (void)fmv_type_b->Append("int"); break;
          case DataValue::DOUBLE_VALUE: (void)fmv_type_b->Append("float"); break;
          case DataValue::STRING_VALUE: (void)fmv_type_b->Append("str"); break;
          default: (void)fmv_type_b->Append("str"); break;
        }
      }
    }
  } // end feature loop

  // Finalize all arrays
  std::shared_ptr<arrow::Array> arr_sequence, arr_peptidoform, arr_modifications;
  std::shared_ptr<arrow::Array> arr_charge, arr_calc_mz, arr_obs_mz, arr_rt;
  std::shared_ptr<arrow::Array> arr_pep, arr_is_decoy, arr_additional_scores;
  std::shared_ptr<arrow::Array> arr_predicted_rt, arr_ref_file, arr_cv_params;
  std::shared_ptr<arrow::Array> arr_scan, arr_ion_mobility;
  std::shared_ptr<arrow::Array> arr_start_im, arr_stop_im;
  std::shared_ptr<arrow::Array> arr_intensities, arr_additional_intensities;
  std::shared_ptr<arrow::Array> arr_pg_acc, arr_anchor, arr_unique, arr_pg_qval;
  std::shared_ptr<arrow::Array> arr_gg_acc, arr_gg_names;
  std::shared_ptr<arrow::Array> arr_scan_ref_file, arr_rt_start, arr_rt_stop;
  std::shared_ptr<arrow::Array> arr_quality, arr_score, arr_score_type;
  std::shared_ptr<arrow::Array> arr_spectrum_ref, arr_feature_mvs;

  status = sequence_builder.Finish(&arr_sequence);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: sequence_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = peptidoform_builder.Finish(&arr_peptidoform);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: peptidoform_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = modifications_builder.Finish(&arr_modifications);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: modifications_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = precursor_charge_builder.Finish(&arr_charge);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: precursor_charge_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = calculated_mz_builder.Finish(&arr_calc_mz);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: calculated_mz_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = observed_mz_builder.Finish(&arr_obs_mz);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: observed_mz_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = rt_builder.Finish(&arr_rt);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: rt_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = pep_builder.Finish(&arr_pep);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: pep_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = is_decoy_builder.Finish(&arr_is_decoy);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: is_decoy_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = additional_scores_builder.Finish(&arr_additional_scores);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: additional_scores_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = predicted_rt_builder.Finish(&arr_predicted_rt);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: predicted_rt_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = reference_file_builder.Finish(&arr_ref_file);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: reference_file_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = cv_params_builder.Finish(&arr_cv_params);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: cv_params_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = scan_builder.Finish(&arr_scan);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: scan_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = ion_mobility_builder.Finish(&arr_ion_mobility);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: ion_mobility_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = start_ion_mobility_builder.Finish(&arr_start_im);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: start_ion_mobility_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = stop_ion_mobility_builder.Finish(&arr_stop_im);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: stop_ion_mobility_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = intensities_builder.Finish(&arr_intensities);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: intensities_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = additional_intensities_builder.Finish(&arr_additional_intensities);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: additional_intensities_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = pg_accessions_builder.Finish(&arr_pg_acc);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: pg_accessions_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = anchor_protein_builder.Finish(&arr_anchor);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: anchor_protein_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = unique_builder.Finish(&arr_unique);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: unique_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = pg_qvalue_builder.Finish(&arr_pg_qval);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: pg_qvalue_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = gg_accessions_builder.Finish(&arr_gg_acc);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: gg_accessions_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = gg_names_builder.Finish(&arr_gg_names);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: gg_names_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = scan_ref_file_builder.Finish(&arr_scan_ref_file);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: scan_ref_file_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = rt_start_builder.Finish(&arr_rt_start);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: rt_start_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = rt_stop_builder.Finish(&arr_rt_stop);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: rt_stop_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = quality_builder.Finish(&arr_quality);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: quality_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = score_builder.Finish(&arr_score);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: score_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = score_type_builder.Finish(&arr_score_type);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: score_type_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = spectrum_reference_builder.Finish(&arr_spectrum_ref);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: spectrum_reference_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }
  status = feature_metavalues_builder.Finish(&arr_feature_mvs);
  if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: feature_metavalues_builder Finish failed: " << status.ToString() << std::endl; return nullptr; }

  // Build schema from registry
  auto schema = ConsensusFeatureExportSchema::schema();

  auto table = arrow::Table::Make(schema, {
    arr_sequence, arr_peptidoform, arr_modifications,
    arr_charge, arr_calc_mz, arr_obs_mz, arr_rt,
    arr_pep, arr_is_decoy, arr_additional_scores,
    arr_predicted_rt, arr_ref_file, arr_cv_params,
    arr_scan, arr_ion_mobility, arr_start_im, arr_stop_im,
    arr_intensities, arr_additional_intensities,
    arr_pg_acc, arr_anchor, arr_unique, arr_pg_qval,
    arr_gg_acc, arr_gg_names, arr_scan_ref_file,
    arr_rt_start, arr_rt_stop,
    arr_quality, arr_score, arr_score_type, arr_spectrum_ref,
    arr_feature_mvs
  });

  // Validate table against registry schema (strict — write path must match exactly)
  auto validation = ArrowSchemaValidation::validate(table, ConsensusFeatureExportSchema::schema());
  if (!validation.valid)
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowExport: Schema validation failed: " << validation.toString() << std::endl;
    return nullptr;
  }

  return table;
}


bool ConsensusMapArrowExport::exportToParquet(
  const ConsensusMap& cmap,
  const String& filename,
  const ParquetWriteConfig& config)
{
  auto table = exportToArrow(cmap);
  if (!table)
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowExport: Failed to create Arrow table" << std::endl;
    return false;
  }

  // Open output file
  auto result = arrow::io::FileOutputStream::Open(filename);
  if (!result.ok())
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowExport: Failed to open file: " << filename << std::endl;
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
  auto status = parquet::arrow::WriteTable(
    *table, arrow::default_memory_pool(), outfile,
    config.row_group_size, writer_props, arrow_props);

  if (!status.ok())
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowExport: Failed to write Parquet: "
                     << status.ToString() << std::endl;
    return false;
  }

  return true;
}

} // namespace OpenMS

#endif // WITH_PARQUET
