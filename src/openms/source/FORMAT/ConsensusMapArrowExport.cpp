// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ConsensusMapArrowExport.h>

#ifdef WITH_PARQUET

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/Scores.h>
#include <OpenMS/CHEMISTRY/ProForma.h>
#include <OpenMS/METADATA/SpectrumLookup.h>

#include <arrow/api.h>
#include <arrow/builder.h>
#include <arrow/io/file.h>
#include <parquet/arrow/writer.h>
#include <parquet/properties.h>

#include <unordered_map>
#include <unordered_set>

namespace OpenMS
{

namespace // anonymous
{
  /// Try to extract a scan number from a native ID string using known CV accessions.
  /// Returns -1 on failure.
  Int extractScan(const String& native_id)
  {
    static const std::vector<String> accessions = {
      "MS:1000768", "MS:1000769", "MS:1000771", "MS:1001508",
      "MS:1000774", "MS:1000777", "MS:1001530",
    };
    for (const auto& acc : accessions)
    {
      try
      {
        Int scan = SpectrumLookup::extractScanNumber(native_id, acc);
        if (scan >= 0) return scan;
      }
      catch (...)
      {
        // accession didn't match this native ID format, try next
      }
    }
    return -1;
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
    (void)spectrum_reference_builder.Append(spec_ref);

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
    (void)reference_file_builder.Append("");
    (void)scan_ref_file_builder.Append("");
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
    std::vector<std::string> protein_accs;
    if (best_hit)
    {
      for (const auto& ev : best_hit->getPeptideEvidences())
      {
        protein_accs.push_back(ev.getProteinAccession());
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

    (void)unique_builder.Append(protein_accs.size() == 1 ? 1 : 0);

    // === pg_global_qvalue ===
    bool found_qval = false;
    for (const auto& acc : protein_accs)
    {
      auto it = pg_qvalue_lookup.find(acc);
      if (it != pg_qvalue_lookup.end())
      {
        (void)pg_qvalue_builder.Append(it->second);
        found_qval = true;
        break;
      }
    }
    if (!found_qval)
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

  (void)sequence_builder.Finish(&arr_sequence);
  (void)peptidoform_builder.Finish(&arr_peptidoform);
  (void)modifications_builder.Finish(&arr_modifications);
  (void)precursor_charge_builder.Finish(&arr_charge);
  (void)calculated_mz_builder.Finish(&arr_calc_mz);
  (void)observed_mz_builder.Finish(&arr_obs_mz);
  (void)rt_builder.Finish(&arr_rt);
  (void)pep_builder.Finish(&arr_pep);
  (void)is_decoy_builder.Finish(&arr_is_decoy);
  (void)additional_scores_builder.Finish(&arr_additional_scores);
  (void)predicted_rt_builder.Finish(&arr_predicted_rt);
  (void)reference_file_builder.Finish(&arr_ref_file);
  (void)cv_params_builder.Finish(&arr_cv_params);
  (void)scan_builder.Finish(&arr_scan);
  (void)ion_mobility_builder.Finish(&arr_ion_mobility);
  (void)start_ion_mobility_builder.Finish(&arr_start_im);
  (void)stop_ion_mobility_builder.Finish(&arr_stop_im);
  (void)intensities_builder.Finish(&arr_intensities);
  (void)additional_intensities_builder.Finish(&arr_additional_intensities);
  (void)pg_accessions_builder.Finish(&arr_pg_acc);
  (void)anchor_protein_builder.Finish(&arr_anchor);
  (void)unique_builder.Finish(&arr_unique);
  (void)pg_qvalue_builder.Finish(&arr_pg_qval);
  (void)gg_accessions_builder.Finish(&arr_gg_acc);
  (void)gg_names_builder.Finish(&arr_gg_names);
  (void)scan_ref_file_builder.Finish(&arr_scan_ref_file);
  (void)rt_start_builder.Finish(&arr_rt_start);
  (void)rt_stop_builder.Finish(&arr_rt_stop);
  (void)quality_builder.Finish(&arr_quality);
  (void)score_builder.Finish(&arr_score);
  (void)score_type_builder.Finish(&arr_score_type);
  (void)spectrum_reference_builder.Finish(&arr_spectrum_ref);
  (void)feature_metavalues_builder.Finish(&arr_feature_mvs);

  // Build schema matching the Python to_feature_arrow() column order
  auto schema = arrow::schema({
    arrow::field("sequence", arrow::utf8()),
    arrow::field("peptidoform", arrow::utf8()),
    arrow::field("modifications", modifications_builder.type()),
    arrow::field("precursor_charge", arrow::int32()),
    arrow::field("calculated_mz", arrow::float32()),
    arrow::field("observed_mz", arrow::float32()),
    arrow::field("rt", arrow::float32()),
    arrow::field("posterior_error_probability", arrow::float64()),
    arrow::field("is_decoy", arrow::int32()),
    arrow::field("additional_scores", additional_scores_builder.type()),
    arrow::field("predicted_rt", arrow::float32()),
    arrow::field("reference_file_name", arrow::utf8()),
    arrow::field("cv_params", arrow::utf8()),
    arrow::field("scan", arrow::utf8()),
    arrow::field("ion_mobility", arrow::float32()),
    arrow::field("start_ion_mobility", arrow::float32()),
    arrow::field("stop_ion_mobility", arrow::float32()),
    arrow::field("intensities", intensities_builder.type()),
    arrow::field("additional_intensities", arrow::utf8()),
    arrow::field("pg_accessions", arrow::list(arrow::utf8())),
    arrow::field("anchor_protein", arrow::utf8()),
    arrow::field("unique", arrow::int32()),
    arrow::field("pg_global_qvalue", arrow::float64()),
    arrow::field("gg_accessions", arrow::list(arrow::utf8())),
    arrow::field("gg_names", arrow::list(arrow::utf8())),
    arrow::field("scan_reference_file_name", arrow::utf8()),
    arrow::field("rt_start", arrow::float32()),
    arrow::field("rt_stop", arrow::float32()),
    // OpenMS-specific columns
    arrow::field("quality", arrow::float32()),
    arrow::field("score", arrow::float64()),
    arrow::field("score_type", arrow::utf8()),
    arrow::field("spectrum_reference", arrow::utf8()),
    arrow::field("feature_metavalues", feature_metavalues_builder.type()),
  });

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
