// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ConsensusMapArrowExport.h>

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
  arrow::StringBuilder run_file_name_builder, id_run_file_name_builder;
  arrow::Int16Builder charge_builder, missed_cleavages_builder;
  arrow::BooleanBuilder is_decoy_builder, unique_builder;
  arrow::FloatBuilder calculated_mz_builder, observed_mz_builder, rt_builder;
  arrow::FloatBuilder mass_error_ppm_builder;
  arrow::FloatBuilder rt_start_builder, rt_stop_builder;
  arrow::FloatBuilder ion_mobility_builder, predicted_rt_builder;
  arrow::FloatBuilder im_start_builder, im_stop_builder;
  arrow::DoubleBuilder pep_builder, pg_qvalue_builder;

  // -- List<utf8> builders for gg_accessions, gg_names --
  auto gg_acc_vb = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder gg_accessions_builder(arrow::default_memory_pool(), gg_acc_vb);

  auto gg_names_vb = std::make_shared<arrow::StringBuilder>();
  arrow::ListBuilder gg_names_builder(arrow::default_memory_pool(), gg_names_vb);

  // -- scan: list<int32> --
  auto scan_val_b = std::make_shared<arrow::Int32Builder>();
  arrow::ListBuilder scan_builder(arrow::default_memory_pool(), scan_val_b);

  // -- Modifications: list<struct{name, accession, positions: list<struct{position: int32, amino_acid, scores: list<struct>}>}> --
  // scores within position
  auto mod_score_name_b = std::make_shared<arrow::StringBuilder>();
  auto mod_score_value_b = std::make_shared<arrow::DoubleBuilder>();
  auto mod_score_hb_b = std::make_shared<arrow::BooleanBuilder>();
  auto mod_score_struct_type = arrow::struct_({
    arrow::field("score_name", arrow::utf8(), /*nullable=*/false),
    arrow::field("score_value", arrow::float64(), /*nullable=*/false),
    arrow::field("higher_better", arrow::boolean())
  });
  auto mod_score_struct_b = std::make_shared<arrow::StructBuilder>(
    mod_score_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{mod_score_name_b, mod_score_value_b, mod_score_hb_b});
  auto mod_scores_list_b = std::make_shared<arrow::ListBuilder>(arrow::default_memory_pool(), mod_score_struct_b);

  auto pos_position_b = std::make_shared<arrow::Int32Builder>();
  auto pos_amino_acid_b = std::make_shared<arrow::StringBuilder>();
  auto position_struct_type = arrow::struct_({
    arrow::field("position", arrow::int32(), /*nullable=*/false),
    arrow::field("amino_acid", arrow::utf8()),
    arrow::field("scores", arrow::list(mod_score_struct_type))
  });
  auto position_struct_b = std::make_shared<arrow::StructBuilder>(
    position_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{pos_position_b, pos_amino_acid_b, mod_scores_list_b});
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

  // -- additional_scores: list<struct{score_name, score_value, higher_better}> --
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

  // -- cv_params: list<struct{cv_name, cv_value}> --
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

  // -- intensities: list<struct{label, intensity}> --
  auto intensity_label_b = std::make_shared<arrow::StringBuilder>();
  auto intensity_val_b = std::make_shared<arrow::FloatBuilder>();
  auto intensity_struct_type = arrow::struct_({
    arrow::field("label", arrow::utf8(), /*nullable=*/false),
    arrow::field("intensity", arrow::float32(), /*nullable=*/false)
  });
  auto intensity_struct_b = std::make_shared<arrow::StructBuilder>(
    intensity_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{intensity_label_b, intensity_val_b});
  arrow::ListBuilder intensities_builder(arrow::default_memory_pool(), intensity_struct_b);

  // -- additional_intensities: list<struct{label, intensities: list<struct{intensity_name, intensity_value}>}> --
  auto ai_iname_b = std::make_shared<arrow::StringBuilder>();
  auto ai_ival_b = std::make_shared<arrow::FloatBuilder>();
  auto ai_entry_type = arrow::struct_({
    arrow::field("intensity_name", arrow::utf8(), /*nullable=*/false),
    arrow::field("intensity_value", arrow::float32(), /*nullable=*/false)
  });
  auto ai_entry_b = std::make_shared<arrow::StructBuilder>(
    ai_entry_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{ai_iname_b, ai_ival_b});
  auto ai_list_b = std::make_shared<arrow::ListBuilder>(arrow::default_memory_pool(), ai_entry_b);
  auto ai_label_b = std::make_shared<arrow::StringBuilder>();
  auto ai_outer_type = arrow::struct_({
    arrow::field("label", arrow::utf8(), /*nullable=*/false),
    arrow::field("intensities", arrow::list(ai_entry_type), /*nullable=*/false)
  });
  auto ai_outer_b = std::make_shared<arrow::StructBuilder>(
    ai_outer_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{ai_label_b, ai_list_b});
  arrow::ListBuilder additional_intensities_builder(arrow::default_memory_pool(), ai_outer_b);

  // -- pg_accessions: list<struct{accession, start, end, pre, post}> --
  auto pga_acc_b = std::make_shared<arrow::StringBuilder>();
  auto pga_start_b = std::make_shared<arrow::Int32Builder>();
  auto pga_end_b = std::make_shared<arrow::Int32Builder>();
  auto pga_pre_b = std::make_shared<arrow::StringBuilder>();
  auto pga_post_b = std::make_shared<arrow::StringBuilder>();
  auto pga_struct_type = arrow::struct_({
    arrow::field("accession", arrow::utf8(), /*nullable=*/false),
    arrow::field("start", arrow::int32()),
    arrow::field("end", arrow::int32()),
    arrow::field("pre", arrow::utf8()),
    arrow::field("post", arrow::utf8())
  });
  auto pga_struct_b = std::make_shared<arrow::StructBuilder>(
    pga_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{pga_acc_b, pga_start_b, pga_end_b, pga_pre_b, pga_post_b});
  arrow::ListBuilder pg_accessions_builder(arrow::default_memory_pool(), pga_struct_b);

  // -- pg_positions: list<struct{protein_accession, start, end}> --
  auto pgp_acc_b = std::make_shared<arrow::StringBuilder>();
  auto pgp_start_b = std::make_shared<arrow::Int32Builder>();
  auto pgp_end_b = std::make_shared<arrow::Int32Builder>();
  auto pgp_struct_type = arrow::struct_({
    arrow::field("protein_accession", arrow::utf8(), /*nullable=*/false),
    arrow::field("start", arrow::int32(), /*nullable=*/false),
    arrow::field("end", arrow::int32(), /*nullable=*/false)
  });
  auto pgp_struct_b = std::make_shared<arrow::StructBuilder>(
    pgp_struct_type, arrow::default_memory_pool(),
    std::vector<std::shared_ptr<arrow::ArrayBuilder>>{pgp_acc_b, pgp_start_b, pgp_end_b});
  arrow::ListBuilder pg_positions_builder(arrow::default_memory_pool(), pgp_struct_b);

  arrow::Status status;
  const Size num_features = cmap.size();

  // Reserve capacity for simple column builders
  #define RESERVE_OR_FAIL(builder) \
    status = builder.Reserve(num_features); \
    if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: " #builder " Reserve failed: " << status.ToString() << std::endl; return nullptr; }

  RESERVE_OR_FAIL(sequence_builder)
  RESERVE_OR_FAIL(peptidoform_builder)
  RESERVE_OR_FAIL(anchor_protein_builder)
  RESERVE_OR_FAIL(run_file_name_builder)
  RESERVE_OR_FAIL(id_run_file_name_builder)
  RESERVE_OR_FAIL(charge_builder)
  RESERVE_OR_FAIL(missed_cleavages_builder)
  RESERVE_OR_FAIL(is_decoy_builder)
  RESERVE_OR_FAIL(unique_builder)
  RESERVE_OR_FAIL(calculated_mz_builder)
  RESERVE_OR_FAIL(observed_mz_builder)
  RESERVE_OR_FAIL(mass_error_ppm_builder)
  RESERVE_OR_FAIL(rt_builder)
  RESERVE_OR_FAIL(rt_start_builder)
  RESERVE_OR_FAIL(rt_stop_builder)
  RESERVE_OR_FAIL(ion_mobility_builder)
  RESERVE_OR_FAIL(predicted_rt_builder)
  RESERVE_OR_FAIL(im_start_builder)
  RESERVE_OR_FAIL(im_stop_builder)
  RESERVE_OR_FAIL(pep_builder)
  RESERVE_OR_FAIL(pg_qvalue_builder)

  #undef RESERVE_OR_FAIL

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
        float calc_mz = static_cast<float>(seq.getMZ(charge));
        (void)calculated_mz_builder.Append(calc_mz);
        float obs_mz = static_cast<float>(cf.getMZ());
        // mass_error_ppm = (observed - calculated) / calculated * 1e6
        if (calc_mz > 0)
        {
          (void)mass_error_ppm_builder.Append((obs_mz - calc_mz) / calc_mz * 1e6f);
        }
        else
        {
          (void)mass_error_ppm_builder.AppendNull();
        }
      }
      else
      {
        (void)calculated_mz_builder.Append(static_cast<float>(cf.getMZ()));
        (void)mass_error_ppm_builder.AppendNull();
      }

      // === modifications ===
      (void)modifications_builder.Append();
      if (seq.isModified())
      {
        // Collect modifications keyed by name to group positions
        // name -> (accession, [(1-based-position, amino_acid_char)])
        std::map<std::string, std::pair<std::string, std::vector<std::pair<int, std::string>>>> mod_map;

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
          mod_map[name].second.push_back({0, ""});
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
            if (mod_map.find(name) == mod_map.end())
            {
              mod_map[name] = {acc_str, {}};
            }
            mod_map[name].second.push_back({static_cast<int>(pos + 1), std::string(residue.getOneLetterCode())});
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
          if (mod_map.find(name) == mod_map.end())
          {
            mod_map[name] = {acc_str, {}};
          }
          mod_map[name].second.push_back({static_cast<int>(seq.size() + 1), ""});
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
          for (const auto& [pos_idx, aa] : acc_positions.second)
          {
            (void)position_struct_b->Append();
            (void)pos_position_b->Append(pos_idx);
            if (aa.empty())
            {
              (void)pos_amino_acid_b->AppendNull();
            }
            else
            {
              (void)pos_amino_acid_b->Append(aa);
            }
            (void)mod_scores_list_b->Append(); // empty scores list
          }
        }
      }
      // else: empty modifications list (Append() already called)
    }
    else
    {
      (void)sequence_builder.Append("");
      (void)peptidoform_builder.Append("");
      (void)calculated_mz_builder.Append(static_cast<float>(cf.getMZ()));
      (void)mass_error_ppm_builder.AppendNull();
      // empty modifications list
      (void)modifications_builder.Append();
    }

    (void)charge_builder.Append(static_cast<int16_t>(cf.getCharge()));
    (void)observed_mz_builder.Append(static_cast<float>(cf.getMZ()));
    (void)rt_builder.Append(static_cast<float>(cf.getRT()));

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
      (void)is_decoy_builder.Append(td.substr(0, 5) == "decoy");
    }
    else
    {
      (void)is_decoy_builder.Append(false);
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

    // === run_file_name ===
    // Derive from the first FeatureHandle's map index -> column header filename
    {
      bool found_run = false;
      for (const auto& fh : cf.getFeatures())
      {
        auto it = column_headers.find(fh.getMapIndex());
        if (it != column_headers.end())
        {
          (void)run_file_name_builder.Append(it->second.filename);
          found_run = true;
          break;
        }
      }
      if (!found_run)
      {
        (void)run_file_name_builder.Append("");
      }
    }

    // === scan (list<int32>) ===
    {
      std::string spec_ref;
      if (!pep_ids.empty())
      {
        // Prefer the dedicated member, fall back to metavalue
        spec_ref = pep_ids[0].getSpectrumReference();
        if (spec_ref.empty() && pep_ids[0].metaValueExists("spectrum_reference"))
        {
          spec_ref = pep_ids[0].getMetaValue("spectrum_reference").toString();
        }
      }

      (void)scan_builder.Append();
      if (!spec_ref.empty())
      {
        Int scan_num = extractScan(spec_ref);
        if (scan_num >= 0)
        {
          (void)scan_val_b->Append(scan_num);
        }
      }
    }

    // === cv_params (empty list for now) ===
    (void)cv_params_builder.AppendNull();

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
      (void)im_start_builder.Append(static_cast<float>(static_cast<double>(cf.getMetaValue("start_ion_mobility"))));
    }
    else
    {
      (void)im_start_builder.AppendNull();
    }

    if (cf.metaValueExists("stop_ion_mobility"))
    {
      (void)im_stop_builder.Append(static_cast<float>(static_cast<double>(cf.getMetaValue("stop_ion_mobility"))));
    }
    else
    {
      (void)im_stop_builder.AppendNull();
    }

    // === missed_cleavages ===
    if (best_hit && best_hit->metaValueExists("missed_cleavages"))
    {
      (void)missed_cleavages_builder.Append(static_cast<int16_t>(static_cast<int>(best_hit->getMetaValue("missed_cleavages"))));
    }
    else
    {
      (void)missed_cleavages_builder.AppendNull();
    }

    // === additional_intensities (empty list for now) ===
    (void)additional_intensities_builder.AppendNull();

    // === id_run_file_name ===
    (void)id_run_file_name_builder.AppendNull();

    // === Protein accessions (now list<struct>) ===
    // Collect ALL peptide evidences (including repeated accessions with different positions)
    // and separately track unique accessions for anchor_protein/unique
    std::vector<std::string> protein_accs;
    std::unordered_set<std::string> seen_accs;
    struct EvidenceInfo { std::string acc; int start; int end; char pre; char post; };
    std::vector<EvidenceInfo> evidences;
    if (best_hit)
    {
      for (const auto& ev : best_hit->getPeptideEvidences())
      {
        const std::string& acc = ev.getProteinAccession();
        // Track unique accessions for anchor_protein / unique
        if (seen_accs.insert(acc).second)
        {
          protein_accs.push_back(acc);
        }
        // Emit every evidence (including repeated accessions at different positions)
        EvidenceInfo ei;
        ei.acc = acc;
        ei.start = ev.getStart();
        ei.end = ev.getEnd();
        ei.pre = ev.getAABefore();
        ei.post = ev.getAAAfter();
        evidences.push_back(ei);
      }
    }

    (void)pg_accessions_builder.Append();
    for (const auto& ei : evidences)
    {
      (void)pga_struct_b->Append();
      (void)pga_acc_b->Append(ei.acc);
      if (ei.start != PeptideEvidence::UNKNOWN_POSITION)
      {
        (void)pga_start_b->Append(ei.start);
      }
      else
      {
        (void)pga_start_b->AppendNull();
      }
      if (ei.end != PeptideEvidence::UNKNOWN_POSITION)
      {
        (void)pga_end_b->Append(ei.end);
      }
      else
      {
        (void)pga_end_b->AppendNull();
      }
      if (ei.pre != PeptideEvidence::UNKNOWN_AA && ei.pre != 0)
      {
        (void)pga_pre_b->Append(std::string(1, ei.pre));
      }
      else
      {
        (void)pga_pre_b->AppendNull();
      }
      if (ei.post != PeptideEvidence::UNKNOWN_AA && ei.post != 0)
      {
        (void)pga_post_b->Append(std::string(1, ei.post));
      }
      else
      {
        (void)pga_post_b->AppendNull();
      }
    }

    // anchor_protein
    if (protein_accs.empty())
    {
      (void)anchor_protein_builder.Append("");
    }
    else
    {
      (void)anchor_protein_builder.Append(protein_accs[0]);
    }

    // unique: true if exactly one unique protein accession
    (void)unique_builder.Append(protein_accs.size() == 1);

    // === pg_global_qvalue ===
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

    // === pg_positions (empty for now — no positional data available at consensus level) ===
    (void)pg_positions_builder.AppendNull();

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

    // === intensities: list<struct{label, intensity}> ===
    (void)intensities_builder.Append();
    for (const auto& fh : cf.getFeatures())
    {
      auto it = column_headers.find(fh.getMapIndex());
      if (it != column_headers.end())
      {
        (void)intensity_struct_b->Append();
        (void)intensity_label_b->Append(it->second.filename);
        (void)intensity_val_b->Append(static_cast<float>(fh.getIntensity()));
      }
    }
  } // end feature loop

  // Finalize all arrays
  #define FINISH_OR_FAIL(builder, arr) \
    status = builder.Finish(&arr); \
    if (!status.ok()) { OPENMS_LOG_ERROR << "ConsensusMapArrowExport: " #builder " Finish failed: " << status.ToString() << std::endl; return nullptr; }

  std::shared_ptr<arrow::Array> arr_sequence, arr_peptidoform, arr_modifications;
  std::shared_ptr<arrow::Array> arr_charge, arr_pep, arr_is_decoy;
  std::shared_ptr<arrow::Array> arr_calc_mz, arr_obs_mz, arr_mass_error_ppm;
  std::shared_ptr<arrow::Array> arr_additional_scores, arr_predicted_rt;
  std::shared_ptr<arrow::Array> arr_run_file, arr_cv_params, arr_scan;
  std::shared_ptr<arrow::Array> arr_rt, arr_ion_mobility, arr_missed_cleavages;
  std::shared_ptr<arrow::Array> arr_intensities, arr_additional_intensities;
  std::shared_ptr<arrow::Array> arr_pg_acc, arr_anchor, arr_unique, arr_pg_qval;
  std::shared_ptr<arrow::Array> arr_pg_positions;
  std::shared_ptr<arrow::Array> arr_im_start, arr_im_stop;
  std::shared_ptr<arrow::Array> arr_gg_acc, arr_gg_names;
  std::shared_ptr<arrow::Array> arr_id_run_file;
  std::shared_ptr<arrow::Array> arr_rt_start, arr_rt_stop;

  FINISH_OR_FAIL(sequence_builder, arr_sequence)
  FINISH_OR_FAIL(peptidoform_builder, arr_peptidoform)
  FINISH_OR_FAIL(modifications_builder, arr_modifications)
  FINISH_OR_FAIL(charge_builder, arr_charge)
  FINISH_OR_FAIL(pep_builder, arr_pep)
  FINISH_OR_FAIL(is_decoy_builder, arr_is_decoy)
  FINISH_OR_FAIL(calculated_mz_builder, arr_calc_mz)
  FINISH_OR_FAIL(observed_mz_builder, arr_obs_mz)
  FINISH_OR_FAIL(mass_error_ppm_builder, arr_mass_error_ppm)
  FINISH_OR_FAIL(additional_scores_builder, arr_additional_scores)
  FINISH_OR_FAIL(predicted_rt_builder, arr_predicted_rt)
  FINISH_OR_FAIL(run_file_name_builder, arr_run_file)
  FINISH_OR_FAIL(cv_params_builder, arr_cv_params)
  FINISH_OR_FAIL(scan_builder, arr_scan)
  FINISH_OR_FAIL(rt_builder, arr_rt)
  FINISH_OR_FAIL(ion_mobility_builder, arr_ion_mobility)
  FINISH_OR_FAIL(missed_cleavages_builder, arr_missed_cleavages)
  FINISH_OR_FAIL(intensities_builder, arr_intensities)
  FINISH_OR_FAIL(additional_intensities_builder, arr_additional_intensities)
  FINISH_OR_FAIL(pg_accessions_builder, arr_pg_acc)
  FINISH_OR_FAIL(anchor_protein_builder, arr_anchor)
  FINISH_OR_FAIL(unique_builder, arr_unique)
  FINISH_OR_FAIL(pg_qvalue_builder, arr_pg_qval)
  FINISH_OR_FAIL(pg_positions_builder, arr_pg_positions)
  FINISH_OR_FAIL(im_start_builder, arr_im_start)
  FINISH_OR_FAIL(im_stop_builder, arr_im_stop)
  FINISH_OR_FAIL(gg_accessions_builder, arr_gg_acc)
  FINISH_OR_FAIL(gg_names_builder, arr_gg_names)
  FINISH_OR_FAIL(id_run_file_name_builder, arr_id_run_file)
  FINISH_OR_FAIL(rt_start_builder, arr_rt_start)
  FINISH_OR_FAIL(rt_stop_builder, arr_rt_stop)

  #undef FINISH_OR_FAIL

  // Build schema from registry (QPX format)
  auto schema = QPXFeatureSchema::schema();

  auto table = arrow::Table::Make(schema, {
    arr_sequence, arr_peptidoform, arr_modifications,
    arr_charge, arr_pep, arr_is_decoy,
    arr_calc_mz, arr_obs_mz, arr_mass_error_ppm,
    arr_additional_scores, arr_predicted_rt,
    arr_run_file, arr_cv_params, arr_scan,
    arr_rt, arr_ion_mobility, arr_missed_cleavages,
    arr_intensities, arr_additional_intensities,
    arr_pg_acc, arr_anchor, arr_unique, arr_pg_qval,
    arr_pg_positions, arr_im_start, arr_im_stop,
    arr_gg_acc, arr_gg_names, arr_id_run_file,
    arr_rt_start, arr_rt_stop,
  });

  // Validate table against registry schema (strict — write path must match exactly)
  auto validation = ArrowSchemaValidation::validate(table, QPXFeatureSchema::schema());
  if (!validation.valid)
  {
    OPENMS_LOG_ERROR << "ConsensusMapArrowExport: Schema validation failed: " << validation.toString() << "\n";
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
