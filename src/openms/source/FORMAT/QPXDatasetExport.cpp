// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/QPXDatasetExport.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>

#include <cctype>
#include <map>
#include <optional>
#include <set>

namespace OpenMS
{
  namespace
  {
    /// Whether a higher value of the score @p name (a QPXDatasetExport::scoreName()) is better; none if not known
    std::optional<bool> higherBetter(const std::string& name)
    {
      if (name == "global_qvalue" || name == "percolator_qvalue" || name == "percolator_pep") { return false; }
      if (name == "percolator_score" || name == "ln_hyperscore") { return true; }
      return std::nullopt;
    }
  }

  std::string QPXDatasetExport::scoreName(const std::string& name)
  {
    static const std::map<std::string, std::string> known = {
      {"q-value", "global_qvalue"},
      {"MS:1001491", "percolator_qvalue"}, {"percolator_q_value", "percolator_qvalue"},
      {"MS:1001492", "percolator_score"}, {"svm", "percolator_score"},
      {"MS:1001493", "percolator_pep"},
      {"ln(hyperscore)", "ln_hyperscore"}, {"ln(hyperscore)_score", "ln_hyperscore"},
      {"Posterior Error Probability", ""}, {"Posterior Error Probability_score", ""},
      {"IsotopeError", "isotope_error"}};
    const auto found = known.find(name);
    if (found != known.end()) { return found->second; }
    std::string out;
    for (Size i = 0; i < name.size(); ++i)
    {
      const unsigned char c = static_cast<unsigned char>(name[i]);
      const unsigned char before = i > 0 ? static_cast<unsigned char>(name[i - 1]) : 0;
      if (std::isupper(c) && (std::islower(before) || std::isdigit(before))) { out += '_'; }
      if (std::isalnum(c)) { out += static_cast<char>(std::tolower(c)); }
      else if (!out.empty() && out.back() != '_') { out += '_'; }
    }
    while (!out.empty() && out.back() == '_') { out.pop_back(); }
    return out;
  }

  std::vector<QPXDatasetExport::Column> QPXDatasetExport::psmColumns(const PeptideIdentificationList& peptide_ids,
                                                                     const std::vector<bool>& has_spectrum,
                                                                     bool isotope_observed_minus_theoretical)
  {
    if (!has_spectrum.empty() && has_spectrum.size() != peptide_ids.size())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "psmColumns: has_spectrum has " + std::to_string(has_spectrum.size()) + " entries for "
        + std::to_string(peptide_ids.size()) + " peptide identifications");
    }

    auto column = [](const std::string& name, Column::Type type) { Column c; c.name = name; c.type = type; return c; };
    Column scores = column(QPXPSMSchema::ADDITIONAL_SCORES, Column::Type::SCORE_LIST);
    Column reference = column("spectrum_reference", Column::Type::UTF8);
    Column isotope = column("isotope_error", Column::Type::INT16);
    Column mappings = column("protein_mappings", Column::Type::STRUCT_LIST);
    mappings.struct_fields = {{"accession", Column::Type::UTF8, false}, {"start", Column::Type::INT32}, {"end", Column::Type::INT32},
                              {"aa_before", Column::Type::UTF8}, {"aa_after", Column::Type::UTF8}};

    // The isotope error as QPX has it (bigbio/qpx#355: the 13C peaks from the monoisotopic to the selected peak, i.e.
    // observed minus theoretical), in the column and in 'additional_scores'
    const int isotope_sign = isotope_observed_minus_theoretical ? 1 : -1;
    // meta values that are no scores (columns of their own, or attributes of the spectrum)
    static const std::set<std::string> not_scores = {"target_decoy", "predicted_RT", "predicted_rt", "ion_mobility", "IM",
                                                     "FAIMS_CV", "scan", "reference_file_name", "id_merge_index"};
    auto residue = [](char aa) -> std::optional<std::string>
    {
      if (aa == PeptideEvidence::N_TERMINAL_AA || aa == PeptideEvidence::C_TERMINAL_AA) { return std::string("-"); }
      if (aa == PeptideEvidence::UNKNOWN_AA) { return std::nullopt; }
      return std::string(1, aa);
    };
    auto position = [](Int p) -> std::optional<Int32>
    {
      if (p == PeptideEvidence::UNKNOWN_POSITION || p < 0) { return std::nullopt; }
      return static_cast<Int32>(p + 1); // 1-based
    };

    for (Size i = 0; i < peptide_ids.size(); ++i)
    {
      const PeptideIdentification& id = peptide_ids[i];
      if (id.getHits().empty()) { continue; } // no QPX row
      const PeptideHit& hit = id.getHits()[0];

      std::vector<Column::Score> row_scores;
      std::set<std::string> seen;
      auto add = [&](const std::string& key, double value, std::optional<bool> higher_better)
      {
        const std::string name = scoreName(key);
        if (name.empty() || !seen.insert(name).second) { return; }
        const std::optional<bool> known = higherBetter(name);
        if (name == "isotope_error" && value != 0.0) { value *= isotope_sign; } // (no -0)
        row_scores.push_back({name, value, known ? known : higher_better});
      };
      add(id.getScoreType(), hit.getScore(), id.isHigherScoreBetter());
      std::vector<std::string> keys;
      hit.getKeys(keys);
      for (const std::string& key : keys)
      {
        if (not_scores.count(key) > 0) { continue; }
        const DataValue& value = hit.getMetaValue(key);
        if (value.valueType() == DataValue::INT_VALUE || value.valueType() == DataValue::DOUBLE_VALUE)
        {
          add(key, static_cast<double>(value), std::nullopt);
        }
      }
      scores.score_lists.push_back(std::move(row_scores));

      // (null: the mz view has no row of its spectrum, or no native ID)
      const bool spectrum = has_spectrum.empty() || has_spectrum[i];
      reference.utf8_values.push_back(spectrum ? id.getSpectrumReference() : std::string());
      reference.null_rows.push_back(!spectrum || id.getSpectrumReference().empty());

      const std::string isotope_key = hit.metaValueExists("isotope_error") ? "isotope_error" : "IsotopeError";
      isotope.int64_values.push_back(hit.metaValueExists(isotope_key) ? isotope_sign * static_cast<int>(hit.getMetaValue(isotope_key)) : 0);
      isotope.null_rows.push_back(!hit.metaValueExists(isotope_key));

      std::vector<Column::StructEntry> row_mappings;
      for (const PeptideEvidence& evidence : hit.getPeptideEvidences())
      {
        row_mappings.push_back({evidence.getProteinAccession(), position(evidence.getStart()), position(evidence.getEnd()),
                                residue(evidence.getAABefore()), residue(evidence.getAAAfter())});
      }
      mappings.struct_lists.push_back(std::move(row_mappings));
    }
    return {std::move(scores), std::move(reference), std::move(isotope), std::move(mappings)};
  }
} // namespace OpenMS
