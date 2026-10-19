// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/EntrapmentEvaluation.h>

#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <set>
#include <string_view>
#include <unordered_set>

namespace OpenMS
{
  namespace
  {
    double combinedFDP(Size accepted, Size entrapment, double te_ratio)
    {
      if (accepted == 0) return 0.0;
      return static_cast<double>(entrapment) / static_cast<double>(accepted) * (1.0 + te_ratio);
    }

    bool startsWith(const std::string& s, const std::string& prefix)
    {
      return s.size() >= prefix.size() && s.compare(0, prefix.size(), prefix) == 0;
    }

    /// whether @p acc carries @p decoy_string as a suffix, possibly followed by "_<n>" (DecoyDatabase with shuffle_decoy_ratio > 1)
    bool hasDecoySuffix(const std::string& acc, const std::string& decoy_string)
    {
      const size_t pos = acc.rfind(decoy_string);
      if (pos == std::string::npos) return false;
      const size_t rest = pos + decoy_string.size();
      if (rest == acc.size()) return true;
      if (acc[rest] != '_' || rest + 1 == acc.size()) return false;
      return std::all_of(acc.begin() + rest + 1, acc.end(), [](char c) { return std::isdigit(static_cast<unsigned char>(c)); });
    }

    /// I and L have the same mass: count them as one residue
    std::string toIL(std::string_view peptide)
    {
      std::string s(peptide);
      std::replace(s.begin(), s.end(), 'I', 'L');
      return s;
    }

    /// distinct digested peptides (I = L) of the target and of the entrapment proteins of @p db; decoys are skipped
    void digestDatabase(const std::vector<FASTAFile::FASTAEntry>& db, const std::string& entrapment_prefix,
                        const std::string& decoy_string, bool decoy_is_prefix, const std::string& enzyme,
                        Size missed_cleavages, Size min_length, Size max_length,
                        std::unordered_set<std::string>& target_peptides, std::unordered_set<std::string>* entrapment_peptides)
    {
      ProteaseDigestion digestion;
      digestion.setEnzyme(enzyme);
      digestion.setMissedCleavages(missed_cleavages);
      std::vector<std::string_view> peptides;
      for (const auto& entry : db)
      {
        const std::string& acc = entry.identifier;
        if (!decoy_string.empty() && (decoy_is_prefix ? startsWith(acc, decoy_string) : hasDecoySuffix(acc, decoy_string))) continue;
        const bool is_entrapment = !entrapment_prefix.empty() && startsWith(acc, entrapment_prefix);
        if (is_entrapment && entrapment_peptides == nullptr) continue;
        peptides.clear();
        digestion.digestUnmodified(entry.sequence, peptides, min_length, max_length);
        auto& bucket = is_entrapment ? *entrapment_peptides : target_peptides;
        for (const auto& p : peptides) bucket.emplace(toIL(p));
      }
    }
  }

  double EntrapmentEvaluation::Result::psmFDP() const
  {
    return combinedFDP(accepted_psms, entrapment_psms, te_ratio);
  }

  double EntrapmentEvaluation::Result::peptideFDP() const
  {
    return combinedFDP(accepted_peptides, entrapment_peptides, te_ratio);
  }

  EntrapmentEvaluation::Kind EntrapmentEvaluation::classify(const PeptideHit& hit, const std::string& entrapment_prefix)
  {
    if (hit.isDecoy()) return Kind::DECOY;
    const std::set<std::string> accessions = hit.extractProteinAccessionsSet();
    if (accessions.empty() || entrapment_prefix.empty()) return Kind::TARGET;
    for (const auto& acc : accessions)
    {
      if (!startsWith(acc, entrapment_prefix)) return Kind::TARGET;
    }
    return Kind::ENTRAPMENT;
  }

  double EntrapmentEvaluation::estimateTERatio(const std::vector<FASTAFile::FASTAEntry>& db,
                                               const std::string& entrapment_prefix,
                                               const std::string& decoy_string,
                                               bool decoy_is_prefix,
                                               const std::string& enzyme,
                                               Size missed_cleavages,
                                               Size min_length,
                                               Size max_length)
  {
    std::unordered_set<std::string> target_peptides, entrapment_peptides;
    digestDatabase(db, entrapment_prefix, decoy_string, decoy_is_prefix, enzyme, missed_cleavages, min_length, max_length,
                   target_peptides, &entrapment_peptides);
    Size n_entrapment = 0;
    for (const auto& p : entrapment_peptides)
    {
      if (!target_peptides.contains(p)) ++n_entrapment;
    }
    if (n_entrapment == 0) return 1.0;
    return static_cast<double>(target_peptides.size()) / static_cast<double>(n_entrapment);
  }

  std::unordered_set<std::string> EntrapmentEvaluation::targetPeptides(const std::vector<FASTAFile::FASTAEntry>& db,
                                                                      const std::string& entrapment_prefix,
                                                                      const std::string& decoy_string,
                                                                      bool decoy_is_prefix,
                                                                      const std::string& enzyme,
                                                                      Size missed_cleavages,
                                                                      Size min_length,
                                                                      Size max_length)
  {
    std::unordered_set<std::string> target_peptides;
    digestDatabase(db, entrapment_prefix, decoy_string, decoy_is_prefix, enzyme, missed_cleavages, min_length, max_length,
                   target_peptides, nullptr);
    return target_peptides;
  }

  double EntrapmentEvaluation::qValue(const PeptideIdentification& pep_id, const PeptideHit& hit)
  {
    if (pep_id.getScoreType().find("q-value") != std::string::npos) return hit.getScore();
    if (hit.metaValueExists("q-value")) return static_cast<double>(hit.getMetaValue("q-value"));
    if (hit.metaValueExists("MS:1001491")) return static_cast<double>(hit.getMetaValue("MS:1001491"));
    throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "No q-value for a peptide hit: the score type is '" + pep_id.getScoreType() +
      "' and the hit has no 'q-value' meta value. Run FalseDiscoveryRate or PercolatorAdapter first.");
  }

  std::map<std::string, EntrapmentEvaluation::Result> EntrapmentEvaluation::evaluate(const PeptideIdentificationList& pep_ids,
                                                                                    double q_threshold,
                                                                                    const std::string& entrapment_prefix,
                                                                                    double te_ratio,
                                                                                    const std::string& class_meta_value,
                                                                                    const std::unordered_set<std::string>* target_peptides)
  {
    std::map<std::string, Result> results;
    // distinct peptides per class, with "is entrapment" (a sequence counts as entrapment only if all its PSMs are)
    std::map<std::string, std::map<std::string, bool>> peptides;

    auto count = [&](const std::string& key, const std::string& sequence, bool is_entrapment)
    {
      Result& r = results[key];
      r.te_ratio = te_ratio;
      ++r.accepted_psms;
      if (is_entrapment) ++r.entrapment_psms;
      auto [it, inserted] = peptides[key].emplace(sequence, is_entrapment);
      if (!inserted) it->second = it->second && is_entrapment;
    };

    results["all"].te_ratio = te_ratio;
    for (const auto& pep_id : pep_ids)
    {
      if (pep_id.getHits().empty()) continue;
      const PeptideHit& hit = pep_id.getHits().front();
      const Kind kind = classify(hit, entrapment_prefix);
      if (kind == Kind::DECOY) continue;
      const double q = qValue(pep_id, hit);
      if (!std::isfinite(q) || q < 0.0 || q > q_threshold) continue; // also rejects NaN

      // an entrapment peptide that equals a target peptide up to I/L cannot be told apart from it
      const bool is_entrapment = (kind == Kind::ENTRAPMENT) &&
        (target_peptides == nullptr || !target_peptides->contains(toIL(hit.getSequence().toUnmodifiedString())));
      const std::string sequence = hit.getSequence().toString();
      count("all", sequence, is_entrapment);
      if (!class_meta_value.empty())
      {
        std::string cls = hit.metaValueExists(class_meta_value) ? hit.getMetaValue(class_meta_value).toString() : std::string("NA");
        if (cls == "all") cls = "class:all"; // "all" is the overall result
        count(cls, sequence, is_entrapment);
      }
    }

    for (auto& [key, r] : results)
    {
      const auto& seqs = peptides[key];
      r.accepted_peptides = seqs.size();
      r.entrapment_peptides = 0;
      for (const auto& [seq, is_entrapment] : seqs)
      {
        if (is_entrapment) ++r.entrapment_peptides;
      }
    }
    return results;
  }
} // namespace OpenMS
