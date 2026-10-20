// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/EntrapmentEvaluation.h>

#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <set>
#include <string_view>
#include <unordered_set>

namespace OpenMS
{
  namespace
  {
    double combinedFDP_(Size accepted, Size entrapment, double te_ratio)
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

  double EntrapmentEvaluation::Result::lowerBoundFDP() const
  {
    if (accepted == 0) return 0.0;
    return static_cast<double>(entrapment) / static_cast<double>(accepted);
  }

  double EntrapmentEvaluation::Result::combinedFDP() const
  {
    return combinedFDP_(accepted, entrapment, te_ratio);
  }

  double EntrapmentEvaluation::Result::pairedFDP() const
  {
    if (!paired) return std::numeric_limits<double>::quiet_NaN();
    if (accepted == 0) return 0.0;
    return (static_cast<double>(entrapment) + static_cast<double>(entrapment_target_not_accepted) +
            2.0 * static_cast<double>(entrapment_above_target)) / static_cast<double>(accepted);
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

  EntrapmentEvaluation::PeptidePairs EntrapmentEvaluation::pairPeptides(const std::vector<FASTAFile::FASTAEntry>& db,
                                                                       const std::string& entrapment_prefix,
                                                                       const std::string& decoy_string,
                                                                       bool decoy_is_prefix,
                                                                       const std::string& enzyme,
                                                                       Size missed_cleavages,
                                                                       Size min_length,
                                                                       Size max_length,
                                                                       Size* mismatched)
  {
    PeptidePairs pairs;
    if (mismatched != nullptr) *mismatched = 0;
    if (entrapment_prefix.empty()) return pairs;

    auto isDecoy = [&](const std::string& acc)
    {
      return !decoy_string.empty() && (decoy_is_prefix ? startsWith(acc, decoy_string) : hasDecoySuffix(acc, decoy_string));
    };
    std::unordered_map<std::string, const std::string*> target_sequences;
    for (const auto& entry : db)
    {
      if (!isDecoy(entry.identifier) && !startsWith(entry.identifier, entrapment_prefix)) target_sequences.emplace(entry.identifier, &entry.sequence);
    }

    ProteaseDigestion digestion;
    digestion.setEnzyme(enzyme);
    digestion.setMissedCleavages(missed_cleavages);
    std::vector<std::string_view> peptides;
    for (const auto& entry : db)
    {
      if (isDecoy(entry.identifier) || !startsWith(entry.identifier, entrapment_prefix)) continue;
      const auto target_it = target_sequences.find(entry.identifier.substr(entrapment_prefix.size()));
      if (target_it == target_sequences.end()) continue;
      const std::string& target = *target_it->second;
      const std::string& entrapment = entry.sequence;
      if (target.size() != entrapment.size())
      {
        if (mismatched != nullptr) ++*mismatched;
        continue;
      }
      peptides.clear();
      digestion.digestUnmodified(target, peptides, min_length, max_length);
      for (const auto& peptide : peptides)
      {
        // the views point into the target sequence
        const Size offset = static_cast<Size>(peptide.data() - target.data());
        std::string target_peptide(peptide);
        std::string entrapment_peptide = entrapment.substr(offset, peptide.size());
        std::string sorted_target = target_peptide, sorted_entrapment = entrapment_peptide;
        std::sort(sorted_target.begin(), sorted_target.end());
        std::sort(sorted_entrapment.begin(), sorted_entrapment.end());
        if (sorted_target != sorted_entrapment)
        {
          if (mismatched != nullptr) ++*mismatched;
          continue;
        }
        std::string entrapment_il = toIL(entrapment_peptide);
        std::string target_il = toIL(target_peptide);
        if (entrapment_il == target_il) continue;
        pairs.emplace(std::move(entrapment_il), std::move(target_il));
      }
    }
    return pairs;
  }

  std::map<std::string, EntrapmentEvaluation::Result> EntrapmentEvaluation::evaluate(const PeptideIdentificationList& pep_ids,
                                                                                    double q_threshold,
                                                                                    const std::string& entrapment_prefix,
                                                                                    double te_ratio,
                                                                                    Level level,
                                                                                    const std::string& class_meta_value,
                                                                                    const std::unordered_set<std::string>* target_peptides,
                                                                                    const PeptidePairs* pairs,
                                                                                    const std::string& score_meta_value)
  {
    const bool paired = (level == Level::PEPTIDE) && (pairs != nullptr);
    std::map<std::string, Result> results;
    auto result = [&](const std::string& key) -> Result&
    {
      Result& r = results[key];
      r.te_ratio = te_ratio;
      r.paired = paired;
      return r;
    };
    result("all");

    auto classOf = [&](const PeptideHit& hit)
    {
      std::string cls = hit.metaValueExists(class_meta_value) ? hit.getMetaValue(class_meta_value).toString() : std::string("NA");
      if (cls.empty()) cls = "NA";
      if (cls == "all") cls = "class:all"; // "all" is the overall result
      return cls;
    };
    // an entrapment hit that equals a target peptide up to I/L cannot be told apart from it
    auto isEntrapment = [&](const PeptideHit& hit, Kind kind)
    {
      return (kind == Kind::ENTRAPMENT) &&
        (target_peptides == nullptr || !target_peptides->contains(toIL(hit.getSequence().toUnmodifiedString())));
    };
    auto accepted = [&](double q) { return std::isfinite(q) && q >= 0.0 && q <= q_threshold; }; // also rejects NaN

    if (level == Level::PSM)
    {
      for (const auto& pep_id : pep_ids)
      {
        if (pep_id.getHits().empty()) continue;
        const PeptideHit& hit = pep_id.getHits().front();
        const Kind kind = classify(hit, entrapment_prefix);
        if (kind == Kind::DECOY) continue;
        if (!accepted(qValue(pep_id, hit))) continue;
        const bool is_entrapment = isEntrapment(hit, kind);
        for (const std::string& key : {std::string("all"), class_meta_value.empty() ? std::string() : classOf(hit)})
        {
          if (key.empty()) continue;
          Result& r = result(key);
          ++r.accepted;
          if (is_entrapment) ++r.entrapment;
        }
      }
      return results;
    }

    // peptide level: one entry per unmodified sequence (I = L), ranked and classed by its best PSM
    struct Peptide
    {
      double q = std::numeric_limits<double>::infinity();
      double score = -std::numeric_limits<double>::infinity();
      bool has_score = false;
      bool entrapment = true; // only if all its PSMs are entrapment
      std::string cls;
    };
    std::unordered_map<std::string, Peptide> peptides;
    for (const auto& pep_id : pep_ids)
    {
      if (pep_id.getHits().empty()) continue;
      const PeptideHit& hit = pep_id.getHits().front();
      const Kind kind = classify(hit, entrapment_prefix);
      if (kind == Kind::DECOY) continue;
      const double q = hit.metaValueExists(Constants::UserParam::PEPTIDE_Q_VALUE) ?
        static_cast<double>(hit.getMetaValue(Constants::UserParam::PEPTIDE_Q_VALUE)) : qValue(pep_id, hit);
      if (!std::isfinite(q) || q < 0.0) continue;
      const bool has_score = !score_meta_value.empty() && hit.metaValueExists(score_meta_value);
      const double score = has_score ? static_cast<double>(hit.getMetaValue(score_meta_value)) : -std::numeric_limits<double>::infinity();

      Peptide& p = peptides[toIL(hit.getSequence().toUnmodifiedString())];
      if (p.cls.empty() || q < p.q || (q == p.q && score > p.score)) // new, or a better PSM
      {
        p.q = q;
        p.score = score;
        p.has_score = has_score;
        p.cls = class_meta_value.empty() ? std::string("all") : classOf(hit);
      }
      p.entrapment = p.entrapment && isEntrapment(hit, kind);
    }

    for (const auto& [sequence, p] : peptides)
    {
      if (!accepted(p.q)) continue;
      // paired counts (Wen et al.): where does the paired target rank?
      bool target_not_accepted = false, above_target = false;
      if (paired && p.entrapment)
      {
        const auto pair_it = pairs->find(sequence);
        if (pair_it == pairs->end())
        {
          above_target = true; // no pair known: count it like the worst case
        }
        else
        {
          const auto target_it = peptides.find(pair_it->second);
          if (target_it == peptides.end() || target_it->second.entrapment || !accepted(target_it->second.q))
          {
            target_not_accepted = true;
          }
          else
          {
            const Peptide& t = target_it->second;
            const bool target_above = (t.q < p.q) || (t.q == p.q && t.has_score && p.has_score && t.score > p.score);
            above_target = !target_above;
          }
        }
      }
      for (const std::string& key : {std::string("all"), class_meta_value.empty() ? std::string() : p.cls})
      {
        if (key.empty()) continue;
        Result& r = result(key);
        ++r.accepted;
        if (p.entrapment) ++r.entrapment;
        if (target_not_accepted) ++r.entrapment_target_not_accepted;
        if (above_target) ++r.entrapment_above_target;
      }
    }
    return results;
  }
} // namespace OpenMS
