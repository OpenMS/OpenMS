// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/IDScoreGetterSetter.h>

#include <OpenMS/DATASTRUCTURES/StringView.h>

using namespace std;

namespace OpenMS
{

  /**
   * @ingroup getScoresFunctions
   * @brief For protein groups. Groups are target if at least one protein is target
   * Therefore it needs an unordered set of decoy accessions to evaluate that.
  */
  void IDScoreGetterSetter::getScores_(
      ScoreToTgtDecLabelPairs &scores_labels,
      const std::vector<ProteinIdentification::ProteinGroup> &grps,
      const std::unordered_set<string> &decoy_accs)
  {
    for (const auto &grp : grps)
    {
      double score = grp.probability;
      bool target = false;
      for (const auto &acc : grp.accessions)
      {
        // In groups, you usually want to check if at least one member is a real target
        if (decoy_accs.find(acc) == decoy_accs.end())
        {
          target = true;
          break;
        }
      }
      scores_labels.emplace_back(score, target);
    }
  }

  inline bool isFirstBetterScore_(double first, double second, bool isHigherBetter)
  {
    if (isHigherBetter) return first > second; else return first < second;
  }

  inline void addToPeptideScoreMap_(
    std::unordered_map<String, ScoreToTgtDecLabelPair>& seq_to_score_labels,
    const PeptideIdentification& id)
  {
    bool higher_better = id.isHigherScoreBetter();
    if (id.getHits().empty())
    {
      return;
    }
    const auto& best_hit = id.getHits()[0];
    double score = best_hit.getScore();
    auto [it, found] = seq_to_score_labels.try_emplace(
      best_hit.getSequence().toUnmodifiedString(),
      score,
      (best_hit.getMetaValue("target_decoy") != DataValue::EMPTY) &&
        (best_hit.getMetaValue("target_decoy").toString().hasPrefix("target")));

    if (found && isFirstBetterScore_(score, it->second.first, higher_better))
    {
      it->second.first = score;
    }
  }

  void IDScoreGetterSetter::fillPeptideScoreMap_(
    std::unordered_map<String, ScoreToTgtDecLabelPair>& seq_to_score_labels,
    const vector<PeptideIdentification>& ids)
  {
    for (auto const & id : ids)
    {
      addToPeptideScoreMap_(seq_to_score_labels, id);
    }
  }

  void IDScoreGetterSetter::fillPeptideScoreMap_(
    std::unordered_map<String, ScoreToTgtDecLabelPair>& seq_to_score_labels,
    ConsensusMap const& map,
    bool include_unassigned = true)
  {
    map.applyFunctionOnPeptideIDs(
      [&seq_to_score_labels](const PeptideIdentification& id){addToPeptideScoreMap_(seq_to_score_labels, id);},
      include_unassigned);
  }

  /**
   * @ingroup getScoresFunctions
   * @brief For protein groups. Groups are target if at least one protein is target
   * Decoy accessions are determined by the decoy substring.
   * Uses the "picked" algorithm. As soon as there was one member which was picked as target over a decoy, the group is counted as target. Otherwise as decoy.
   */
  void IDScoreGetterSetter::getPickedProteinGroupScores_(
      const std::unordered_map<String, ScoreToTgtDecLabelPair>& picked_scores,
      ScoreToTgtDecLabelPairs& scores_labels,
      const vector<ProteinIdentification::ProteinGroup>& grps,
      const String& decoy_string,
      bool decoy_prefix)
  {
    for (const auto& grp : grps)
    {
      bool decoy_picked = false;
      for (const auto& acc : grp.accessions)
      {
        auto [isDecoy, tgt_accession] = removeDecoyStringIfPresent_(acc, decoy_string, decoy_prefix);
        const double tgt_proportion = picked_scores.at(tgt_accession).second;
        // The problem is here, that in theory, the matching pair can be in different groups
        //  therefore a group cannot really be matched one-to-one. So we say:
        //  If at least one (single) target was picked over a decoy in the group, the group
        //  is a target. This could in theory mean, that two targets in different groups
        //  are counted as +1, while a group containing both decoy-partners is counted as
        //  a single 0
        if (!isDecoy && tgt_proportion > 0.) // target was picked on single protein level
        {
          scores_labels.emplace_back(grp.probability, 1.0);
          break;
        }
        else if (isDecoy && tgt_proportion == 0)
        {
          decoy_picked = true;
        }
      }
      // if for none of the proteins the target version was picked, add as decoy
      if (decoy_picked) scores_labels.emplace_back(grp.probability, 0.0);
      // TODO I think we need an unordered_set to check which proteins were picked already
      //  and add skip groups where every protein was picked already.
    }
  }

  pair<bool,String> IDScoreGetterSetter::removeDecoyStringIfPresent_(const String& acc, const String& decoy_string, bool decoy_prefix)
  {
    if (decoy_prefix && acc.hasPrefix(decoy_string))
    {
      return {true ,acc.suffix(acc.size() - decoy_string.size())};
    }
    else if (acc.hasSuffix(decoy_string))
    {
      return {true, acc.prefix(acc.size() - decoy_string.size())};
    }
    else
    {
      return {false, acc};
    }
  }

  void IDScoreGetterSetter::getPickedProteinScores_(
      std::unordered_map<String, ScoreToTgtDecLabelPair>& picked_scores,
      const ProteinIdentification& id,
      const String& decoy_string,
      bool decoy_prefix)
  {
    for (const auto& hit : id.getHits())
    {
      checkTDAnnotation_(hit);
      StringView tgt_accession(hit.getAccession());
      bool target = getTDLabel_(hit);
      if (!target)
      {
        if (decoy_prefix) //TODO double-check hasSuffix/Prefix? Ignore TD Metavalue?
        {
          tgt_accession = tgt_accession.substr(decoy_string.size(),-1);
        }
        else
        {
          tgt_accession = tgt_accession.substr(0,tgt_accession.size()-decoy_string.size());
        }
      }
      auto[it, inserted] = picked_scores.try_emplace(tgt_accession.getString(), hit.getScore(), target);
      if (!inserted)
      {
        if ((id.isHigherScoreBetter() && (hit.getScore() > it->second.first)) ||
        (!id.isHigherScoreBetter() && (hit.getScore() < it->second.first)))
        {
          it->second = {hit.getScore(), target};
        }
        else if (hit.getScore() == it->second.first)
        {
          it->second = {hit.getScore(), true}; //prefer targets. Alternative: put 0.5
        }
      }
    }
  }

  /*static void getPickedProteinScores_(
      ScoreToTgtDecLabelPairs& scores_labels,
      const std::vector<ProteinIdentification::ProteinGroup> &grps,
      const std::unordered_set<std::string> &decoy_accs,
      const String& decoy_string,
      bool prefix)
  {

    //TODO potential algorithm: Create a winner set based on single protein scores
    // Iff a group contains at least one winner, add the group with its group score to the
    // vector (for input to group FDR).
    // Otherwise I feel like groups would block/steal too many singles/small groups
    // On the other hand, with aggregational inference groups and singles will have the same scores anyway
    std::unordered_map<String, std::pair<double, double>> picked_scores;
    for (const auto& grp : grps)
    {
      StringView tgt_accession(grp.accessions);
      bool target = getTDLabel_(hit);
      if (!target)
      {
        if (decoy_prefix)
        {
          tgt_accession = tgt_accession.substr(decoy_string.size(),-1);
        }
        else
        {
          tgt_accession = tgt_accession.substr(0,tgt_accession.size()-decoy_string.size());
        }
      }
      auto[it, inserted] = picked_scores.try_emplace(tgt_accession.getString(), hit.getScore(), target);
      if (!inserted)
      {
        if ((id.isHigherScoreBetter() && (hit.getScore() > it->second.first)) ||
        (!id.isHigherScoreBetter() && (hit.getScore() < it->second.first)))
        {
          it->second = {hit.getScore(), target};
        }
        else if (hit.getScore() == it->second.first)
        {
          it->second = {hit.getScore(), true}; //prefer targets
        }
      }
    }

    scores_labels.reserve(picked_scores.size());
    for(auto& kv : picked_scores)
    {
      scores_labels.emplace_back(std::move(kv.second));
    }
  }*/

  /** @ingroup setScoresFunctions
  * @brief For protein groups. Unaffected by keep_decoy_proteins. Always keeps all for now @todo.
  * score_type and higher_better unused since ProteinGroups do not carry that information.
  * You have to assume that groups will always have the same scores as the ProteinHits
  */
  void IDScoreGetterSetter::setScores_(const map<double, double> &scores_to_FDR,
                                      vector <ProteinIdentification::ProteinGroup> &grps,
                                      const string & /*score_type*/,
                                      bool /*higher_better*/)
  {
    for (auto &grp : grps)
    {
      grp.probability = (scores_to_FDR.lower_bound(grp.probability)->second);
    }
  }
  void IDScoreGetterSetter::setPeptideScoresFromMap_(std::unordered_map<String, ScoreToTgtDecLabelPair> const& seq_to_fdr,
                                                     vector<PeptideIdentification>& ids,
                                                     std::string const& score_type,
                                                     bool keep_decoys)
  {
    for (auto& id : ids)
    {
      if (id.getHits().empty())
      {
        continue;
      }
      auto& best_hit = id.getHits()[0];
      if (!keep_decoys && (best_hit.getMetaValue("target_decoy") == DataValue::EMPTY || best_hit.getMetaValue("target_decoy") == "decoy"))
      {
        id.setHits({});
        continue;
      }
      const auto seq = best_hit.getSequence().toUnmodifiedString();
      auto it = seq_to_fdr.find(seq);
      const auto& old_score_type = id.getScoreType();
      if (it != seq_to_fdr.end())
      {
        best_hit.setMetaValue(old_score_type, best_hit.getScore());
        best_hit.setScore(it->second.first);
        id.setScoreType(score_type);
      }
      else
      {
        OPENMS_LOG_ERROR << "Error: No FDR found for " + seq + "." << std::endl;
        continue;
      }
    }
  }

  void IDScoreGetterSetter::setPeptideScoresFromMap_(std::unordered_map<String, ScoreToTgtDecLabelPair> const& seq_to_fdr,
                                                     ConsensusMap& map,
                                                     std::string const& score_type,
                                                     bool keep_decoys,
                                                     bool include_unassigned)
  {
    for (auto& f : map)
    {
      setPeptideScoresFromMap_(seq_to_fdr, f.getPeptideIdentifications(), score_type, keep_decoys);
    }
    if (include_unassigned)
    {
      setPeptideScoresFromMap_(seq_to_fdr, map.getUnassignedPeptideIdentifications(), score_type, keep_decoys);
    }
  }

} // namespace std
