// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
        // In groups you usually want to check if at least one member is a real target
        if (decoy_accs.find(acc) == decoy_accs.end())
        {
          target = true;
          break;
        }
      }
      scores_labels.emplace_back(score, target);
    }
  }


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
} // namespace std
