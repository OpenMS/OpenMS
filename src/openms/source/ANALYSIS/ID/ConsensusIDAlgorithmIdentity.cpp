// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Sven Nahnsen, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmIdentity.h>
#include <OpenMS/CONCEPT/LogStream.h>

using namespace std;

namespace OpenMS
{
  ConsensusIDAlgorithmIdentity::ConsensusIDAlgorithmIdentity()
  {
    setName("ConsensusIDAlgorithmIdentity"); // DefaultParamHandler
  }


  void ConsensusIDAlgorithmIdentity::preprocess_(
    vector<PeptideIdentification>& ids)
  {
    // check score types and orientations:
    bool higher_better = ids[0].isHigherScoreBetter();
    set<String> score_types;

    for (vector<PeptideIdentification>::iterator pep_it = ids.begin();
         pep_it != ids.end(); ++pep_it)
    {
      if (pep_it->isHigherScoreBetter() != higher_better)
      {
        // scores with different orientations definitely aren't comparable:
        String hi_lo = higher_better ? "higher/lower" : "lower/higher";
        String msg = "Score types '" + ids[0].getScoreType() + "' and '" +
          pep_it->getScoreType() + "' have different orientations (" + hi_lo +
          " is better) and cannot be compared meaningfully.";
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      msg, higher_better ? "false" : "true");
      }
      score_types.insert(pep_it->getScoreType());
    }

    if (score_types.size() > 1)
    {
      String types;
      types.concatenate(score_types.begin(), score_types.end(), "'/'");
      OPENMS_LOG_WARN << "Warning: Different score types for peptide hits found ('"
               << types << "'). If the scores are not comparable, "
               << "results will be meaningless." << endl;
    }
  }


  void ConsensusIDAlgorithmIdentity::apply_(vector<PeptideIdentification>& ids,
                                            const map<String, String>& se_info,
                                            SequenceGrouping& results)
  {
    preprocess_(ids);

    // group peptide hits by sequence:
    for (vector<PeptideIdentification>::iterator pep_it = ids.begin();
         pep_it != ids.end(); ++pep_it)
    {
      String score_type = pep_it->getScoreType();
      auto se = se_info.find(pep_it->getIdentifier());
      if (se != se_info.end())
      {
        score_type = se->second + "_" + score_type;
      }

      for (vector<PeptideHit>::iterator hit_it = pep_it->getHits().begin();
           hit_it != pep_it->getHits().end(); ++hit_it)
      {
        const AASequence& seq = hit_it->getSequence();
        auto pos = results.find(seq);
        if (pos == results.end()) // new sequence
        {
          auto ev = hit_it->getPeptideEvidences();
          results[seq] = HitInfo{
              hit_it->getCharge(),
              {hit_it->getScore()},
              {score_type},
              hit_it->getMetaValue("target_decoy").toString(),
              {std::make_move_iterator(ev.begin()), std::make_move_iterator(ev.end())},
              0.,
              0.
          };
        }
        else // previously seen sequence
        {
          compareChargeStates_(pos->second.charge, hit_it->getCharge(),
                               pos->first);
          pos->second.scores.emplace_back(hit_it->getScore());
          pos->second.types.emplace_back(score_type);
          for (const auto& ev : hit_it->getPeptideEvidences())
          {
            pos->second.evidence.emplace(ev);
          }
        }
      }
    }

    // calculate score and support, and update results with them:
    bool higher_better = ids[0].isHigherScoreBetter();
    Size n_other_ids = (count_empty_ ? number_of_runs_ : ids.size()) - 1;
    for (SequenceGrouping::iterator res_it = results.begin(); 
         res_it != results.end(); ++res_it)
    {
      double score = getAggregateScore_(res_it->second.scores, higher_better);
      // if 'count_empty' is false, 'n_other_ids' may be zero, in which case
      // we define the support to be one to avoid a NaN:
      double support = 1.0;
      if (n_other_ids > 0) // the normal case
      {
        support = (res_it->second.scores.size() - 1.0) / n_other_ids;
      }
      res_it->second.final_score = score;
      res_it->second.support = support;
    }
  }

} // namespace OpenMS
