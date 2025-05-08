// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Sven Nahnsen, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmRanks.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <numeric> // for "accumulate"

using namespace std;

namespace OpenMS
{
  ConsensusIDAlgorithmRanks::ConsensusIDAlgorithmRanks()
  {
    setName("ConsensusIDAlgorithmRanks"); // DefaultParamHandler
  }


  void ConsensusIDAlgorithmRanks::preprocess_(
    vector<PeptideIdentification>& ids)
  {
    // The idea here is that each peptide hit (sequence) gets assigned a score
    // from each ID run, based on its rank in the list of search results.
    // The best hit of a run will receive score 0, the second best 1, etc. up to
    // a score of N = considered_hits - 1 for the last hit (if there are that
    // many). A hit that was not observed in an ID run receives a score of
    // N + 1 = considered_hits from that run. In the end the scores for each
    // hit (sequence) are averaged and normalized to the range from 0
    // (exclusive, worst) to 1 (inclusive, best).

    current_number_of_runs_ = ((number_of_runs_ > 0) ? 
                               number_of_runs_ : ids.size());
    current_considered_hits_ = considered_hits_;
    bool set_considered_hits = (considered_hits_ == 0);

    for (vector<PeptideIdentification>::iterator pep_it = ids.begin();
         pep_it != ids.end(); ++pep_it)
    {
      pep_it->assignRanks();
      for (vector<PeptideHit>::iterator hit_it = pep_it->getHits().begin();
           hit_it != pep_it->getHits().end(); ++hit_it)
      {
        // give each hit a score based on the search rank (counting from 0):
        hit_it->setScore(hit_it->getRank() - 1);
      }
      pep_it->setScoreType("ConsensusID_ranks");
      pep_it->setHigherScoreBetter(true); // not true now, but after normalizing

      // if "considered_hits" wasn't set, we find the max. number of hits:
      if (set_considered_hits &&
          (pep_it->getHits().size() > current_considered_hits_))
      {
        current_considered_hits_ = pep_it->getHits().size();
      }
    }
  }


  double ConsensusIDAlgorithmRanks::getAggregateScore_(vector<double>& scores,
                                                       bool /* higher_better */)
  {
    double sum_scores = accumulate(scores.begin(), scores.end(), 0.0);
    // add score contributions equivalent to "not found":
    sum_scores += ((current_number_of_runs_ - scores.size()) * 
                   current_considered_hits_);
    // normalize to range 0-1:
    return 1.0 - sum_scores / (current_considered_hits_ * 
                               current_number_of_runs_);
  }

} // namespace OpenMS
