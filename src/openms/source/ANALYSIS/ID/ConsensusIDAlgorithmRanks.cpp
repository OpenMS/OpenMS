// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Hendrik Weisser $
// $Authors: Sven Nahnsen, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmRanks.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <cmath>
#include <numeric> // for "accumulate"

using namespace std;

#define DEBUG_ID_CONSENSUS
#undef  DEBUG_ID_CONSENSUS

namespace OpenMS
{
  ConsensusIDAlgorithmRanks::ConsensusIDAlgorithmRanks()
  {
    setName("ConsensusIDAlgorithmRanks"); // DefaultParamHandler

    defaults_.setValue("number_of_runs", 0, "The number of runs used as input. If not given, the number of input identifications is used.");
    defaults_.setMinInt("number_of_runs", 0);

    defaultsToParam_();
  }


  void ConsensusIDAlgorithmRanks::updateMembers_()
  {
    ConsensusIDAlgorithmIdentity::updateMembers_();

    number_of_runs_ = param_.getValue("number_of_runs");
  }


  void ConsensusIDAlgorithmRanks::apply_(vector<PeptideIdentification>& ids)
  {
    // The idea here is that each peptide hit (sequence) gets assigned a score
    // from each ID run, based on its rank in the list of search results.
    // The best hit of a run will receive score 0, the second best 1, etc. up to
    // a score of N = considered_hits - 1 for the last hit (if there are that
    // many). A hit that was not observed in an ID run receives a score of
    // N + 1 = considered_hits from that run. In the end the scores for each
    // hit (sequence) are averaged and normalized to the range from 0
    // (exclusive, worst) to 1 (inclusive, best).

    Size number_of_runs = (number_of_runs_ > 0) ? number_of_runs_ : ids.size();
    Size considered_hits = considered_hits_;
    bool set_considered_hits = (considered_hits == 0);

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
      pep_it->setScoreType("rank");

      // if "considered_hits" wasn't set, we find the max. number of hits:
      if (set_considered_hits && (pep_it->getHits().size() > considered_hits))
      {
        considered_hits = pep_it->getHits().size();
      }
    }

    SequenceGrouping grouping;
    groupHits_(ids, grouping);

    ids.clear();
    ids.resize(1);
    // the original score type doesn't matter here, so don't include it:
    ids[0].setScoreType("Consensus_ranks");
    ids[0].setHigherScoreBetter(true);
    for (SequenceGrouping::iterator group_it = grouping.begin(); 
         group_it != grouping.end(); ++group_it)
    {
      PeptideHit hit;
      hit.setSequence(group_it->first);
      double sum_scores = accumulate(group_it->second.second.begin(),
                                     group_it->second.second.end(), 0.0);
      // add score contributions equivalent to "not found":
      sum_scores += ((number_of_runs - group_it->second.second.size()) *
                     considered_hits);
      // normalize to range 0-1:
      hit.setScore(1.0 - sum_scores / (considered_hits * number_of_runs));
      hit.setCharge(group_it->second.first);
      ids[0].insertHit(hit);
#ifdef DEBUG_ID_CONSENSUS
      LOG_DEBUG << " - Output hit: " << hit.getSequence() << " "
                << hit.getScore() << endl;
#endif
    }
  }

} // namespace OpenMS
