// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
