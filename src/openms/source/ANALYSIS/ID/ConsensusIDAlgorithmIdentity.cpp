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

#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmIdentity.h>
#include <OpenMS/CONCEPT/LogStream.h>

using namespace std;

#define DEBUG_ID_CONSENSUS
#undef  DEBUG_ID_CONSENSUS

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
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                      msg, higher_better ? "false" : "true");
      }
      score_types.insert(pep_it->getScoreType());
    }

    if (score_types.size() > 1)
    {
      String types;
      types.concatenate(score_types.begin(), score_types.end(), "'/'");
      LOG_WARN << "Warning: Different score types for peptide hits found ('"
               << types << "'). If the scores are not comparable, "
               << "results will be meaningless." << endl;
    }
  }


  void ConsensusIDAlgorithmIdentity::filter_(SequenceGrouping& grouping,
                                             Size n_ids)
  {
    SequenceGrouping filtered;
    for (SequenceGrouping::iterator it = grouping.begin(); it != grouping.end();
         ++it)
    {
      Size n_scores = it->second.second.size();
      if (count_empty_) n_ids = number_of_runs_;
      if ((n_scores - 1) / float(n_ids - 1) >= min_support_)
      {
        filtered.insert(*it);
      }
    }
    grouping = filtered;
  }


  void ConsensusIDAlgorithmIdentity::apply_(vector<PeptideIdentification>& ids)
  {
    preprocess_(ids);

    // group peptide hits by sequence:
    SequenceGrouping grouping;
    for (vector<PeptideIdentification>::iterator pep_it = ids.begin();
         pep_it != ids.end(); ++pep_it)
    {
      for (vector<PeptideHit>::iterator hit_it = pep_it->getHits().begin();
           hit_it != pep_it->getHits().end(); ++hit_it)
      {
        const AASequence& seq = hit_it->getSequence();
        SequenceGrouping::iterator pos = grouping.find(seq);
        if (pos == grouping.end()) // new sequence
        {
          grouping[seq] = make_pair(hit_it->getCharge(), 
                                    vector<double>(1, hit_it->getScore()));
        }
        else // previously seen sequence
        { 
          compareChargeStates_(pos->second.first, hit_it->getCharge(),
                               pos->first);
          pos->second.second.push_back(hit_it->getScore());
        }
      }
    }

    // filter by number of identifications:
    if (min_support_ > 0.0) filter_(grouping, ids.size());

    String score_type = ids[0].getScoreType();
    bool higher_better = ids[0].isHigherScoreBetter();
    ids.clear();
    ids.resize(1);
    ids[0].setScoreType(score_type);
    ids[0].setHigherScoreBetter(higher_better);
    for (SequenceGrouping::iterator group_it = grouping.begin(); 
         group_it != grouping.end(); ++group_it)
    {
      PeptideHit hit;
      hit.setSequence(group_it->first);
      hit.setCharge(group_it->second.first);
      double score = getAggregateScore_(group_it->second.second, higher_better);
      hit.setScore(score);
      ids[0].insertHit(hit);
#ifdef DEBUG_ID_CONSENSUS
      LOG_DEBUG << " - Output hit: " << hit.getSequence() << " "
                << hit.getScore() << endl;
#endif
    }
  }

} // namespace OpenMS
