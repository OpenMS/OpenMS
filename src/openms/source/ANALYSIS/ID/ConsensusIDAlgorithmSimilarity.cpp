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

#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmSimilarity.h>
#include <OpenMS/CONCEPT/LogStream.h>

using namespace std;

namespace OpenMS
{
  ConsensusIDAlgorithmSimilarity::ConsensusIDAlgorithmSimilarity()
  {
    setName("ConsensusIDAlgorithmSimilarity"); // DefaultParamHandler
  }


  void ConsensusIDAlgorithmSimilarity::apply_(
    vector<PeptideIdentification>& ids,
    const map<String, String>& se_info,
    SequenceGrouping& results)
  {
    for (vector<PeptideIdentification>::iterator id = ids.begin();
         id != ids.end(); ++id)
    {
      if (id->getScoreType() != "Posterior Error Probability" &&
          id->getScoreType() != "pep" &&
          id->getScoreType() != "MS:1001493")
      {
        String msg = "Score type must be 'Posterior Error Probability'";
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      msg, id->getScoreType());
      }
    }

    for (vector<PeptideIdentification>::iterator id1 = ids.begin();
         id1 != ids.end(); ++id1)
    {
      String score_type = id1->getScoreType();
      auto se = se_info.find(id1->getIdentifier());
      if (se != se_info.end())
      {
        score_type = se->second + "_" + score_type;
      }

      for (vector<PeptideHit>::iterator hit1 = id1->getHits().begin();
           hit1 != id1->getHits().end(); ++hit1)
      {
        // have we scored this sequence already? if yes, skip:
        SequenceGrouping::iterator pos = results.find(hit1->getSequence());
        if (pos != results.end())
        { 
          compareChargeStates_(pos->second.charge, hit1->getCharge(),
                               pos->first);
          pos->second.scores.emplace_back(hit1->getScore());
          pos->second.types.emplace_back(id1->getScoreType());
          for (const auto& ev : hit1->getPeptideEvidences())
          {
            pos->second.evidence.emplace(ev);
          }
          continue;
        }
        
        // similarity scores and PEPs of best matches for all ID runs:
        vector<pair<double, double> > best_matches;
        best_matches.reserve(ids.size() - 1);
        for (vector<PeptideIdentification>::iterator id2 = ids.begin();
             id2 != ids.end(); ++id2)
        {
          if (id1 == id2) continue;
          
          // similarity scores and PEPs of all matches in current ID run
          // (to get the best match, we look for highest similarity, breaking
          // ties by better PEP - so we need to transform PEP so higher scores
          // are better, same as similarity):
          vector<pair<double, double> > current_matches;
          current_matches.reserve(id2->getHits().size());
          for (vector<PeptideHit>::iterator hit2 = id2->getHits().begin();
               hit2 != id2->getHits().end(); ++hit2)
          {
            double sim_score = getSimilarity_(hit1->getSequence(),
                                              hit2->getSequence());
            // use "1 - PEP" so higher scores are better (for "max_element"):
            current_matches.push_back(make_pair(sim_score,
                                                1.0 - hit2->getScore()));
          }
          best_matches.push_back(*max_element(current_matches.begin(),
                                              current_matches.end()));
        }
        double score = hit1->getScore();
        double sum_sim = 1.0; // sum of similarity scores
        for (vector<pair<double, double> >::iterator it = best_matches.begin();
             it != best_matches.end(); ++it)
        {
          score += it->first * (1.0 - it->second); // undo "1 - PEP" transform
          sum_sim += it->first;
        }
        score /= (sum_sim * sum_sim);

        double support = 0.;
        // normalize similarity score to range 0-1:
        Size n_other_ids = (count_empty_ ?
                            number_of_runs_ - 1 : best_matches.size());
        if (n_other_ids == 0) // only one ID run -> similarity is ill-defined
        {
          support = double(!count_empty_); // 0 or 1 depending on parameter
        }
        else
        {
          support = (sum_sim - 1.0) / n_other_ids;
        }

        auto ev = hit1->getPeptideEvidences();
        // don't filter based on "min_score_" yet, so we don't recompute results
        // for the same peptide sequence:
        results[hit1->getSequence()] =
            {
              hit1->getCharge(),
              {hit1->getScore()},
              {score_type},
              hit1->getMetaValue("target_decoy").toString(),
              {std::make_move_iterator(ev.begin()), std::make_move_iterator(ev.end())},
              score,
              support
            };
      }
    }
  }

} // namespace OpenMS
