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

#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmSimilarity.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <cmath>
#include <numeric> // for "accumulate"

using namespace std;

#define DEBUG_ID_CONSENSUS
#undef  DEBUG_ID_CONSENSUS

namespace OpenMS
{
  ConsensusIDAlgorithmSimilarity::ConsensusIDAlgorithmSimilarity()
  {
    setName("ConsensusIDAlgorithmSimilarity"); // DefaultParamHandler
  }


  void ConsensusIDAlgorithmSimilarity::apply_(
    vector<PeptideIdentification>& ids)
  {
    for (vector<PeptideIdentification>::iterator id = ids.begin();
         id != ids.end(); ++id)
    {
      if (id->getScoreType() != "Posterior Error Probability")
      {
        String msg = "Score type must be 'Posterior Error Probablity'";
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                      msg, id->getScoreType());
      }
    }

    // mapping: peptide sequence -> (charge, [consensus score, similarity])
    SequenceGrouping results;
    
    for (vector<PeptideIdentification>::iterator id1 = ids.begin();
         id1 != ids.end(); ++id1)
    {
      for (vector<PeptideHit>::iterator hit1 = id1->getHits().begin();
           hit1 != id1->getHits().end(); ++hit1)
      {
        // have we scored this sequence already? if yes, skip:
        SequenceGrouping::iterator pos = results.find(hit1->getSequence());
        if (pos != results.end())
        { 
          compareChargeStates_(pos->second.first, hit1->getCharge(),
                               pos->first);
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

        vector<double> scores(2, score);
        // normalize similarity score to range 0-1:
        if (best_matches.empty()) // only one ID run -> similarity is undefined
        {
          scores[1] = -1.0;
        }
        else
        {
          scores[1] = (sum_sim - 1.0) / best_matches.size();
        }

        results[hit1->getSequence()] = make_pair(hit1->getCharge(), scores);
      }
    }

    ids.clear();
    ids.resize(1);
    ids[0].setScoreType("Posterior Error Probability");
    ids[0].setHigherScoreBetter(false);
    for (SequenceGrouping::iterator res_it = results.begin(); 
         res_it != results.end(); ++res_it)
    {
      PeptideHit hit;
      hit.setSequence(res_it->first);
      hit.setScore(res_it->second.second[0]);
      hit.setCharge(res_it->second.first);
      hit.setMetaValue("consensus_similarity", res_it->second.second[1]);
      ids[0].insertHit(hit);
#ifdef DEBUG_ID_CONSENSUS
      LOG_DEBUG << " - Output hit: " << hit.getSequence() << " "
                << hit.getScore() << endl;
#endif
    }
  }

} // namespace OpenMS
