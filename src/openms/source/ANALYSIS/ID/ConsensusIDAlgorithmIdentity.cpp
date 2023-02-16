// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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

    for (PeptideIdentification& pep : ids)
    {
      if (pep.isHigherScoreBetter() != higher_better)
      {
        // scores with different orientations definitely aren't comparable:
        String hi_lo = higher_better ? "higher/lower" : "lower/higher";
        String msg = "Score types '" + ids[0].getScoreType() + "' and '" +
          pep.getScoreType() + "' have different orientations (" + hi_lo +
          " is better) and cannot be compared meaningfully.";
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      msg, higher_better ? "false" : "true");
      }
      score_types.insert(pep.getScoreType());
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
    for (PeptideIdentification& pep : ids)
    {
      String score_type = pep.getScoreType();
      auto se = se_info.find(pep.getIdentifier());
      if (se != se_info.end())
      {
        score_type = se->second + "_" + score_type;
      }

      for (PeptideHit& hit : pep.getHits())
      {
        const AASequence& seq = hit.getSequence();
        auto pos = results.find(seq);
        if (pos == results.end()) // new sequence
        {
          auto ev = hit.getPeptideEvidences();
          results[seq] = HitInfo{
              hit.getCharge(),
              {hit.getScore()},
              {score_type},
              hit.getMetaValue("target_decoy").toString(),
              {std::make_move_iterator(ev.begin()), std::make_move_iterator(ev.end())},
              0.,
              0.
          };
        }
        else // previously seen sequence
        {
          compareChargeStates_(pos->second.charge, hit.getCharge(),
                               pos->first);
          pos->second.scores.emplace_back(hit.getScore());
          pos->second.types.emplace_back(score_type);
          for (const auto& ev : hit.getPeptideEvidences())
          {
            pos->second.evidence.emplace(ev);
          }
        }
      }
    }

    // calculate score and support, and update results with them:
    bool higher_better = ids[0].isHigherScoreBetter();
    Size n_other_ids = (count_empty_ ? number_of_runs_ : ids.size()) - 1;
    for (auto& res : results)
    {
      double score = getAggregateScore_(res.second.scores, higher_better);
      // if 'count_empty' is false, 'n_other_ids' may be zero, in which case
      // we define the support to be one to avoid a NaN:
      double support = 1.0;
      if (n_other_ids > 0) // the normal case
      {
        support = (res.second.scores.size() - 1.0) / n_other_ids;
      }
      res.second.final_score = score;
      res.second.support = support;
    }
  }

} // namespace OpenMS
