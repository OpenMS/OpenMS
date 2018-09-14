// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_ID_SCOREDPROCESSINGRESULT_H
#define OPENMS_METADATA_ID_SCOREDPROCESSINGRESULT_H

#include <OpenMS/METADATA/ID/DataProcessingStep.h>
#include <OpenMS/METADATA/ID/ScoreType.h>

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    /// Base class for ID data with scores and processing steps (and meta info)
    struct ScoredProcessingResult: public MetaInfoInterface
    {
      ScoreList scores;

      // @TODO: use a "boost::multi_index_container" here for efficient look-up?
      std::vector<ProcessingStepRef> processing_step_refs;

      /// Merge in data from another objects
      ScoredProcessingResult& operator+=(const ScoredProcessingResult& other)
      {
        // merge processing steps:
        for (auto step_ref : other.processing_step_refs)
        {
          if (std::find(processing_step_refs.begin(),
                        processing_step_refs.end(), step_ref) ==
              processing_step_refs.end())
          {
            processing_step_refs.push_back(step_ref);
          }
        }
        // merge scores:
        for (auto score_pair : other.scores)
        {
          // @TODO: should we overwrite scores?
          if (std::find(scores.begin(), scores.end(), score_pair) ==
              scores.end())
          {
            scores.push_back(score_pair);
          }
        }
        // merge meta info:
        std::vector<UInt> keys;
        other.getKeys(keys);
        for (const UInt key : keys)
        {
          // @TODO: should we overwrite meta values?
          if (!metaValueExists(key))
          {
            setMetaValue(key, other.getMetaValue(key));
          }
        }

        return *this;
      }

      std::pair<double, bool> getScore(ScoreTypeRef score_ref) const
      {
        // give priority to "later" scores in the list:
        for (ScoreList::const_reverse_iterator it = scores.rbegin();
             it != scores.rend(); ++it)
        {
          if (it->first == score_ref) return std::make_pair(it->second, true);
        }
        return std::make_pair(std::numeric_limits<double>::quiet_NaN(), false);
      }

    protected:
      explicit ScoredProcessingResult(
        const ScoreList& scores = ScoreList(),
        const std::vector<ProcessingStepRef>& processing_step_refs =
        std::vector<ProcessingStepRef>()):
        scores(scores), processing_step_refs(processing_step_refs)
      {
      }

      ScoredProcessingResult(const ScoredProcessingResult& other) = default;
    };

  }
}

#endif
