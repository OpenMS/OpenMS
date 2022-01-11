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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/ID/ProcessingStep.h>
#include <OpenMS/METADATA/ID/ScoreType.h>

#include <boost/optional.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/sequenced_index.hpp>

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    /*!
      A processing step that was applied to a data item, possibly with associated scores.
    */
    struct AppliedProcessingStep
    {
      /*!
        @brief (Optional) reference to the processing step

        If there are only scores, the processing step may be missing.
       */
      std::optional<ProcessingStepRef> processing_step_opt;

      /// Map of scores and their types
      std::map<ScoreTypeRef, double> scores;

      /// Constructor
      explicit AppliedProcessingStep(
        const std::optional<ProcessingStepRef>& processing_step_opt =
        std::nullopt, const std::map<ScoreTypeRef, double>& scores =
        std::map<ScoreTypeRef, double>()):
        processing_step_opt(processing_step_opt), scores(scores)
      {
      }

      /// Equality operator (needed for multi-index container)
      bool operator==(const AppliedProcessingStep& other) const
      {
        return ((processing_step_opt == other.processing_step_opt) &&
                (scores == other.scores));
      }

      /*!
        @brief Return scores in order of priority (primary first).

        The order is defined in the @p ProcessingSoftware referenced by the processing step (if available).
        Scores not listed there are included at the end of the output.

        @param primary_only Only return the primary score (ignoring any others)?
      */
      std::vector<std::pair<ScoreTypeRef, double>>
      getScoresInOrder(bool primary_only = false) const
      {
        std::vector<std::pair<ScoreTypeRef, double>> result;
        std::set<ScoreTypeRef> scores_done;

        if (processing_step_opt)
        {
          ProcessingSoftwareRef sw_ref = (*processing_step_opt)->software_ref;
          for (ScoreTypeRef score_ref : sw_ref->assigned_scores)
          {
            auto pos = scores.find(score_ref);
            if (pos != scores.end())
            {
              result.push_back(*pos);
              if (primary_only) return result;
              scores_done.insert(score_ref);
            }
          }
        }
        for (const auto& pair: scores)
        {
          if (!scores_done.count(pair.first))
          {
            result.push_back(pair);
            if (primary_only) return result;
          }
        }
        return result;
      }
    };

    // we want to keep track of the processing steps in sequence (order of
    // application), but also ensure there are no duplicate steps:
    typedef boost::multi_index_container<
      AppliedProcessingStep,
      boost::multi_index::indexed_by<
        boost::multi_index::sequenced<>,
        boost::multi_index::ordered_unique<
          boost::multi_index::member<
            AppliedProcessingStep, std::optional<ProcessingStepRef>,
            &AppliedProcessingStep::processing_step_opt>>>
      > AppliedProcessingSteps;

  }
}
