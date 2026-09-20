// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/ID/IDDataContainer.h>

#include <OpenMS/METADATA/ID/ProcessingStep.h>
#include <OpenMS/METADATA/ID/ScoreType.h>


#include <optional>

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

        @param[in] primary_only Only return the primary score (ignoring any others)?
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
          if (!scores_done.contains(pair.first))
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
    using AppliedProcessingSteps = IDDataContainer<AppliedProcessingStep, std::optional<ProcessingStepRef>, std::optional<ProcessingStepRef>>;
    extern template class OPENMS_DLLAPI IDDataContainer<AppliedProcessingStep, std::optional<ProcessingStepRef>, std::optional<ProcessingStepRef>>;
  }
}
