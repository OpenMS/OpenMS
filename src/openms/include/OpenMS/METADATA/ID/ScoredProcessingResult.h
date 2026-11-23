// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/ID/AppliedProcessingStep.h>
#include <map>
#include <optional>
#include <tuple>
#include <utility>

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    /// Base class for ID data with scores and processing steps (and meta info)
    struct OPENMS_DLLAPI ScoredProcessingResult: public MetaInfoInterface
    {
      AppliedProcessingSteps steps_and_scores;

      /// Return the applied processing steps (incl. scores) as a set ordered by processing step reference (option)
      AppliedProcessingSteps::OrderedView getStepsAndScoresByStep()
      {
        return steps_and_scores.get<1>();
      }

      /// Return the applied processing steps (incl. scores) as a set ordered by processing step reference (option) - const variant
      AppliedProcessingSteps::ConstOrderedView getStepsAndScoresByStep() const
      {
        return steps_and_scores.get<1>();
      }

      /*!
        @brief Add an applied processing step

        If the step already exists, scores are merged (existing ones updated).
      */
      void addProcessingStep(const AppliedProcessingStep& step);

      /// Add a processing step (and associated scores, if any)
      void addProcessingStep(ProcessingStepRef step_ref,
                             const std::map<ScoreTypeRef, double>& scores =
                             std::map<ScoreTypeRef, double>());

      /// Add a score (possibly connected to a processing step)
      void addScore(ScoreTypeRef score_type, double score,
                    const std::optional<ProcessingStepRef>&
                    processing_step_opt = std::nullopt);

      /// Merge in data from another object
      ScoredProcessingResult& merge(const ScoredProcessingResult& other);

      /*!
        @brief Look up a score by score type

        All processing steps are considered, in "most recent first" order.

        @return A pair: score (or NaN), success indicator
      */
      std::pair<double, bool> getScore(ScoreTypeRef score_ref) const;

      /*!
        @brief Look up a score by score type and processing step

        @return A pair: score (or NaN), success indicator
      */
      std::pair<double, bool> getScore(ScoreTypeRef score_ref,
                                       std::optional<ProcessingStepRef>
                                       processing_step_opt) const;

      /*!
        @brief Look up a score and associated processing step by score type

        All processing steps are considered, in "most recent first" order.

        @return A tuple: score (or NaN), processing step reference (option), success indicator
      */
      std::tuple<double, std::optional<ProcessingStepRef>, bool>
      getScoreAndStep(ScoreTypeRef score_ref) const;

      /*!
        @brief Get the (primary) score from the most recent processing step

        This will fail if no scores have been assigned.

        @return A tuple: score (or NaN), score type reference (option), success indicator
        */
      std::tuple<double, std::optional<ScoreTypeRef>, bool>
      getMostRecentScore() const;

      /// Return the number of scores associated with this result
      Size getNumberOfScores() const;

    protected:
      /// Constructor
      explicit ScoredProcessingResult(
        const AppliedProcessingSteps& steps_and_scores =
        AppliedProcessingSteps()):
        steps_and_scores(steps_and_scores)
      {
      }

      /// Copy c'tor
      ScoredProcessingResult(const ScoredProcessingResult&) = default;
    };

  }
}
