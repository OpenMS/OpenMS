// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/ID/ScoredProcessingResult.h>

#include <limits>
#include <ranges>

namespace OpenMS::IdentificationDataInternal
{
  void ScoredProcessingResult::addProcessingStep(const AppliedProcessingStep& step)
  {
    auto step_pos =
      steps_and_scores.get<1>().find(step.processing_step_opt);
    if (step_pos == steps_and_scores.get<1>().end()) // new step
    {
      steps_and_scores.push_back(step);
    }
    else // existing step - add or update scores
    {
      steps_and_scores.get<1>().modify(
        step_pos, [&](AppliedProcessingStep& old_step)
        {
          for (const auto& pair : step.scores)
          {
            old_step.scores[pair.first] = pair.second;
          }
        });
    }
  }

  void ScoredProcessingResult::addProcessingStep(ProcessingStepRef step_ref,
                                                 const std::map<ScoreTypeRef, double>& scores)
  {
    AppliedProcessingStep applied(step_ref, scores);
    addProcessingStep(applied);
  }

  void ScoredProcessingResult::addScore(ScoreTypeRef score_type, double score,
                                        const std::optional<ProcessingStepRef>& processing_step_opt)
  {
    AppliedProcessingStep applied(processing_step_opt);
    applied.scores[score_type] = score;
    addProcessingStep(applied);
  }

  ScoredProcessingResult& ScoredProcessingResult::merge(const ScoredProcessingResult& other)
  {
    // merge applied processing steps and scores:
    for (const auto& step : other.steps_and_scores)
    {
      addProcessingStep(step);
    }
    // merge meta info - existing entries may be overwritten:
    addMetaValues(other);

    return *this;
  }

  std::pair<double, bool> ScoredProcessingResult::getScore(ScoreTypeRef score_ref) const
  {
    std::tuple<double, std::optional<ProcessingStepRef>, bool> result =
      getScoreAndStep(score_ref);
    return std::make_pair(std::get<0>(result), std::get<2>(result));
  }

  std::pair<double, bool> ScoredProcessingResult::getScore(ScoreTypeRef score_ref,
                                                           std::optional<ProcessingStepRef> processing_step_opt) const
  {
    auto step_pos = steps_and_scores.get<1>().find(processing_step_opt);
    if (step_pos != steps_and_scores.get<1>().end())
    {
      auto score_pos = step_pos->scores.find(score_ref);
      if (score_pos != step_pos->scores.end())
      {
        return std::make_pair(score_pos->second, true);
      }
    }
    // not found:
    return std::make_pair(std::numeric_limits<double>::quiet_NaN(), false);
  }

  std::tuple<double, std::optional<ProcessingStepRef>, bool>
  ScoredProcessingResult::getScoreAndStep(ScoreTypeRef score_ref) const
  {
    // give priority to scores from later processing steps:
    for (const auto& step : std::views::reverse(steps_and_scores))
    {
      auto pos = step.scores.find(score_ref);
      if (pos != step.scores.end())
      {
        return std::make_tuple(pos->second, step.processing_step_opt, true);
      }
    }
    // not found:
    return std::make_tuple(std::numeric_limits<double>::quiet_NaN(),
                           std::nullopt, false);
  }

  std::tuple<double, std::optional<ScoreTypeRef>, bool>
  ScoredProcessingResult::getMostRecentScore() const
  {
    // check steps starting with most recent:
    for (const auto& step : std::views::reverse(steps_and_scores))
    {
      auto top_score = step.getScoresInOrder(true);
      if (!top_score.empty())
      {
        return std::make_tuple(top_score[0].second, top_score[0].first,
                               true);
      }
    }
    return std::make_tuple(std::numeric_limits<double>::quiet_NaN(),
                           std::nullopt, false); // no score available
  }

  Size ScoredProcessingResult::getNumberOfScores() const
  {
    Size counter = 0;
    for (const auto& step : steps_and_scores)
    {
      counter += step.scores.size();
    }
    return counter;
  }
} // namespace OpenMS::IdentificationDataInternal
