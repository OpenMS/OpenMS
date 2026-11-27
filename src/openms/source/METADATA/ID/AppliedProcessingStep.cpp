// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/ID/AppliedProcessingStep.h>

#include <set>

namespace OpenMS::IdentificationDataInternal
{
  bool AppliedProcessingStep::operator==(const AppliedProcessingStep& other) const
  {
    return ((processing_step_opt == other.processing_step_opt) &&
            (scores == other.scores));
  }

  std::vector<std::pair<ScoreTypeRef, double>>
  AppliedProcessingStep::getScoresInOrder(bool primary_only) const
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
} // namespace OpenMS::IdentificationDataInternal
