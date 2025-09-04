// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <unordered_map>

using namespace std;
namespace OpenMS
{

  IDScoreSwitcherAlgorithm::IDScoreSwitcherAlgorithm() :
      IDScoreSwitcherAlgorithm::DefaultParamHandler("IDScoreSwitcherAlgorithm")
  {
    defaults_.setValue("new_score", "", "Name of the meta value to use as the new score");
    defaults_.setValue("new_score_orientation", "", "Orientation of the new score (are higher or lower values better?)");
    defaults_.setValidStrings("new_score_orientation", {"lower_better","higher_better"});
    defaults_.setValue("new_score_type", "", "Name to use as the type of the new score (default: same as 'new_score')");
    defaults_.setValue("old_score", "", "Name to use for the meta value storing the old score (default: old score type)");
    defaults_.setValue("proteins", "false", "Apply to protein scores instead of PSM scores");
    defaults_.setValidStrings("proteins", {"true","false"});
    defaultsToParam_();
    updateMembers_();
  }

  void IDScoreSwitcherAlgorithm::updateMembers_()
  {
    new_score_ = param_.getValue("new_score").toString();
    new_score_type_ = param_.getValue("new_score_type").toString();
    old_score_ = param_.getValue("old_score").toString();
    higher_better_ = (param_.getValue("new_score_orientation").toString() ==
                      "higher_better");

    if (new_score_type_.empty()) new_score_type_ = new_score_;
  }

  std::vector<String> IDScoreSwitcherAlgorithm::getScoreNames()
  {
    std::vector<String> names;
    for (auto i : type_to_str_)
    {
      const std::set<String>& n = i.second;
      for (auto j : n)
      {
        names.push_back(j);
      }
    }
    return names;
  }

  IDScoreSwitcherAlgorithm::ScoreSearchResult IDScoreSwitcherAlgorithm::findScoreType(const PeptideIdentification& id, ScoreType score_type)
  {
    ScoreSearchResult result;
    
    // First check if main score is already of the requested score type using existing infrastructure
    const String& main_score_type = id.getScoreType();
    result.is_main_score_type = isScoreType(main_score_type, score_type);
    
    // If main score is not of the requested type, look for it in meta values
    if (!result.is_main_score_type && !id.getHits().empty())
    {
      const auto& first_hit = id.getHits()[0];
      const std::set<String>& score_types = type_to_str_.at(score_type);
      
      // Search for scores of the requested type in meta values using the existing score type collection
      for (const String& score_name : score_types)
      {
        if (first_hit.metaValueExists(score_name))
        {
          result.meta_value_name = score_name;
          break;
        }
        // Also check for "_score" suffix variant
        String score_name_with_suffix = score_name + "_score";
        if (first_hit.metaValueExists(score_name_with_suffix))
        {
          result.meta_value_name = score_name_with_suffix;
          break;
        }
      }
    }
    
    return result;
  }


} // namespace OpenMS
