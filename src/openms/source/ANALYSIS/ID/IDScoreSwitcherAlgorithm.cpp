// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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


} // namespace OpenMS
