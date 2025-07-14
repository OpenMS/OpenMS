// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/PeptideIdentificationList.h>

namespace OpenMS
{
  const String& PeptideIdentificationList::getScoreType() const
  {
    return score_type_;
  }

  void PeptideIdentificationList::setScoreType(const String& type)
  {
    score_type_ = type;
  }

  bool PeptideIdentificationList::isHigherScoreBetter() const
  {
    return higher_score_better_;
  }

  void PeptideIdentificationList::setHigherScoreBetter(bool value)
  {
    higher_score_better_ = value;
  }

  String PeptideIdentificationList::getEffectiveScoreType() const
  {
    if (!score_type_.empty())
    {
      return score_type_;
    }
    if (!empty())
    {
      return front().getScoreType();
    }
    return "";
  }

  bool PeptideIdentificationList::getEffectiveHigherScoreBetter() const
  {
    if (!score_type_.empty())
    {
      return higher_score_better_;
    }
    if (!empty())
    {
      return front().isHigherScoreBetter();
    }
    return false;
  }
}