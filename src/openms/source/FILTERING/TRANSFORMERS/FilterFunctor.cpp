// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

#include <OpenMS/FILTERING/TRANSFORMERS/ComplementFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/GoodDiffFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/IntensityBalanceFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NeutralLossDiffFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/IsotopeDiffFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/TICFilter.h>

namespace OpenMS
{
  FilterFunctor::FilterFunctor() :
    DefaultParamHandler("FilterFunctor")
  {
  }

  FilterFunctor::FilterFunctor(const FilterFunctor & source) = default;

  FilterFunctor & FilterFunctor::operator=(const FilterFunctor & source)
  {
    if (this != &source)
    {
      DefaultParamHandler::operator=(source);
    }
    return *this;
  }

  FilterFunctor::~FilterFunctor() = default;

}
