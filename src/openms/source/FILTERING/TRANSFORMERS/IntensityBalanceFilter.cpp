// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/IntensityBalanceFilter.h>

using namespace std;

namespace OpenMS
{
  IntensityBalanceFilter::IntensityBalanceFilter() :
    FilterFunctor()
  {
    check_defaults_ = false;
    setName(IntensityBalanceFilter::getProductName());
    defaultsToParam_();
  }

  IntensityBalanceFilter::IntensityBalanceFilter(const IntensityBalanceFilter & source) :
    FilterFunctor(source)
  {
    check_defaults_ = false;
  }

  IntensityBalanceFilter & IntensityBalanceFilter::operator=(const IntensityBalanceFilter & source)
  {
    if (this != &source)
    {
      FilterFunctor::operator=(source);
    }
    return *this;
  }

  IntensityBalanceFilter::~IntensityBalanceFilter() = default;

}
