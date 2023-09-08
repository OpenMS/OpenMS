// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/TICFilter.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

using namespace std;

namespace OpenMS
{
  TICFilter::TICFilter() :
    FilterFunctor()
  {
    setName(TICFilter::getProductName());
    defaults_.setValue("window", 5, "Windowing parameter which defines the windows size");
    defaultsToParam_();
  }

  TICFilter::TICFilter(const TICFilter & source) = default;

  TICFilter & TICFilter::operator=(const TICFilter & source)
  {
    if (this != &source)
    {
      FilterFunctor::operator=(source);
    }
    return *this;
  }

  TICFilter::~TICFilter() = default;

}
