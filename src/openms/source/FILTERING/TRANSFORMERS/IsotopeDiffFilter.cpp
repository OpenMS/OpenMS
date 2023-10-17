// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/IsotopeDiffFilter.h>

using namespace std;

namespace OpenMS
{
  // Bern have a different one

  IsotopeDiffFilter::IsotopeDiffFilter() :
    FilterFunctor()
  {
    setName(IsotopeDiffFilter::getProductName());
    //value from Bioinformatics, Bern 2004
    defaults_.setValue("tolerance", 0.37, "Tolerance value defined by Bern et al.");
    defaultsToParam_();
  }

  IsotopeDiffFilter::IsotopeDiffFilter(const IsotopeDiffFilter & source) = default;

  IsotopeDiffFilter & IsotopeDiffFilter::operator=(const IsotopeDiffFilter & source)
  {
    if (this != &source)
    {
      FilterFunctor::operator=(source);
    }
    return *this;
  }

  IsotopeDiffFilter::~IsotopeDiffFilter() = default;

}
