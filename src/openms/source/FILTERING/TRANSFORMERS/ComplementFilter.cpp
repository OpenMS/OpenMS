// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/ComplementFilter.h>

using namespace std;

namespace OpenMS
{
  ComplementFilter::ComplementFilter() :
    FilterFunctor()
  {
    setName(ComplementFilter::getProductName());
    //value from Bioinformatics, Bern 2004
    defaults_.setValue("tolerance", 0.37, "Tolerance value as defined by Bern et al.");
    defaultsToParam_();
  }

  ComplementFilter::ComplementFilter(const ComplementFilter & source) = default;

  ComplementFilter & ComplementFilter::operator=(const ComplementFilter & source)
  {
    if (this != &source)
    {
      FilterFunctor::operator=(source);
    }
    return *this;
  }

  ComplementFilter::~ComplementFilter() = default;

}
