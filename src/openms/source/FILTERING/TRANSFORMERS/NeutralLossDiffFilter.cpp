// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/NeutralLossDiffFilter.h>

using namespace std;

namespace OpenMS
{
  NeutralLossDiffFilter::NeutralLossDiffFilter() :
    FilterFunctor()
  {
    setName(NeutralLossDiffFilter::getProductName());
    //value from Bioinformatics, Bern 2004
    defaults_.setValue("tolerance", 0.37, "Tolerance value defined by Bern et al.");
    defaultsToParam_();
  }

  NeutralLossDiffFilter::NeutralLossDiffFilter(const NeutralLossDiffFilter & source) = default;

  NeutralLossDiffFilter & NeutralLossDiffFilter::operator=(const NeutralLossDiffFilter & source)
  {
    if (this != &source)
    {
      FilterFunctor::operator=(source);
    }
    return *this;
  }

  NeutralLossDiffFilter::~NeutralLossDiffFilter() = default;

}
