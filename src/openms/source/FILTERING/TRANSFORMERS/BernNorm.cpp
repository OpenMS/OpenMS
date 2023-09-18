// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/FILTERING/TRANSFORMERS/BernNorm.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

using namespace std;

namespace OpenMS
{
  BernNorm::BernNorm() :
    DefaultParamHandler("BernNorm")
  {
    // values from the paper
    // they should be good for GoodDiff and Complements
    // IsotopeDiffs needs lower peaks
    defaults_.setValue("C1", 28.0, "C1 value of the normalization.", {"advanced"});
    defaults_.setValue("C2", 400.0, "C2 value of the normalization.", {"advanced"});
    defaults_.setValue("threshold", 0.1, "Threshold of the Bern et al. normalization."); // i.e. what is a significant peak
    defaultsToParam_();
    c1_ = 28.0;
    c2_ = 400.0;
    th_ = 0.1;
  }

  BernNorm::~BernNorm() = default;

  BernNorm::BernNorm(const BernNorm & source) :
    DefaultParamHandler(source)
  {
  }

  BernNorm & BernNorm::operator=(const BernNorm & source)
  {
    if (this != &source)
    {
      DefaultParamHandler::operator=(source);
    }
    return *this;
  }

  void BernNorm::filterPeakSpectrum(PeakSpectrum & spectrum)
  {
    filterSpectrum(spectrum);
  }

  void BernNorm::filterPeakMap(PeakMap & exp)
  {
    for (PeakMap::Iterator it = exp.begin(); it != exp.end(); ++it)
    {
      filterSpectrum(*it);
    }
  }

}
