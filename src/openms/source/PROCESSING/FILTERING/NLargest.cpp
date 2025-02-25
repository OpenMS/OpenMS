// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/PROCESSING/FILTERING/NLargest.h>

using namespace std;

namespace OpenMS
{

  NLargest::NLargest() :
    DefaultParamHandler("NLargest")
  {
    init_();
  }

  NLargest::NLargest(UInt n) :
    DefaultParamHandler("NLargest")
  {
    init_();
    // after initialising with the default value, use the provided n
    param_.setValue("n", n);
    updateMembers_();
  }

  void NLargest::init_()
  {
    defaults_.setValue("n", 200, "The number of peaks to keep");
    defaultsToParam_();
  }

  NLargest::~NLargest() = default;

  NLargest::NLargest(const NLargest & source) :
    DefaultParamHandler(source)
  {
    updateMembers_();
  }

  NLargest & NLargest::operator=(const NLargest & source)
  {
    if (this != &source)
    {
      DefaultParamHandler::operator=(source);
      updateMembers_();
    }
    return *this;
  }

  void NLargest::filterPeakSpectrum(PeakSpectrum & spectrum)
  {
    filterSpectrum(spectrum);
  }

  void NLargest::filterPeakMap(PeakMap & exp)
  {
    for (PeakMap::Iterator it = exp.begin(); it != exp.end(); ++it)
    {
      filterSpectrum(*it);
    }
  }

  void NLargest::updateMembers_()
  {
    peakcount_ = (UInt)param_.getValue("n");
  }

}
