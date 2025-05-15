// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/PROCESSING/FILTERING/WindowMower.h>

using namespace std;

namespace OpenMS
{
  WindowMower::WindowMower() :
    DefaultParamHandler("WindowMower")
  {
    defaults_.setValue("windowsize", 50.0, "The size of the sliding window along the m/z axis.");
    defaults_.setValue("peakcount", 2, "The number of peaks that should be kept.");
    defaults_.setValue("movetype", "slide", "Whether sliding window (one peak steps) or jumping window (window size steps) should be used.");
    defaults_.setValidStrings("movetype", {"slide","jump"});
    defaultsToParam_();
  }

  WindowMower::~WindowMower() = default;

  WindowMower::WindowMower(const WindowMower & source) :
    DefaultParamHandler(source)
  {
  }

  WindowMower & WindowMower::operator=(const WindowMower & source)
  {
    if (this != &source)
    {
      DefaultParamHandler::operator=(source);
    }
    return *this;
  }

  void WindowMower::filterPeakSpectrum(PeakSpectrum & spectrum)
  {
    bool sliding = param_.getValue("movetype").toString() == "slide" ? true : false;

    if (sliding)
    {
      filterPeakSpectrumForTopNInSlidingWindow(spectrum);
    }
    else
    {
      filterPeakSpectrumForTopNInJumpingWindow(spectrum);
    }
  }

  void WindowMower::filterPeakMap(PeakMap & exp)
  {
    bool sliding = param_.getValue("movetype").toString() == "slide" ? true : false;
    for (PeakMap::Iterator it = exp.begin(); it != exp.end(); ++it)
    {
      if (sliding)
      {
        filterPeakSpectrumForTopNInSlidingWindow(*it);
      }
      else
      {
        filterPeakSpectrumForTopNInJumpingWindow(*it);
      }
    }
  }

}
