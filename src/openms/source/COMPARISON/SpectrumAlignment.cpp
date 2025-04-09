// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/COMPARISON/SpectrumAlignment.h>

using namespace std;

namespace OpenMS
{
  SpectrumAlignment::SpectrumAlignment() :
    DefaultParamHandler("SpectrumAlignment")
  {
    defaults_.setValue("tolerance", 0.3, "Defines the absolute (in Da) or relative (in ppm) tolerance");
    defaults_.setValue("is_relative_tolerance", "false", "If true, the 'tolerance' is interpreted as ppm-value");
    defaults_.setValidStrings("is_relative_tolerance", {"true","false"});
    defaultsToParam_();
  }

  SpectrumAlignment::SpectrumAlignment(const SpectrumAlignment & source) = default;

  SpectrumAlignment::~SpectrumAlignment() = default;

  SpectrumAlignment & SpectrumAlignment::operator=(const SpectrumAlignment & source)
  {
    if (this != &source)
    {
      DefaultParamHandler::operator=(source);
    }
    return *this;
  }

}
