// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/NeutralLossMarker.h>

using namespace std;

namespace OpenMS
{

  NeutralLossMarker::NeutralLossMarker() :
    PeakMarker()
  {
    setName(NeutralLossMarker::getProductName());
    defaults_.setValue("marks", 1, "How often a peak must be marked to be reported");
    defaults_.setValue("tolerance", 0.2, "Tolerance in m/z direction");
    defaultsToParam_();
  }

  NeutralLossMarker::NeutralLossMarker(const NeutralLossMarker & source) = default;

  NeutralLossMarker::~NeutralLossMarker() = default;

  NeutralLossMarker & NeutralLossMarker::operator=(const NeutralLossMarker & source)
  {
    if (this != &source)
    {
      PeakMarker::operator=(source);
    }
    return *this;
  }

}
