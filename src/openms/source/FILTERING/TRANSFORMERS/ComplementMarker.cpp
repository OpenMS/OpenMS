// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/ComplementMarker.h>

using namespace std;

namespace OpenMS
{
  ComplementMarker::ComplementMarker() :
    PeakMarker()
  {
    setName(ComplementMarker::getProductName());
    defaults_.setValue("tolerance", 1.0, "Tolerance value as defined by Bern et al.");
    defaults_.setValue("marks", 1, "How often a peak needs to be marked to be returned");
    defaultsToParam_();
  }

  ComplementMarker::ComplementMarker(const ComplementMarker & source) = default;

  ComplementMarker::~ComplementMarker() = default;

  ComplementMarker & ComplementMarker::operator=(const ComplementMarker & source)
  {
    if (this != &source)
    {
      PeakMarker::operator=(source);
    }
    return *this;
  }

}
