// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/BaseSuperimposer.h>

namespace OpenMS
{

  BaseSuperimposer::BaseSuperimposer() :
    DefaultParamHandler("BaseSuperimposer"),
    ProgressLogger()
  {
  }

  BaseSuperimposer::~BaseSuperimposer() = default;

}
