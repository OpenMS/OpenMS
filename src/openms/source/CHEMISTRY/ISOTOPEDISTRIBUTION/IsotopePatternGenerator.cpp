// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Nikos Patikos $
// $Authors: Nikos Patikos $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopePatternGenerator.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <cmath>
#include <fstream>

using namespace std;

namespace OpenMS
{
  IsotopePatternGenerator::IsotopePatternGenerator(double probability_cutoff) :
    min_prob_(probability_cutoff)
  {
  }

  IsotopePatternGenerator::IsotopePatternGenerator() :
    min_prob_(1e-15)
  {
  }
  
  IsotopePatternGenerator::~IsotopePatternGenerator() = default;

}
