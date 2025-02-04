// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/FEATUREFINDER/MultiplexFilteringProfile.h>

using namespace std;

namespace OpenMS
{
  MultiplexSatelliteProfile::MultiplexSatelliteProfile(float rt, double mz, float intensity) :
    rt_(rt), mz_(mz), intensity_(intensity)
  {
  }

  float MultiplexSatelliteProfile::getRT() const
  {
    return rt_;
  }
  
  double MultiplexSatelliteProfile::getMZ() const
  {
    return mz_;
  }
  
  float MultiplexSatelliteProfile::getIntensity() const
  {
    return intensity_;
  }
  
}
