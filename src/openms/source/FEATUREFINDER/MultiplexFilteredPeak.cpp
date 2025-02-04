// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FEATUREFINDER/MultiplexSatelliteCentroided.h>
#include <OpenMS/FEATUREFINDER/MultiplexSatelliteProfile.h>
#include <OpenMS/FEATUREFINDER/MultiplexFilteredPeak.h>

#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

namespace OpenMS
{
  MultiplexFilteredPeak::MultiplexFilteredPeak(double mz, float rt, size_t mz_idx, size_t rt_idx) :
    mz_(mz), rt_(rt), mz_idx_(mz_idx), rt_idx_(rt_idx)
  {
  }

  double MultiplexFilteredPeak::getMZ() const
  {
    return mz_;
  }

  float MultiplexFilteredPeak::getRT() const
  {
    return rt_;
  }

  size_t MultiplexFilteredPeak::getMZidx() const
  {
    return mz_idx_;
  }

  size_t MultiplexFilteredPeak::getRTidx() const
  {
    return rt_idx_;
  }

  void MultiplexFilteredPeak::addSatellite(size_t rt_idx, size_t mz_idx, size_t pattern_idx)
  {
    satellites_.insert(std::make_pair(pattern_idx, MultiplexSatelliteCentroided(rt_idx, mz_idx)));
  }
  
  void MultiplexFilteredPeak::addSatellite(const MultiplexSatelliteCentroided& satellite, size_t pattern_idx)
  {
    satellites_.insert(std::make_pair(pattern_idx, satellite));
  }
  
  void MultiplexFilteredPeak::addSatelliteProfile(float rt, double mz, float intensity, size_t pattern_idx)
  {
    satellites_profile_.insert(std::make_pair(pattern_idx, MultiplexSatelliteProfile(rt, mz, intensity)));
  }
  
  void MultiplexFilteredPeak::addSatelliteProfile(const MultiplexSatelliteProfile& satellite, size_t pattern_idx)
  {
    satellites_profile_.insert(std::make_pair(pattern_idx, satellite));
  }
  
  bool MultiplexFilteredPeak::checkSatellite(size_t rt_idx, size_t mz_idx) const
  {
    for (const auto &satellite_it : satellites_)
    {
      if (((satellite_it.second).getRTidx() == rt_idx) && ((satellite_it.second).getMZidx() == mz_idx))
      {
        return true;
      }
    }
    
    return false;
  }
  
  const std::multimap<size_t, MultiplexSatelliteCentroided >& MultiplexFilteredPeak::getSatellites() const
  {
    return satellites_;
  }
  
  const std::multimap<size_t, MultiplexSatelliteProfile >& MultiplexFilteredPeak::getSatellitesProfile() const
  {
    return satellites_profile_;
  }
  
  size_t MultiplexFilteredPeak::size() const
  {
    return satellites_.size();
  }
  
  size_t MultiplexFilteredPeak::sizeProfile() const
  {
    return satellites_profile_.size();
  }
}
