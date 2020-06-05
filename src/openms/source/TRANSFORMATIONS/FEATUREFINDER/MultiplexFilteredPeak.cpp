// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexSatelliteCentroided.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexSatelliteProfile.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteredPeak.h>

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
