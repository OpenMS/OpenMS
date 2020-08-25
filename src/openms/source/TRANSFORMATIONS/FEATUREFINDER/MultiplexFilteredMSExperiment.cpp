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
#include <OpenMS/KERNEL/FeatureHandle.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexSatelliteCentroided.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteredPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteredMSExperiment.h>

#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

namespace OpenMS
{

  MultiplexFilteredMSExperiment::MultiplexFilteredMSExperiment()
  {
  }

  void MultiplexFilteredMSExperiment::addPeak(const MultiplexFilteredPeak& peak)
  {
    result_.push_back(peak);
  }

  MultiplexFilteredPeak MultiplexFilteredMSExperiment::getPeak(size_t i) const
  {
    return result_[i];
  }

  double MultiplexFilteredMSExperiment::getMZ(size_t i) const
  {
    return result_[i].getMZ();
  }

  vector<double> MultiplexFilteredMSExperiment::getMZ() const
  {
    vector<double> mz;
    mz.resize(result_.size());
    for (size_t i = 0; i < result_.size(); ++i)
    {
      mz[i] = result_[i].getMZ();
    }
    return mz;
  }

  double MultiplexFilteredMSExperiment::getRT(size_t i) const
  {
    return result_[i].getRT();
  }

  vector<double> MultiplexFilteredMSExperiment::getRT() const
  {
    vector<double> rt;
    rt.resize(result_.size());
    for (size_t i = 0; i < result_.size(); ++i)
    {
      rt[i] = result_[i].getRT();
    }
    return rt;
  }

  size_t MultiplexFilteredMSExperiment::size() const
  {
    return result_.size();
  }
  
  void MultiplexFilteredMSExperiment::writeDebugOutput(const MSExperiment& exp_picked, String debug_out) const
  {
    ConsensusMap map;
    
    // loop over peaks
    for (std::vector<MultiplexFilteredPeak>::const_iterator it_peak = result_.begin(); it_peak < result_.end(); ++it_peak)
    {
      ConsensusFeature consensus;
      
      consensus.setRT(it_peak->getRT());
      consensus.setMZ(it_peak->getMZ());
      consensus.setIntensity(1.0);
      consensus.setCharge(1);
      consensus.setQuality(1.0);
      
      std::multimap<size_t, MultiplexSatelliteCentroided > satellites = it_peak->getSatellites();
      int count = 0;
      
      // loop over satellites
      for (const auto &it_satellite : satellites)
      {
        // find indices of the peak
        size_t rt_idx = (it_satellite.second).getRTidx();
        size_t mz_idx = (it_satellite.second).getMZidx();
        
        // find peak itself
        MSExperiment::ConstIterator it_rt = exp_picked.begin();
        std::advance(it_rt, rt_idx);
        MSSpectrum::ConstIterator it_mz = it_rt->begin();
        std::advance(it_mz, mz_idx);
        
        FeatureHandle feature_handle;
        
        feature_handle.setRT(it_rt->getRT());
        feature_handle.setMZ(it_mz->getMZ());
        feature_handle.setIntensity(1.0);
        feature_handle.setCharge(it_satellite.first);
        feature_handle.setMapIndex(count);
        
        consensus.insert(feature_handle);
        map.getColumnHeaders()[count].size++;
        
        // give the maps some names (irrelevant for the debug output)
        ConsensusMap::ColumnHeader& col = map.getColumnHeaders()[count];
        std::stringstream ss;
        ss << "satellite_" << count;
        col.label = ss.str();
        col.filename = "satellites";
       
        ++count;
      }
      
      map.push_back(consensus);
    }
    
    map.sortByPosition();
    map.applyMemberFunction(&UniqueIdInterface::setUniqueId);
    map.setExperimentType("label-free");

    ConsensusXMLFile file;
    file.store(debug_out, map);
  }

}
