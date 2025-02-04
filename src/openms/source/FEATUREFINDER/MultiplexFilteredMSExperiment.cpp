// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/FeatureHandle.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FEATUREFINDER/MultiplexSatelliteCentroided.h>
#include <OpenMS/FEATUREFINDER/MultiplexFilteredPeak.h>
#include <OpenMS/FEATUREFINDER/MultiplexFilteredMSExperiment.h>

#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

namespace OpenMS
{

  MultiplexFilteredMSExperiment::MultiplexFilteredMSExperiment() = default;

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
  
  void MultiplexFilteredMSExperiment::writeDebugOutput(const MSExperiment& exp_picked, const String& debug_out) const
  {
    ConsensusMap map;
    
    // loop over peaks
    for (const MultiplexFilteredPeak& peak : result_)
    {
      ConsensusFeature consensus;
      
      consensus.setRT(peak.getRT());
      consensus.setMZ(peak.getMZ());
      consensus.setIntensity(1.0);
      consensus.setCharge(1);
      consensus.setQuality(1.0);
      
      std::multimap<size_t, MultiplexSatelliteCentroided > satellites = peak.getSatellites();
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

    FileHandler().storeConsensusFeatures(debug_out, map);
  }

}
