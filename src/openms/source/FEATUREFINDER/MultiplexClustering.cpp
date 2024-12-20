// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/FEATUREFINDER/MultiplexFiltering.h>
#include <OpenMS/FEATUREFINDER/MultiplexClustering.h>
#include <OpenMS/FEATUREFINDER/PeakWidthEstimator.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/ML/CLUSTERING/GridBasedClustering.h>

#include <QtCore/QDir>

using namespace std;

namespace OpenMS
{

  MultiplexClustering::MultiplexClustering(const MSExperiment& exp_profile, const MSExperiment& exp_picked, const std::vector<std::vector<PeakPickerHiRes::PeakBoundary> >& boundaries, double rt_typical) :
    rt_typical_(rt_typical)
  {
    if (exp_picked.size() != boundaries.size())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Centroided data and the corresponding list of peak boundaries do not contain same number of spectra.");
    }
    
    // ranges of the experiment
    double mz_min = exp_profile.getMinMZ();
    double mz_max = exp_profile.getMaxMZ();
    double rt_min = exp_profile.getMinRT();
    double rt_max = exp_profile.getMaxRT();
    
    // extend the grid by a small absolute margin
    double mz_margin = 1e-2;
    double rt_margin = 1e-2;
    mz_min -= mz_margin; 
    mz_max += mz_margin; 
    rt_min -= rt_margin; 
    rt_max += rt_margin;
    
    // generate grid spacing
    PeakWidthEstimator estimator(exp_picked, boundaries);
    // We assume that the jitter of the peak centres are less than <scaling> times the peak width.
    // This factor ensures that two neighbouring peaks at the same RT cannot be in the same cluster.
    double scaling = 0.4;
    for (double mz = mz_min; mz < mz_max; mz = mz + scaling * estimator.getPeakWidth(mz))
    {
      grid_spacing_mz_.push_back(mz);
    }
    grid_spacing_mz_.push_back(mz_max);

    for (double rt = rt_min; rt < rt_max; rt = rt + rt_typical)
    {
      grid_spacing_rt_.push_back(rt);
    }
    grid_spacing_rt_.push_back(rt_max);

    // determine RT scaling
    std::vector<double> mz;
    
    for (const MSSpectrum& spec : exp_picked)
    {
      MSSpectrum::ConstIterator it_mz;
      for (it_mz = spec.begin(); it_mz != spec.end(); ++it_mz)
      {
        mz.push_back(it_mz->getMZ());
      }
    }
    std::sort(mz.begin(), mz.end());
    rt_scaling_ = estimator.getPeakWidth(mz[(int) mz.size() / 2]) / rt_typical_;

  }

  MultiplexClustering::MultiplexClustering(const MSExperiment& exp, double mz_tolerance, bool mz_tolerance_unit, double rt_typical) :
    rt_typical_(rt_typical)
  {
    std::cout << "Peaks: " << exp.getSize() << std::endl;
    // ranges of the experiment
    double mz_min = exp.getMinMZ();
    double mz_max = exp.getMaxMZ();
    double rt_min = exp.getMinRT();
    double rt_max = exp.getMaxRT();

    if (!RangeMZ(1e-12, 1.0e12).containsMZ({mz_min, mz_max}) ||
        !RangeRT(-1.0e12, 1.0e12).containsRT({rt_min, rt_max}) ) 
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
        "MinMZ,MaxMZ,MinRT,MaxRT values outside of sensible value ranges. Are they uninitialized? (" 
        + String(mz_min) + "/" + String(mz_max) + "/" + String(rt_min) + "/" + String(rt_max));
    }

    std::cout << "Peaks: " << exp.getSize() << std::endl;

    // extend the grid by a small absolute margin
    double mz_margin = 1e-2;
    double rt_margin = 1e-2;
    mz_min -= mz_margin; 
    mz_max += mz_margin; 
    rt_min -= rt_margin; 
    rt_max += rt_margin;
    
    // output all variables
    #ifdef DEBUG
     std::cout << "mz_min: " << mz_min << std::endl;
     std::cout << "mz_max: " << mz_max << std::endl;
     std::cout << "rt_min: " << rt_min << std::endl;
     std::cout << "rt_max: " << rt_max << std::endl;
     std::cout << "rt_typical: " << rt_typical << std::endl;
     std::cout << "mz_tolerance: " << mz_tolerance << std::endl;
     std::cout << "mz_tolerance_unit: " << mz_tolerance_unit << std::endl;
     std::cout << "rt_typical: " << rt_typical << std::endl;
     std::cout << "factor:" << 1.0 + mz_tolerance/1000000 << std::endl;
    #endif
    
    // generate grid spacing
    // We assume that the jitter of the peak centres are less than <scaling> times the user specified m/z tolerance.
    double scaling = 1.0;
    
    if (mz_tolerance_unit)
    {
      for (double mz = mz_min; mz < mz_max; mz *= 1.0 + scaling * mz_tolerance/1000000)
      {
        grid_spacing_mz_.push_back(mz);
      }
    }
    else
    {
      for (double mz = mz_min; mz < mz_max; mz += scaling * mz_tolerance)
      {
        grid_spacing_mz_.push_back(mz);
      }
    }
    grid_spacing_mz_.push_back(mz_max);

    for (double rt = rt_min; rt < rt_max; rt += rt_typical)
    {
      grid_spacing_rt_.push_back(rt);
    }
    grid_spacing_rt_.push_back(rt_max);

    // determine RT scaling
    std::vector<double> mz;
    for (const MSSpectrum& spec : exp)
    {
      for (const Peak1D& p : spec)
      {
        mz.push_back(p.getMZ());
      }
    }
    std::sort(mz.begin(), mz.end());
    if (mz_tolerance_unit)
    {
      rt_scaling_ = (mz[(int) mz.size() / 2] * mz_tolerance/1000000) / rt_typical_;
    }
    else
    {
      rt_scaling_ = mz_tolerance / rt_typical_;
    }
  }

  std::vector<std::map<int, GridBasedCluster> > MultiplexClustering::cluster(const std::vector<MultiplexFilteredMSExperiment>& filter_results)
  {
    // progress logger
    unsigned progress = 0;
    startProgress(0, filter_results.size(), "clustering filtered LC-MS data");
      
    std::vector<std::map<int, GridBasedCluster> > cluster_results;

    // loop over patterns i.e. cluster each of the corresponding filter results
    for (unsigned i = 0; i < filter_results.size(); ++i)
    {
      setProgress(++progress);
        
      GridBasedClustering<MultiplexDistance> clustering(MultiplexDistance(rt_scaling_), filter_results[i].getMZ(), filter_results[i].getRT(), grid_spacing_mz_, grid_spacing_rt_);
      clustering.cluster();
      //clustering.extendClustersY();
      cluster_results.push_back(clustering.getResults());
    }
    OPENMS_LOG_INFO << "Cluster results: " << cluster_results.size() << std::endl;

    // output cluster results
    for (unsigned i = 0; i < cluster_results.size(); ++i)
    {
      for (std::map<int, GridBasedCluster>::const_iterator it = cluster_results[i].begin(); it != cluster_results[i].end(); ++it)
      {
        OPENMS_LOG_INFO << "Cluster " << it->first << " contains " << it->second.getPoints().size() << " points." << std::endl;
      }
    }
    endProgress();

    return cluster_results;
  }

  MultiplexClustering::MultiplexDistance::MultiplexDistance(double rt_scaling)
  : rt_scaling_(rt_scaling)
  {
  }
  
  MultiplexClustering::MultiplexDistance::MultiplexDistance()
  : rt_scaling_(1)
  {
  }
  
  double MultiplexClustering::MultiplexDistance::operator()(const Point& p1, const Point& p2) const
  {
      return sqrt((p1.getX() - p2.getX())*(p1.getX() - p2.getX()) + rt_scaling_ * rt_scaling_ * (p1.getY() - p2.getY())*(p1.getY() - p2.getY()));
  }

}
