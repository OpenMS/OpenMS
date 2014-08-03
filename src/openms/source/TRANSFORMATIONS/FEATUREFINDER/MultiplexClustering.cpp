// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilterResult.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SplinePackage.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SplineSpectrum.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFiltering.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexClustering.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/MATH/STATISTICS/LinearRegression.h>
#include <OpenMS/COMPARISON/CLUSTERING/MultiplexCluster.h>
#include <OpenMS/COMPARISON/CLUSTERING/GridClustering.h>

#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

namespace OpenMS
{

  MultiplexClustering::MultiplexClustering(MSExperiment<Peak1D> exp_profile, MSExperiment<Peak1D> exp_picked, std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries, double rt_typical, double rt_minimum, bool debug) :
    rt_typical_(rt_typical), rt_minimum_(rt_minimum), debug_(debug)
  {
    if (exp_picked.size() != boundaries.size())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Centroided data and the corresponding list of peak boundaries do not contain same number of spectra.");
    }
    
    // ranges of the experiment
    double mz_min = exp_profile.getMinMZ();
    double mz_max = exp_profile.getMaxMZ();
    double rt_min = exp_profile.getMinRT();
    double rt_max = exp_profile.getMaxRT();
    
    // generate grid spacing
    PeakWidthEstimator estimator(exp_picked, boundaries);
    // We assume that the jitter of the peak centres are less than <scaling> times the peak width.
    // This factor ensures that two neighbouring peaks at the same RT cannot be in the same cluster.
    double scaling = 0.2;
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
    MSExperiment<Peak1D>::Iterator it_rt;
    for (it_rt = exp_picked.begin(); it_rt < exp_picked.end(); ++it_rt)
    {
      MSSpectrum<Peak1D>::Iterator it_mz;
      for (it_mz = it_rt->begin(); it_mz < it_rt->end(); ++it_mz)
      {
        mz.push_back(it_mz->getMZ());
      }
    }
    std::sort(mz.begin(), mz.end());
    // RT scaling = peak width at the median of the m/z distribuation / RT threshold
    rt_scaling_ = estimator.getPeakWidth(mz[(int) mz.size() / 2]) / rt_typical_;

  }

  std::vector<std::map<int, MultiplexCluster> > MultiplexClustering::cluster(std::vector<MultiplexFilterResult> filter_results)
  {
    std::vector<std::map<int, MultiplexCluster> > cluster_results;

    // loop over patterns i.e. cluster each of the corresponding filter results
    for (unsigned i = 0; i < filter_results.size(); ++i)
    {
        
      GridClustering<MultiplexDistance> clustering(MultiplexDistance(rt_scaling_), filter_results[i].getMZ(), filter_results[i].getRT(), grid_spacing_mz_, grid_spacing_rt_);
      clustering.cluster();
      //clustering.extendClustersY();
      clustering.removeSmallClustersY(rt_minimum_);
      cluster_results.push_back(clustering.getResults());

      // debug output
      vector<DebugPoint> debug_clustered;
      if (debug_)
      {
        std::map<int, MultiplexCluster> cluster_result = clustering.getResults();
        MultiplexFilterResult filter_result = filter_results[i];

        int cluster_id = 0;
        for (std::map<int, MultiplexCluster>::iterator it = cluster_result.begin(); it != cluster_result.end(); ++it)
        {
          std::vector<int> points = (it->second).getPoints();
          for (std::vector<int>::iterator it2 = points.begin(); it2 != points.end(); ++it2)
          {
            DebugPoint data_point;
            data_point.rt = filter_result.getRT(*it2);
            data_point.mz = filter_result.getMZ(*it2);
            data_point.cluster = cluster_id;
            debug_clustered.push_back(data_point);
          }
          ++cluster_id;
        }
        writeDebug(debug_clustered, i);
      }

    }

    return cluster_results;
  }

  MultiplexClustering::PeakWidthEstimator::PeakWidthEstimator(MSExperiment<Peak1D> exp_picked, std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries)
  {
    if (exp_picked.size() != boundaries.size())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Centroided data and the corresponding list of peak boundaries do not contain same number of spectra.");
    }

    std::map<double, double> mz_peak_width;
    MSExperiment<Peak1D>::Iterator it_rt;
    vector<vector<PeakPickerHiRes::PeakBoundary> >::const_iterator it_rt_boundaries;
    for (it_rt = exp_picked.begin(), it_rt_boundaries = boundaries.begin();
         it_rt < exp_picked.end() && it_rt_boundaries < boundaries.end();
         ++it_rt, ++it_rt_boundaries)
    {
      MSSpectrum<Peak1D>::Iterator it_mz;
      vector<PeakPickerHiRes::PeakBoundary>::const_iterator it_mz_boundary;
      for (it_mz = it_rt->begin(), it_mz_boundary = it_rt_boundaries->begin();
           it_mz < it_rt->end(), it_mz_boundary < it_rt_boundaries->end();
           ++it_mz, ++it_mz_boundary)
      {
        mz_peak_width[it_mz->getMZ()] = (*it_mz_boundary).mz_max - (*it_mz_boundary).mz_min;
      }
    }

    std::vector<double> mz;
    std::vector<double> peak_width;
    for (std::map<double, double>::iterator it = mz_peak_width.begin(); it != mz_peak_width.end(); ++it)
    {
      mz.push_back(it->first);
      peak_width.push_back(it->second);
    }
    mz_min_ = mz.front();
    mz_max_ = mz.back();

    Math::LinearRegression linreg;
    linreg.computeRegression(0.95, mz.begin(), mz.end(), peak_width.begin());
    slope_ = linreg.getSlope();
    intercept_ = linreg.getIntercept();
  }

  double MultiplexClustering::PeakWidthEstimator::getPeakWidth(double mz)
  {
    double width;

    if (mz < mz_min_)
    {
      width = slope_ * mz_min_ + intercept_;
    }
    else if (mz > mz_max_)
    {
      width = slope_ * mz_max_ + intercept_;
    }
    else
    {
      width = slope_ * mz + intercept_;
    }

    if (width < 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Estimated peak width is negative.", "");
    }

    return width;
  }

  MultiplexClustering::MultiplexDistance::MultiplexDistance(double rt_scaling)
  : rt_scaling_(rt_scaling)
  {
  }
  
  MultiplexClustering::MultiplexDistance::MultiplexDistance()
  : rt_scaling_(1)
  {
  }
  
  double MultiplexClustering::MultiplexDistance::operator()(Point p1, Point p2)
  {
      return sqrt((p1.getX() - p2.getX())*(p1.getX() - p2.getX()) + rt_scaling_ * rt_scaling_ * (p1.getY() - p2.getY())*(p1.getY() - p2.getY()));
  }

  void MultiplexClustering::writeDebug(vector<DebugPoint> points, int pattern) const
  {
    // fill consensus map
    ConsensusMap map;
    for (std::vector<DebugPoint>::const_iterator it = points.begin(); it != points.end(); ++it)
    {
      ConsensusFeature consensus;
      consensus.setRT((*it).rt);
      consensus.setMZ((*it).mz);
      consensus.setIntensity((*it).cluster);
      consensus.setCharge(1);          // dummy
      consensus.setMetaValue("color", getColour((*it).cluster));
      consensus.setMetaValue("Cluster ID", (*it).cluster);

      FeatureHandle feature;
      feature.setRT((*it).rt);
      feature.setMZ((*it).mz);
      feature.setUniqueId((*it).cluster);
      consensus.insert(feature);

      map.getFileDescriptions()[0].size++;
      map.push_back(consensus);
    }

    ConsensusMap::FileDescription& desc = map.getFileDescriptions()[0];
    desc.filename = "debug";
    desc.label = "Cluster";

    map.sortByPosition();
    map.applyMemberFunction(&UniqueIdInterface::setUniqueId);
    map.setExperimentType("multiplex");

    // write consensus file
    ConsensusXMLFile file;
    String file_name = "debug_clustered_";
    file_name = file_name + pattern + ".consensusXML";
    file.store(file_name, map);
  }

  String MultiplexClustering::getColour(int c) const
  {
    // 35 + 15 HTML colors
    static const String colours[] =
    {
      "#FAEBD7", "#7FFFD4", "#FFE4C4", "#8A2BE2", "#A52A2A",
      "#DEB887", "#5F9EA0", "#7FFF00", "#D2691E", "#FF7F50",
      "#6495ED", "#DC143C", "#00008B", "#008B8B", "#B8860B",
      "#A9A9A9", "#006400", "#BDB76B", "#8B008B", "#556B2F",
      "#FF8C00", "#8B0000", "#E9967A", "#8FBC8F", "#483D8B",
      "#2F4F4F", "#00CED1", "#9400D3", "#FF1493", "#B22222",
      "#FFD700", "#DAA520", "#008000", "#ADFF2F", "#C71585",

      "#00FFFF", "#000000", "#0000FF", "#FF00FF", "#008000",
      "#808080", "#00FF00", "#800000", "#000080", "#808000",
      "#800080", "#FF0000", "#C0C0C0", "#008080", "#FFFF00",
    };

    return colours[c % 50];
  }

}
