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
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFiltering.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexClustering.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/PeakWidthEstimator.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/COMPARISON/CLUSTERING/GridBasedClustering.h>

#include<QDir>

using namespace std;

namespace OpenMS
{

  MultiplexClustering::MultiplexClustering(const MSExperiment& exp_profile, const MSExperiment& exp_picked, const std::vector<std::vector<PeakPickerHiRes::PeakBoundary> >& boundaries, double rt_typical, double rt_minimum) :
    rt_typical_(rt_typical), rt_minimum_(rt_minimum)
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
    MSExperiment::ConstIterator it_rt;
    for (it_rt = exp_picked.begin(); it_rt != exp_picked.end(); ++it_rt)
    {
      MSSpectrum::ConstIterator it_mz;
      for (it_mz = it_rt->begin(); it_mz != it_rt->end(); ++it_mz)
      {
        mz.push_back(it_mz->getMZ());
      }
    }
    std::sort(mz.begin(), mz.end());
    rt_scaling_ = estimator.getPeakWidth(mz[(int) mz.size() / 2]) / rt_typical_;

  }

  MultiplexClustering::MultiplexClustering(const MSExperiment& exp, double mz_tolerance, bool mz_tolerance_unit, double rt_typical, double rt_minimum) :
    rt_typical_(rt_typical), rt_minimum_(rt_minimum)
  {
    // ranges of the experiment
    double mz_min = exp.getMinMZ();
    double mz_max = exp.getMaxMZ();
    double rt_min = exp.getMinRT();
    double rt_max = exp.getMaxRT();
    
    // extend the grid by a small absolute margin
    double mz_margin = 1e-2;
    double rt_margin = 1e-2;
    mz_min -= mz_margin; 
    mz_max += mz_margin; 
    rt_min -= rt_margin; 
    rt_max += rt_margin;
    
    // generate grid spacing
    // We assume that the jitter of the peak centres are less than <scaling> times the user specified m/z tolerance.
    double scaling = 1.0;
    
    if (mz_tolerance_unit)
    {
      for (double mz = mz_min; mz < mz_max; mz = mz * (1 + scaling * mz_tolerance/1000000))
      {
        grid_spacing_mz_.push_back(mz);
      }
    }
    else
    {
      for (double mz = mz_min; mz < mz_max; mz = mz + scaling * mz_tolerance)
      {
        grid_spacing_mz_.push_back(mz);
      }
    }
    grid_spacing_mz_.push_back(mz_max);

    for (double rt = rt_min; rt < rt_max; rt = rt + rt_typical)
    {
      grid_spacing_rt_.push_back(rt);
    }
    grid_spacing_rt_.push_back(rt_max);

    // determine RT scaling
    std::vector<double> mz;
    MSExperiment::ConstIterator it_rt;
    for (it_rt = exp.begin(); it_rt != exp.end(); ++it_rt)
    {
      MSSpectrum::ConstIterator it_mz;
      for (it_mz = it_rt->begin(); it_mz != it_rt->end(); ++it_mz)
      {
        mz.push_back(it_mz->getMZ());
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
  
  double MultiplexClustering::MultiplexDistance::operator()(Point p1, Point p2)
  {
      return sqrt((p1.getX() - p2.getX())*(p1.getX() - p2.getX()) + rt_scaling_ * rt_scaling_ * (p1.getY() - p2.getY())*(p1.getY() - p2.getY()));
  }

}
