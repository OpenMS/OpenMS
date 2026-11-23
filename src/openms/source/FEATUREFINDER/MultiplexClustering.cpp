// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/FEATUREFINDER/MultiplexFiltering.h>
#include <OpenMS/FEATUREFINDER/MultiplexClustering.h>
#include <OpenMS/FEATUREFINDER/PeakWidthEstimator.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/ML/CLUSTERING/GridBasedClustering.h>

#include <algorithm>
#include <exception>

using namespace std;

namespace OpenMS
{
  namespace
  {
    /// the m/z of the peak at position size/2 of all peaks of @p exp, sorted by m/z (the upper median)
    double medianMZ(const MSExperiment& exp)
    {
      std::vector<double> mz;
      mz.reserve(exp.getSize());
      for (const MSSpectrum& spec : exp)
      {
        for (const Peak1D& peak : spec)
        {
          mz.push_back(peak.getMZ());
        }
      }
      // only the element at size/2 is needed, not the whole order
      std::nth_element(mz.begin(), mz.begin() + mz.size() / 2, mz.end());
      return mz[mz.size() / 2];
    }
  }

  MultiplexClustering::MultiplexClustering(const MSExperiment& exp_profile, const MSExperiment& exp_picked, const std::vector<std::vector<PeakPickerHiRes::PeakBoundary> >& boundaries, double rt_typical) :
    rt_typical_(rt_typical)
  {
    if (exp_picked.size() != boundaries.size())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Centroided data and the corresponding list of peak boundaries do not contain same number of spectra.");
    }
    
    // ranges of the experiment
    double mz_min = exp_profile.spectrumRanges().getMinMZ();
    double mz_max = exp_profile.spectrumRanges().getMaxMZ();
    double rt_min = exp_profile.spectrumRanges().getMinRT();
    double rt_max = exp_profile.spectrumRanges().getMaxRT();

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
    rt_scaling_ = estimator.getPeakWidth(medianMZ(exp_picked)) / rt_typical_;

  }

  MultiplexClustering::MultiplexClustering(const MSExperiment& exp, double mz_tolerance, bool mz_tolerance_unit, double rt_typical) :
    rt_typical_(rt_typical)
  {
    // ranges of the experiment
    double mz_min = exp.spectrumRanges().byMSLevel(1).getMinMZ();
    double mz_max = exp.spectrumRanges().byMSLevel(1).getMaxMZ();
    double rt_min = exp.spectrumRanges().byMSLevel(1).getMinRT();
    double rt_max = exp.spectrumRanges().byMSLevel(1).getMaxRT();

    if (!RangeMZ(0.0, 1.0e12).containsMZ({mz_min, mz_max}) ||
        !RangeRT(-1.0e12, 1.0e12).containsRT({rt_min, rt_max}) ) 
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MinMZ,MaxMZ,MinRT,MaxRT values outside of sensible value ranges. Are they uninitialized? (" + StringUtils::toStr(mz_min) + "/" + StringUtils::toStr(mz_max) + "/" + StringUtils::toStr(rt_min) + "/" + StringUtils::toStr(rt_max));
    }
    
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
    if (mz_tolerance_unit)
    {
      rt_scaling_ = (medianMZ(exp) * mz_tolerance/1000000) / rt_typical_;
    }
    else
    {
      rt_scaling_ = mz_tolerance / rt_typical_;
    }

  }

  std::vector<std::map<int, GridBasedCluster> > MultiplexClustering::cluster(const std::vector<MultiplexFilteredMSExperiment>& filter_results)
  {
    /*
      GridBasedClustering merges a cluster only with clusters centred in its own or a neighbouring grid cell, and the
      centre of a cluster (the mean of its points) stays within the m/z range of its points. Points of a filter result
      that are separated by empty grid columns (m/z) therefore never meet, and the merges within each group of them are
      exactly those of the clustering of all points of the filter result. Each group is clustered on its own, the groups
      of all patterns on several threads, and their clusters are put together under the indices of their points.
      A filter result with a point outside of the grid is clustered as a whole, which throws as before.
    */
    // Groups are separated by at least this many empty columns. One would do; the others cover the rounding of cluster
    // centres, which can lie a few ulps outside of the m/z range of their points.
    const size_t min_empty_columns = 3;

    struct Task
    {
      size_t pattern;             ///< index of the filter result
      std::vector<int> points;    ///< indices of its points (ascending; all of them if whole)
      bool whole;                 ///< all points of the filter result, on the whole grid
      size_t first_column;        ///< range of the grid columns (m/z) that the clustering of the group uses
      size_t last_column;
    };
    std::vector<Task> tasks;
    std::vector<std::vector<double>> mz_of(filter_results.size());
    std::vector<std::vector<double>> rt_of(filter_results.size());
    for (size_t pattern = 0; pattern < filter_results.size(); ++pattern)
    {
      std::vector<double>& mz = mz_of[pattern];
      std::vector<double>& rt = rt_of[pattern];
      mz = filter_results[pattern].getMZ();
      rt = filter_results[pattern].getRT();
      if (mz.empty())
      {
        continue;
      }
      bool in_grid = true;
      for (size_t i = 0; i < mz.size(); ++i)
      {
        in_grid = in_grid && mz[i] >= grid_spacing_mz_.front() && mz[i] <= grid_spacing_mz_.back() && rt[i] >= grid_spacing_rt_.front() && rt[i] <= grid_spacing_rt_.back();
      }
      if (!in_grid)
      {
        std::vector<int> all(mz.size());
        for (size_t i = 0; i < all.size(); ++i) { all[i] = static_cast<int>(i); }
        tasks.push_back({pattern, std::move(all), true, 0, 0});
        continue;
      }
      // the grid column of each point (as ClusteringGrid::getIndex(): the first spacing greater than m/z, at least 1)
      std::vector<std::pair<size_t, int>> column_point(mz.size());
      for (size_t i = 0; i < mz.size(); ++i)
      {
        column_point[i] = {static_cast<size_t>(std::upper_bound(grid_spacing_mz_.begin(), grid_spacing_mz_.end(), mz[i]) - grid_spacing_mz_.begin()), static_cast<int>(i)};
      }
      std::sort(column_point.begin(), column_point.end());
      for (size_t begin = 0; begin < column_point.size();)
      {
        size_t end = begin + 1;
        while (end < column_point.size() && column_point[end].first - column_point[end - 1].first <= min_empty_columns)
        {
          ++end;
        }
        Task task{pattern, {}, false, 0, 0};
        task.points.reserve(end - begin);
        for (size_t k = begin; k < end; ++k)
        {
          task.points.push_back(column_point[k].second);
        }
        std::sort(task.points.begin(), task.points.end());
        // one more column on each side for centres rounded beyond their points (the first column of a point is >= 1)
        task.first_column = column_point[begin].first >= 2 ? column_point[begin].first - 2 : 0;
        task.last_column = std::min(column_point[end - 1].first + 1, grid_spacing_mz_.size() - 1);
        tasks.push_back(std::move(task));
        begin = end;
      }
    }

    // progress logger
    unsigned progress = 0;
    startProgress(0, tasks.size(), "clustering filtered LC-MS data");

    const SignedSize n_tasks = static_cast<SignedSize>(tasks.size());
    std::vector<std::map<int, GridBasedCluster>> task_results(n_tasks);
    // An exception is rethrown after the loop, the one of the first group that threw.
    std::vector<std::exception_ptr> errors(n_tasks);

    #pragma omp parallel for schedule(dynamic, 1)
    for (SignedSize t = 0; t < n_tasks; ++t)
    {
      const Task& task = tasks[t];
      try
      {
        if (task.whole)
        {
          GridBasedClustering<MultiplexDistance> clustering(MultiplexDistance(rt_scaling_), mz_of[task.pattern], rt_of[task.pattern], grid_spacing_mz_, grid_spacing_rt_);
          clustering.cluster();
          //clustering.extendClustersY();
          task_results[t] = clustering.getResults();
        }
        else
        {
          std::vector<double> mz, rt;
          mz.reserve(task.points.size());
          rt.reserve(task.points.size());
          for (int point : task.points)
          {
            mz.push_back(mz_of[task.pattern][point]);
            rt.push_back(rt_of[task.pattern][point]);
          }
          // the columns of the group, at the same positions relative to each other as in the whole grid
          std::vector<double> grid_spacing_mz(grid_spacing_mz_.begin() + task.first_column, grid_spacing_mz_.begin() + task.last_column + 1);
          GridBasedClustering<MultiplexDistance> clustering(MultiplexDistance(rt_scaling_), mz, rt, std::move(grid_spacing_mz), grid_spacing_rt_);
          clustering.cluster();
          task_results[t] = clustering.getResults();
        }
      }
      catch (...)
      {
        errors[t] = std::current_exception();
      }

      #pragma omp critical (MultiplexClustering_progress)
      setProgress(++progress);
    }

    endProgress();

    for (const std::exception_ptr& error : errors)
    {
      if (error)
      {
        std::rethrow_exception(error);
      }
    }

    // the clusters of each filter result, under the indices of their points in it
    std::vector<std::map<int, GridBasedCluster> > cluster_results(filter_results.size());
    for (SignedSize t = 0; t < n_tasks; ++t)
    {
      const Task& task = tasks[t];
      std::map<int, GridBasedCluster>& clusters = cluster_results[task.pattern];
      if (task.whole)
      {
        clusters = std::move(task_results[t]);
        continue;
      }
      for (const auto& [index, cluster] : task_results[t])
      {
        std::vector<int> points;
        points.reserve(cluster.getPoints().size());
        for (int point : cluster.getPoints())
        {
          points.push_back(task.points[point]);
        }
        clusters.emplace(task.points[index],
                         GridBasedCluster(cluster.getCentre(), cluster.getBoundingBox(), std::move(points), cluster.getPropertyA(), std::vector<int>(cluster.getPropertiesB())));
      }
    }

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
