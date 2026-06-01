// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger, Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/MRMRTNormalizer.h>
#include <OpenMS/ML/REGRESSION/LinearRegression.h>
#include <OpenMS/CONCEPT/LogStream.h> // OPENMS_LOG_DEBUG
#include <OpenMS/ML/RANSAC/RANSAC.h> // RANSAC algorithm

#include <numeric>
#include <boost/math/special_functions/erf.hpp>

namespace OpenMS
{

  std::vector<std::pair<double, double> > MRMRTNormalizer::removeOutliersRANSAC(
      const std::vector<std::pair<double, double> >& pairs, double rsq_limit,
      double coverage_limit, size_t max_iterations, double max_rt_threshold, size_t sampling_size)
  {
    size_t n = sampling_size;
    size_t k = (size_t)max_iterations;
    double t = max_rt_threshold*max_rt_threshold;
    size_t d = (size_t)(coverage_limit*pairs.size());

    if (n < 5)
    {
      throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "UnableToFit-LinearRegression-RTNormalizer", "WARNING: RANSAC: " + 
          String(n) + " sampled RT peptides is below limit of 5 peptides required for the RANSAC outlier detection algorithm.");
    }

    if (pairs.size() < 30)
    {
      throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "UnableToFit-LinearRegression-RTNormalizer", "WARNING: RANSAC: " + 
          String(pairs.size()) + " input RT peptides is below limit of 30 peptides required for the RANSAC outlier detection algorithm.");
    }

    Math::RANSAC<Math::RansacModelLinear> r;
    std::vector<std::pair<double, double> > new_pairs = r.ransac(pairs, n, k, t, d);
    double bestrsq = Math::RansacModelLinear::rm_rsq_impl(new_pairs.begin(), new_pairs.end());

    if (bestrsq < rsq_limit)
    {
      throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "UnableToFit-LinearRegression-RTNormalizer", "WARNING: rsq: " +
          String(bestrsq) + " is below limit of " +
          String(rsq_limit) +
          ". Validate assays for RT-peptides and adjust the limit for rsq or coverage.");
    }

    if (new_pairs.size() < d)
    {
      throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "UnableToFit-LinearRegression-RTNormalizer", "WARNING: number of data points: " +
          String(new_pairs.size()) +
          " is below limit of " + String(d) +
          ". Validate assays for RT-peptides and adjust the limit for rsq or coverage.");
    }

    return new_pairs;
  }

  int MRMRTNormalizer::jackknifeOutlierCandidate_(const std::vector<double>& x, const std::vector<double>& y)
  {
    // Returns candidate outlier: A linear regression and rsq is calculated for
    // the data points with one removed pair. The combination resulting in
    // highest rsq is considered corresponding to the outlier candidate. The
    // corresponding iterator position is then returned.
    std::vector<double> x_tmp, y_tmp, rsq_tmp;

    for (Size i = 0; i < x.size(); i++)
    {
      x_tmp = x;
      y_tmp = y;
      x_tmp.erase(x_tmp.begin() + i);
      y_tmp.erase(y_tmp.begin() + i);

      Math::LinearRegression lin_reg;
      lin_reg.computeRegression(0.95, x_tmp.begin(), x_tmp.end(), y_tmp.begin());

      rsq_tmp.push_back(lin_reg.getRSquared());
    }
    return max_element(rsq_tmp.begin(), rsq_tmp.end()) - rsq_tmp.begin();
  }

  int MRMRTNormalizer::residualOutlierCandidate_(const std::vector<double>& x, const std::vector<double>& y)
  {
    // Returns candidate outlier: A linear regression and residuals are calculated for
    // the data points. The one with highest residual error is selected as the outlier candidate. The
    // corresponding iterator position is then returned.
    Math::LinearRegression lin_reg;
    lin_reg.computeRegression(0.95, x.begin(), x.end(), y.begin());

    std::vector<double> residuals;

    for (Size i = 0; i < x.size(); i++)
    {
      double residual = fabs(y[i] - (lin_reg.getIntercept() + (lin_reg.getSlope() * x[i])));
      residuals.push_back(residual);
    }

    return max_element(residuals.begin(), residuals.end()) - residuals.begin();
  }

  std::vector<std::pair<double, double> > MRMRTNormalizer::removeOutliersIterative(
      const std::vector<std::pair<double, double> >& pairs, double rsq_limit,
      double coverage_limit, bool use_chauvenet, const std::string& method)
  {
    if (pairs.size() < 3)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Need at least 3 data points to remove outliers for the regression.");
    }

    // Removes outliers from vector of pairs until upper rsq and lower coverage limits are reached.
    std::vector<double> x, y;
    double confidence_interval = 0.95;

    std::vector<std::pair<double, double> > pairs_corrected;

    for (auto it = pairs.begin(); it != pairs.end(); ++it)
    {
      x.push_back(it->first);
      y.push_back(it->second);
      OPENMS_LOG_DEBUG << "RT Normalization pairs: " << it->first << " : " << it->second << std::endl;
    }

    double rsq;
    rsq = 0;

    while (x.size() >= coverage_limit * pairs.size() && rsq < rsq_limit)
    {
      Math::LinearRegression lin_reg;
      lin_reg.computeRegression(confidence_interval, x.begin(), x.end(), y.begin());

      rsq = lin_reg.getRSquared();

      if (rsq < rsq_limit)
      {
        std::vector<double> residuals;

        // calculate residuals
        for (auto it = pairs.begin(); it != pairs.end(); ++it)
        {
          double intercept = lin_reg.getIntercept();
          double slope = (double)lin_reg.getSlope();
          residuals.push_back(fabs(it->second - (intercept + it->first * slope)));
          OPENMS_LOG_DEBUG << " RT Normalization residual is " << residuals.back() << std::endl;
        }

        int pos;

        if (method == "iter_jackknife")
        {
          // get candidate outlier: removal of which datapoint results in best rsq?
          pos = jackknifeOutlierCandidate_(x, y);
        }
        else if (method == "iter_residual")
        {
          // get candidate outlier: removal of datapoint with largest residual?
          pos = residualOutlierCandidate_(x, y);
        }
        else
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            String("Method ") + method + " is not a valid method for removeOutliersIterative");
        }

        // remove if residual is an outlier according to Chauvenet's criterion
        // or if testing is turned off
        OPENMS_LOG_DEBUG << " Got outlier candidate " << pos << "(" << x[pos] << " / " << y[pos] << std::endl;
        if (!use_chauvenet || chauvenet(residuals, pos))
        {
          x.erase(x.begin() + pos);
          y.erase(y.begin() + pos);
        }
        else
        {
          break;
        }
      }
      else
      {
        break;
      }
    }

    if (rsq < rsq_limit)
    {
      // If the rsq is below the limit, this is an indication that something went wrong!
      throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "UnableToFit-LinearRegression-RTNormalizer", "WARNING: rsq: " +
          String(rsq) + " is below limit of " +
          String(rsq_limit) +
          ". Validate assays for RT-peptides and adjust the limit for rsq or coverage.");
    }

    for (Size i = 0; i < x.size(); i++)
    {
      pairs_corrected.emplace_back(x[i], y[i]);
    }

#ifdef DEBUG_MRMRTNORMALIZER
    std::cout << "=======STARTPOINTS=======" << std::endl;
    for (std::vector<std::pair<double, double> >::iterator it = pairs_corrected.begin(); it != pairs_corrected.end(); ++it)
    {
      std::cout << it->first << "\t" << it->second << std::endl;
    }
    std::cout << "=======ENDPOINTS=======" << std::endl;
#endif

    return pairs_corrected;
  }

  bool MRMRTNormalizer::chauvenet(const std::vector<double>& residuals, int pos)
  {
    double criterion = 1.0 / (2 * residuals.size());
    double prob = MRMRTNormalizer::chauvenet_probability(residuals, pos);

    OPENMS_LOG_DEBUG << " Chauvinet testing " << prob << " < " << criterion << std::endl;
    return prob < criterion;

  }

  double MRMRTNormalizer::chauvenet_probability(const std::vector<double>& residuals, int pos)
  {
    double mean = std::accumulate(residuals.begin(), residuals.end(), 0.0) / residuals.size();
    double stdev = std::sqrt(
        std::inner_product(residuals.begin(), residuals.end(), residuals.begin(), 0.0
          ) / residuals.size() - mean * mean);

    double d = fabs(residuals[pos] - mean) / stdev;
    d /= pow(2.0, 0.5);
    double prob = std::erfc(d);

    return prob;
  }

  bool MRMRTNormalizer::computeBinnedCoverage(const std::pair<double,double> & rtRange, 
      const std::vector<std::pair<double, double> > & pairs, int nrBins, 
      int minPeptidesPerBin, int minBinsFilled)
  {
    std::vector<int> binCounter(nrBins, 0);
    for (std::vector<std::pair<double, double> >::const_iterator pair_it = pairs.begin(); pair_it != pairs.end(); ++pair_it)
    {
      double normRT = (pair_it->second - rtRange.first) / (rtRange.second - rtRange.first); // compute a value between [0,1)
      normRT *= nrBins;
      int bin = (int)normRT;
      if (bin >= nrBins)
      {
        // this should never happen, but just to make sure
        std::cerr << "MRMRTNormalizer::computeBinnedCoverage : computed bin was too large (" << 
          bin << "), setting it to the maximum of " << nrBins - 1 << std::endl;
        bin = nrBins - 1;
      }
      binCounter[ bin ]++;
    }

    int binsFilled = 0;
    for (Size i = 0; i < binCounter.size(); i++)
    {
      OPENMS_LOG_DEBUG <<" In bin " << i << " out of " << binCounter.size() << 
        " we have " << binCounter[i] << " peptides " << std::endl;
      if (binCounter[i] >= minPeptidesPerBin) 
      {
        binsFilled++;
      }
    }

    return (binsFilled >= minBinsFilled);
  }

}

