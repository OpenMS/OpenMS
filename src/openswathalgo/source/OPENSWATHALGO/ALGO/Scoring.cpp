// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/OPENSWATHALGO/ALGO/Scoring.h>
#include <OpenMS/OPENSWATHALGO/Macros.h>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <Eigen/Core>
#include <Eigen/Dense>

namespace OpenSwath::Scoring
{
    namespace detail
    {
      // Internal helper to avoid deprecation warnings between our own functions
      static void normalize_sum_impl(double* x, unsigned int n)
      {
        Eigen::Map<Eigen::VectorXd> v(x, n);
        double s = v.sum();
        if (s != 0.0) { v /= s; }
      }
    }

    void normalize_sum(double x[], unsigned int n)
    {
      detail::normalize_sum_impl(x, n);
    }

    void normalize_sum(std::vector<double>& x)
    {
      detail::normalize_sum_impl(x.data(), static_cast<unsigned int>(x.size()));
    }

    double NormalizedManhattanDist(double x[], double y[], int n)
    {
      OPENSWATH_PRECONDITION(n > 0, "Need at least one element");

      detail::normalize_sum_impl(x, n);
      detail::normalize_sum_impl(y, n);
      Eigen::Map<Eigen::VectorXd> vx(x, n);
      Eigen::Map<Eigen::VectorXd> vy(y, n);
      return (vx - vy).cwiseAbs().sum() / n;
    }

    double NormalizedManhattanDist(std::vector<double>& x, std::vector<double>& y)
    {
      OPENSWATH_PRECONDITION(x.size() == y.size() && !x.empty(), "Both vectors need to have the same non-zero length");
      detail::normalize_sum_impl(x.data(), static_cast<unsigned int>(x.size()));
      detail::normalize_sum_impl(y.data(), static_cast<unsigned int>(y.size()));
      Eigen::Map<Eigen::VectorXd> vx(x.data(), x.size());
      Eigen::Map<Eigen::VectorXd> vy(y.data(), y.size());
      return (vx - vy).cwiseAbs().sum() / x.size();
    }

    double RootMeanSquareDeviation(double x[], double y[], int n)
    {
      OPENSWATH_PRECONDITION(n > 0, "Need at least one element");

      Eigen::Map<const Eigen::VectorXd> vx(x, n);
      Eigen::Map<const Eigen::VectorXd> vy(y, n);
      return std::sqrt((vx - vy).squaredNorm() / n);
    }

    double RootMeanSquareDeviation(const std::vector<double>& x, const std::vector<double>& y)
    {
      OPENSWATH_PRECONDITION(x.size() == y.size() && !x.empty(), "Both vectors need to have the same non-zero length");
      Eigen::Map<const Eigen::VectorXd> vx(x.data(), x.size());
      Eigen::Map<const Eigen::VectorXd> vy(y.data(), y.size());
      return std::sqrt((vx - vy).squaredNorm() / x.size());
    }

    double SpectralAngle(double x[], double y[], int n)
    {
      OPENSWATH_PRECONDITION(n > 0, "Need at least one element");

      Eigen::Map<const Eigen::VectorXd> vx(x, n);
      Eigen::Map<const Eigen::VectorXd> vy(y, n);
      double x_len = vx.norm();
      double y_len = vy.norm();

      // normalise, avoiding a divide by zero. See unit tests for what happens
      // when one of the vectors has a length of zero.
      double denominator = x_len * y_len;
      double theta = (denominator == 0) ? 0.0 : vx.dot(vy) / denominator;

      // clip to range [-1, 1] to save acos blowing up
      theta = std::max(-1.0, std::min(1.0, theta));

      return std::acos(theta);
    }

    double SpectralAngle(const std::vector<double>& x, const std::vector<double>& y)
    {
      OPENSWATH_PRECONDITION(x.size() == y.size() && !x.empty(), "Both vectors need to have the same non-zero length");
      Eigen::Map<const Eigen::VectorXd> vx(x.data(), x.size());
      Eigen::Map<const Eigen::VectorXd> vy(y.data(), y.size());
      double x_len = vx.norm();
      double y_len = vy.norm();
      double denominator = x_len * y_len;
      double theta = (denominator == 0) ? 0.0 : vx.dot(vy) / denominator;
      theta = std::max(-1.0, std::min(1.0, theta));
      return std::acos(theta);
    }

    XCorrArrayType::const_iterator xcorrArrayGetMaxPeak(const XCorrArrayType& array)
    {
      OPENSWATH_PRECONDITION(array.data.size() > 0, "Cannot get highest apex from empty array.");

      XCorrArrayType::const_iterator max_it = array.begin();
      double max = array.begin()->second;
      for (XCorrArrayType::const_iterator it = array.begin(); it != array.end(); ++it)
      {
        if (it->second > max)
        {
          max = it->second;
          max_it = it;
        }
      }
      return max_it;
    }

    void standardize_data(std::vector<double>& data)
    {
      if (data.empty())
      {
        return;
      }

      Eigen::Map<Eigen::VectorXd> v(data.data(), data.size());
      double mean = v.mean();
      double stdev = std::sqrt((v.array() - mean).square().mean());

      if (mean == 0 && stdev == 0)
      {
        return; // all data is zero
      }
      if (stdev == 0)
      {
        stdev = 1; // all data is equal
      }
      v.array() = (v.array() - mean) / stdev;
    }

    XCorrArrayType normalizedCrossCorrelation(std::vector<double>& data1,
                                              std::vector<double>& data2, const int maxdelay, const int lag = 1)
    {
      OPENSWATH_PRECONDITION(data1.size() != 0 && data1.size() == data2.size(), "Both data vectors need to have the same length");

      // normalize the data
      standardize_data(data1);
      standardize_data(data2);
      return normalizedCrossCorrelationPost(data1, data2, maxdelay, lag);
    }

    XCorrArrayType normalizedCrossCorrelationPost(std::vector<double>& normalized_data1,
                                                  std::vector<double>& normalized_data2, const int maxdelay, const int lag = 1)
    {
      XCorrArrayType result = calculateCrossCorrelation(normalized_data1, normalized_data2, maxdelay, lag);

      for (XCorrArrayType::iterator it = result.begin(); it != result.end(); ++it)
      {
        it->second /= normalized_data1.size();
      }
      return result;
    }

    XCorrArrayType calculateCrossCorrelation(const std::vector<double>& data1,
                                             const std::vector<double>& data2, const int maxdelay, const int lag)
    {
      OPENSWATH_PRECONDITION(data1.size() == data2.size(), "Both data vectors need to have the same length");

      XCorrArrayType result;
      result.data.reserve( (size_t)std::ceil((2*maxdelay + 1) / lag));
      const int datasize = static_cast<int>(data1.size());
      const double* __restrict d1 = data1.data();
      const double* __restrict d2 = data2.data();

      for (int delay = -maxdelay; delay <= maxdelay; delay += lag)
      {
        const int start = std::max(0, -delay);
        const int end = std::min(datasize, datasize - delay);
        const int len = end - start;

        double sxy = 0.0;
        if (len > 0)
        {
          Eigen::Map<const Eigen::VectorXd> sub1(d1 + start, len);
          Eigen::Map<const Eigen::VectorXd> sub2(d2 + start + delay, len);
          sxy = sub1.dot(sub2);
        }
        result.data.push_back(std::make_pair(delay, sxy));
      }
      return result;
    }

    XCorrArrayType calcxcorr_legacy_mquest_(std::vector<double>& data1,
                                            std::vector<double>& data2, bool normalize)
    {
      OPENSWATH_PRECONDITION(!data1.empty() && data1.size() == data2.size(), "Both data vectors need to have the same length");
      int datasize = static_cast<int>(data1.size());
      int maxdelay = datasize;
      int lag = 1;

      Eigen::Map<Eigen::VectorXd> map1(data1.data(), datasize);
      Eigen::Map<Eigen::VectorXd> map2(data2.data(), datasize);

      double mean1 = map1.mean();
      double mean2 = map2.mean();
      double denominator = 1.0;

      if (normalize)
      {
        double sqsum1 = (map1.array() - mean1).square().sum();
        double sqsum2 = (map2.array() - mean2).square().sum();
        denominator = std::sqrt(sqsum1 * sqsum2);
      }

      denominator = (denominator > 0) ? (1.0 / denominator) : 0.0;

      XCorrArrayType result;
      result.data.reserve( (size_t)std::ceil((2*maxdelay + 1.0) / lag));

      for (int delay = -maxdelay; delay <= maxdelay; delay += lag)
      {
        double sxy = 0.0;
        int start = std::max(0, -delay);
        int end = std::min(datasize, datasize - delay);
        int len = end - start;

        if (len > 0)
        {
          Eigen::Map<Eigen::VectorXd> sub1(data1.data() + start, len);
          Eigen::Map<Eigen::VectorXd> sub2(data2.data() + start + delay, len);

          if (normalize)
          {
            sxy = (sub1.array() - mean1).matrix().dot((sub2.array() - mean2).matrix());
          }
          else
          {
            sxy = sub1.dot(sub2);
          }
        }

        if (denominator > 0)
        {
          result.data.emplace_back(delay, sxy * denominator);
        }
        else
        {
          result.data.emplace_back(delay, 0);
        }
      }
      return result;
    }

    unsigned int computeAndAppendRank(const std::vector<double>& v_temp, std::vector<unsigned int>& ranks_out)
    {
      std::vector<unsigned int> ranks{};
      ranks.resize(v_temp.size());
      std::iota(ranks.begin(), ranks.end(), 0);
      std::sort(ranks.begin(), ranks.end(),
                [&v_temp](unsigned int i, unsigned int j) { return v_temp[i] < v_temp[j]; });
      ranks_out.resize(v_temp.size());
      double x = 0;
      unsigned int y = 0;
      for(unsigned int i = 0; i < ranks.size();++i)
      {
        if(v_temp[ranks[i]] != x)
        {
          x = v_temp[ranks[i]];
          y = i;
        }
        ranks_out[ranks[i]] = y;
      }
      return y;
    }

    std::vector<unsigned int> computeRankVector(const std::vector<std::vector<double>>& intensity, std::vector<std::vector<unsigned int>>& ranks)
    {
      unsigned int pre_rank_size = ranks.size();
      ranks.resize(pre_rank_size + intensity.size());
      std::vector<unsigned int> max_rank_vec(intensity.size());
      for (std::size_t i = 0; i < intensity.size(); i++)
      {
        max_rank_vec[i] = computeAndAppendRank(intensity[i], ranks[pre_rank_size + i]);
      }
      return max_rank_vec;
    }

    double rankedMutualInformation(std::vector<unsigned int>& ranked_data1, std::vector<unsigned int>& ranked_data2, const unsigned int max_rank1, const unsigned int max_rank2)
    {
      OPENSWATH_PRECONDITION(ranked_data1.size() != 0 && ranked_data1.size() == ranked_data2.size(), "Both data vectors need to have the same length");

      unsigned int inputVectorlength = ranked_data1.size();
      unsigned int firstNumStates = max_rank1 + 1;
      unsigned int secondNumStates = max_rank2 + 1;
      std::vector<double> firstStateCounts(firstNumStates,0);
      std::vector<double> secondStateCounts(secondNumStates,0);
      std::unordered_map<pos2D, double, pair_hash> jointStateCounts{};

      for (unsigned int i = 0; i < inputVectorlength; i++) {
        firstStateCounts[ranked_data1[i]] += 1;
        secondStateCounts[ranked_data2[i]] += 1;
        jointStateCounts[std::make_pair(ranked_data1[i], ranked_data2[i])] += 1;
      }

      double mutualInformation = 0.0;
      for (const auto &[pos, jointStateCount_val]: jointStateCounts) {
        mutualInformation += jointStateCount_val * log(jointStateCount_val / firstStateCounts[pos.first] / secondStateCounts[pos.second]);
      }

      mutualInformation /= inputVectorlength;
      mutualInformation += log(inputVectorlength);
      mutualInformation /= log(2.0);

      return mutualInformation;
    }
}      //namespace OpenMS  // namespace Scoring
