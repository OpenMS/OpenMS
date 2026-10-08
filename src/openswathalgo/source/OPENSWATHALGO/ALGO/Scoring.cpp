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
#include <cstdint>
#include <Eigen/Core>
#include <Eigen/Dense>

namespace OpenSwath::Scoring
{
    void normalize_sum(std::vector<double>& x)
    {
      Eigen::Map<Eigen::VectorXd> v(x.data(), x.size());
      double s = v.sum();
      if (s != 0.0) { v /= s; }
    }

    double NormalizedManhattanDist(std::vector<double>& x, std::vector<double>& y)
    {
      OPENSWATH_PRECONDITION(x.size() == y.size() && !x.empty(), "Both vectors need to have the same non-zero length");
      normalize_sum(x);
      normalize_sum(y);
      Eigen::Map<Eigen::VectorXd> vx(x.data(), x.size());
      Eigen::Map<Eigen::VectorXd> vy(y.data(), y.size());
      return (vx - vy).cwiseAbs().sum() / x.size();
    }

    double RootMeanSquareDeviation(const std::vector<double>& x, const std::vector<double>& y)
    {
      OPENSWATH_PRECONDITION(x.size() == y.size() && !x.empty(), "Both vectors need to have the same non-zero length");
      Eigen::Map<const Eigen::VectorXd> vx(x.data(), x.size());
      Eigen::Map<const Eigen::VectorXd> vy(y.data(), y.size());
      return std::sqrt((vx - vy).squaredNorm() / x.size());
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

      const unsigned int inputVectorlength = ranked_data1.size();
      const unsigned int firstNumStates = max_rank1 + 1;
      const unsigned int secondNumStates = max_rank2 + 1;
      std::vector<unsigned int> firstStateCounts(firstNumStates, 0);
      std::vector<unsigned int> secondStateCounts(secondNumStates, 0);
      std::vector<std::uint64_t> jointStates;
      jointStates.reserve(inputVectorlength);

      for (unsigned int i = 0; i < inputVectorlength; i++)
      {
        const unsigned int first_rank = ranked_data1[i];
        const unsigned int second_rank = ranked_data2[i];
        ++firstStateCounts[first_rank];
        ++secondStateCounts[second_rank];
        jointStates.push_back(static_cast<std::uint64_t>(first_rank) * secondNumStates + second_rank);
      }

      std::sort(jointStates.begin(), jointStates.end());

      double mutualInformation = 0.0;
      for (std::size_t i = 0; i < jointStates.size();)
      {
        const std::uint64_t joint_state = jointStates[i];
        std::size_t j = i + 1;
        while (j < jointStates.size() && jointStates[j] == joint_state)
        {
          ++j;
        }

        const double jointStateCount_val = static_cast<double>(j - i);
        const unsigned int first_rank = static_cast<unsigned int>(joint_state / secondNumStates);
        const unsigned int second_rank = static_cast<unsigned int>(joint_state % secondNumStates);
        mutualInformation += jointStateCount_val *
          log(jointStateCount_val / firstStateCounts[first_rank] / secondStateCounts[second_rank]);
        i = j;
      }

      mutualInformation /= inputVectorlength;
      mutualInformation += log(inputVectorlength);
      mutualInformation /= log(2.0);

      return mutualInformation;
    }
}      //namespace OpenMS  // namespace Scoring
