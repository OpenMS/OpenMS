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

namespace OpenSwath::Scoring
{
    void normalize_sum(double x[], unsigned int n)
    { 
      double sumx = std::accumulate(&x[0], &x[0] + n, 0.0);
      if (sumx == 0.0)
      { // avoid divide by zero below
        return;
      }                           
      auto inverse_sum = 1 / sumx; // precompute inverse since division is expensive!
      for (unsigned int i = 0; i < n; ++i)
      {
        x[i] *= inverse_sum;
      }
    }

    double NormalizedManhattanDist(double x[], double y[], int n)
    {
      OPENSWATH_PRECONDITION(n > 0, "Need at least one element");

      double delta_ratio_sum = 0;
      normalize_sum(x, n);
      normalize_sum(y, n);
      for (int i = 0; i < n; i++)
      {
        delta_ratio_sum += std::fabs(x[i] - y[i]);
      }
      return delta_ratio_sum / n;
    }

    double RootMeanSquareDeviation(double x[], double y[], int n)
    {
      OPENSWATH_PRECONDITION(n > 0, "Need at least one element");

      double result = 0;
      for (int i = 0; i < n; i++)
      {

        result += (x[i] - y[i]) * (x[i] - y[i]);
      }
      return std::sqrt(result / n);
    }

    double SpectralAngle(double x[], double y[], int n)
    {
      OPENSWATH_PRECONDITION(n > 0, "Need at least one element");

      double dotprod = 0;
      double x_len = 0;
      double y_len = 0;
      for (int i = 0; i < n; i++)
      {
        dotprod += x[i] * y[i];
        x_len += x[i] * x[i];
        y_len += y[i] * y[i];
      }
      x_len = std::sqrt(x_len);
      y_len = std::sqrt(y_len);

      // normalise, avoiding a divide by zero. See unit tests for what happens
      // when one of the vectors has a length of zero.
      double denominator = x_len * y_len;
      double theta = (denominator == 0) ? 0.0 : dotprod / denominator;

      // clip to range [-1, 1] to save acos blowing up
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

      // subtract the mean and divide by the standard deviation
      double mean = std::accumulate(data.begin(), data.end(), 0.0) / (double) data.size();
      double sqsum = 0;
      for (std::vector<double>::iterator it = data.begin(); it != data.end(); ++it)
      {
        sqsum += (*it - mean) * (*it - mean);
      }
      double stdev = sqrt(sqsum / data.size()); // standard deviation

      if (mean == 0 && stdev == 0)
      {
        return; // all data is zero
      }
      if (stdev == 0)
      {
        stdev = 1; // all data is equal
      }
      for (std::size_t i = 0; i < data.size(); i++)
      {
        data[i] = (data[i] - mean) / stdev;
      }
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
      int datasize = static_cast<int>(data1.size());
      int i, j, delay;

      for (delay = -maxdelay; delay <= maxdelay; delay = delay + lag)
      {
        double sxy = 0;
        for (i = 0; i < datasize; ++i)
        {
          j = i + delay;
          if (j < 0 || j >= datasize)
          {
            continue;
          }
          sxy += (data1[i]) * (data2[j]);
        }
        result.data.push_back(std::make_pair(delay, sxy));
      }
      return result;
    }

    XCorrArrayType calcxcorr_legacy_mquest_(std::vector<double>& data1,
                                            std::vector<double>& data2, bool normalize)
    {
      OPENSWATH_PRECONDITION(!data1.empty() && data1.size() == data2.size(), "Both data vectors need to have the same length");
      int maxdelay = static_cast<int>(data1.size());
      int lag = 1;

      double mean1 = std::accumulate(data1.begin(), data1.end(), 0.) / (double)data1.size();
      double mean2 = std::accumulate(data2.begin(), data2.end(), 0.) / (double)data2.size();
      double denominator = 1.0;
      int datasize = static_cast<int>(data1.size());
      int i, j, delay;

      // Normalized cross-correlation = subtract the mean and divide by the standard deviation
      if (normalize)
      {
        double sqsum1 = 0;
        double sqsum2 = 0;
        for (std::vector<double>::iterator it = data1.begin(); it != data1.end(); ++it)
        {
          sqsum1 += (*it - mean1) * (*it - mean1);
        }

        for (std::vector<double>::iterator it = data2.begin(); it != data2.end(); ++it)
        {
          sqsum2 += (*it - mean2) * (*it - mean2);
        }
        // sigma_1 * sigma_2 * n
        denominator = sqrt(sqsum1 * sqsum2);
      }
      //avoids division in the for loop
      denominator = 1/denominator;
      XCorrArrayType result;
      result.data.reserve( (size_t)std::ceil((2*maxdelay + 1) / lag));
      int cnt = 0;
      for (delay = -maxdelay; delay <= maxdelay; delay = delay + lag, cnt++)
      {
        double sxy = 0;
        for (i = 0; i < datasize; i++)
        {
          j = i + delay;
          if (j < 0 || j >= datasize)
          {
            continue;
          }
          if (normalize)
          {
            sxy += (data1[i] - mean1) * (data2[j] - mean2);
          }
          else
          {
            sxy += (data1[i]) * (data2[j]);
          }
        }

        if (denominator > 0)
        {
          result.data.emplace_back(delay, sxy*denominator);
        }
        else
        {
          // e.g. if all datapoints are zero
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
