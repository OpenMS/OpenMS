// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Clemens Groepl, Johannes Junker, Mathias Walzer, Chris Bielow $
// --------------------------------------------------------------------------
#pragma once

#include <vector>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Types.h>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <numeric>

namespace OpenMS
{

  namespace Math
  {

    /**
      @brief Helper function checking if two iterators are not equal

      @exception Exception::InvalidRange is thrown if the range is NULL

      @ingroup MathFunctionsStatistics
    */
    template <typename IteratorType>
    static void checkIteratorsNotNULL(IteratorType begin, IteratorType end)
    {
      if (begin == end)
      {
        throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
      }
    }

    /**
       @brief Helper function checking if two iterators are equal
       
       @exception Exception::InvalidRange is thrown if the iterators are not equal

       @ingroup MathFunctionsStatistics
    */
    template <typename IteratorType>
    static void checkIteratorsEqual(IteratorType begin, IteratorType end)
    {
      if (begin != end)
      {
        throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
      }
    }

    /**
       @brief Helper function checking if an iterator and a co-iterator both have a next element
       
       @exception Exception::InvalidRange is thrown if the iterator do not end simultaneously
       
       @ingroup MathFunctionsStatistics
    */
    template <typename IteratorType1, typename IteratorType2>
    static void checkIteratorsAreValid(
      IteratorType1 begin_b, IteratorType1 end_b,
      IteratorType2 begin_a, IteratorType2 end_a)
    {
      if (begin_b != end_b && begin_a == end_a)
      {
        throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
      }
    }

    /**
       @brief Calculates the sum of a range of values
       
       @ingroup MathFunctionsStatistics
    */
    template <typename IteratorType>
    static double sum(IteratorType begin, IteratorType end)
    {
      return std::accumulate(begin, end, 0.0);
    }

    /**
       @brief Calculates the mean of a range of values
      
       @exception Exception::InvalidRange is thrown if the range is NULL
       
       @ingroup MathFunctionsStatistics
    */
    template <typename IteratorType>
    static double mean(IteratorType begin, IteratorType end)
    {
      checkIteratorsNotNULL(begin, end);
      return sum(begin, end) / std::distance(begin, end);
    }

    /**
       @brief Calculates the median of a range of values
       
       @param begin Start of range
       @param end End of range (past-the-end iterator)
       @param sorted Is the range already sorted? If not, it will be sorted.
    @return Median (as floating point, since we need to support average of middle values)
       @exception Exception::InvalidRange is thrown if the range is NULL

       @ingroup MathFunctionsStatistics
    */
    template <typename IteratorType>
    static double median(IteratorType begin, IteratorType end, 
                         bool sorted = false)
    {
      checkIteratorsNotNULL(begin, end);
      if (!sorted)
      {
        std::sort(begin, end);
      }
      
      Size size = std::distance(begin, end);
      if (size % 2 == 0) // even size => average two middle values
      {
        IteratorType it1 = begin;
        std::advance(it1, size / 2 - 1);
        IteratorType it2 = it1;
        std::advance(it2, 1);
        return (*it1 + *it2) / 2.0;
      }
      else
      {
        IteratorType it = begin;
        std::advance(it, (size - 1) / 2);
        return *it;
      }
    }

  
    /** 
      @brief median absolute deviation (MAD)

      Computes the MAD, defined as 

      MAD = median( | x_i - median(x) | ) for a vector x with indices i in [1,n].

      Sortedness of the input is not required (nor does it provide a speedup).
      For efficiency, you must provide the median separately, in order to avoid potentially duplicate efforts (usually one
      computes the median anyway externally).
      
      @param begin Start of range
      @param end End of range (past-the-end iterator)
      @param median_of_numbers The precomputed median of range @p begin - @p end.
      @return the MAD

      @ingroup MathFunctionsStatistics

    */
    template <typename IteratorType>
    double MAD(IteratorType begin, IteratorType end, double median_of_numbers)
    {
      std::vector<double> diffs;
      diffs.reserve(std::distance(begin, end));
      for (IteratorType it = begin; it != end; ++it)
      {
        diffs.push_back(fabs(*it - median_of_numbers));
      }
      return median(diffs.begin(), diffs.end(), false);
    }
    
    /**
      @brief mean absolute deviation (MeanAbsoluteDeviation)

      Computes the MeanAbsoluteDeviation, defined as

      MeanAbsoluteDeviation = mean( | x_i - mean(x) | ) for a vector x with indices i in [1,n].

      For efficiency, you must provide the mean separately, in order to avoid potentially duplicate efforts (usually one
      computes the mean anyway externally).
      
      @param begin Start of range
      @param end End of range (past-the-end iterator)
      @param mean_of_numbers The precomputed mean of range @p begin - @p end.
      @return the MeanAbsoluteDeviation

      @ingroup MathFunctionsStatistics

    */
    template <typename IteratorType>
    double MeanAbsoluteDeviation(IteratorType begin, IteratorType end, double mean_of_numbers)
    {
      double mean_value {0};
      for (IteratorType it = begin; it != end; ++it)
      {
        mean_value += fabs(*it - mean_of_numbers);
      }
      return mean_value / std::distance(begin, end);
    }

    /**
       @brief Calculates the first quantile of a range of values
       
       The range is divided into half and the median for the first half is returned.

       @param begin Start of range
       @param end End of range (past-the-end iterator)
       @param sorted Is the range already sorted? If not, it will be sorted.
       
       @exception Exception::InvalidRange is thrown if the range is NULL
       
       @ingroup MathFunctionsStatistics
    */
    template <typename IteratorType>
    static double quantile1st(IteratorType begin, IteratorType end, 
                              bool sorted = false)
    {
      checkIteratorsNotNULL(begin, end);

      if (!sorted)
      {
        std::sort(begin, end);
      }

      Size size = std::distance(begin, end);
      if (size % 2 == 0)
      {
        return median(begin, begin + (size/2)-1, true); //-1 to exclude median values
      }
      return median(begin, begin + (size/2), true);
    }

    /**
       @brief Calculates the third quantile of a range of values

       The range is divided into half and the median for the second half is returned.

       @param begin Start of range
       @param end End of range (past-the-end iterator)
       @param sorted Is the range already sorted? If not, it will be sorted.

       @exception Exception::InvalidRange is thrown if the range is NULL

       @ingroup MathFunctionsStatistics
    */
    template <typename IteratorType>
    static double quantile3rd(
      IteratorType begin, IteratorType end, bool sorted = false)
    {
      checkIteratorsNotNULL(begin, end);
      if (!sorted)
      {
        std::sort(begin, end);
      }

      Size size = std::distance(begin, end);
      return median(begin + (size/2)+1, end, true); //+1 to exclude median values
    }

    /**
       @brief Calculates the variance of a range of values

  The @p mean can be provided explicitly to save computation time. If left at default, it will be computed internally.

       @exception Exception::InvalidRange is thrown if the range is empty

       @ingroup MathFunctionsStatistics
    */
    template <typename IteratorType>
    static double variance(IteratorType begin, IteratorType end,
                           double mean = std::numeric_limits<double>::max())
    {
      checkIteratorsNotNULL(begin, end);
      double sum_value = 0.0;
      if (mean == std::numeric_limits<double>::max())
      {
        mean = Math::mean(begin, end);
      }
      for (IteratorType iter=begin; iter!=end; ++iter)
      {
        double diff = *iter - mean;
        sum_value += diff * diff;
      }
      return sum_value / (std::distance(begin, end)-1);
    }

    /**
       @brief Calculates the standard deviation of a range of values.

       The @p mean can be provided explicitly to save computation time. If left at default, it will be computed internally.

       @exception Exception::InvalidRange is thrown if the range is empty

       @ingroup MathFunctionsStatistics
    */
    template <typename IteratorType>
    static double sd(IteratorType begin, IteratorType end,
                     double mean = std::numeric_limits<double>::max())
    {
      checkIteratorsNotNULL(begin, end);
      return std::sqrt( variance(begin, end, mean) );
    }

    /**
       @brief Calculates the absolute deviation of a range of values

       @exception Exception::InvalidRange is thrown if the range is empty

       @ingroup MathFunctionsStatistics
    */
    template <typename IteratorType>
    static double absdev(IteratorType begin, IteratorType end,
                         double mean = std::numeric_limits<double>::max())
    {
      checkIteratorsNotNULL(begin, end);
      double sum_value = 0.0;
      if (mean == std::numeric_limits<double>::max())
      {
        mean = Math::mean(begin, end);
      }
      for (IteratorType iter=begin; iter!=end; ++iter)
      {
        sum_value += *iter - mean;
      }
      return sum_value / std::distance(begin, end);
    }

    /**
       @brief Calculates the covariance of two ranges of values.

       Note that the two ranges must be of equal size.

       @exception Exception::InvalidRange is thrown if the range is empty

       @ingroup MathFunctionsStatistics
    */
    template <typename IteratorType1, typename IteratorType2>
    static double covariance(IteratorType1 begin_a, IteratorType1 end_a,
                             IteratorType2 begin_b, IteratorType2 end_b)
    {
      //no data or different lengths
      checkIteratorsNotNULL(begin_a, end_a);

      double sum_value = 0.0;
      double mean_a = Math::mean(begin_a, end_a);
      double mean_b = Math::mean(begin_b, end_b);
      IteratorType1 iter_a = begin_a;
      IteratorType2 iter_b = begin_b;
      for (; iter_a != end_a; ++iter_a, ++iter_b)
      {
        /* assure both ranges have the same number of elements */
        checkIteratorsAreValid(begin_b, end_b, begin_a, end_a);
        sum_value += (*iter_a - mean_a) * (*iter_b - mean_b);
      }
      /* assure both ranges have the same number of elements */
      checkIteratorsEqual(iter_b, end_b);
      Size n = std::distance(begin_a, end_a);
      return sum_value / (n-1);
    }

    /**
       @brief Calculates the mean square error for the values in [begin_a, end_a) and [begin_b, end_b)

       Calculates the mean square error for the data given by the two iterator ranges.

       @exception Exception::InvalidRange is thrown if the iterator ranges are not of the same length or empty.

       @ingroup MathFunctionsStatistics
    */
    template <typename IteratorType1, typename IteratorType2>
    static double meanSquareError(IteratorType1 begin_a, IteratorType1 end_a,
                                  IteratorType2 begin_b, IteratorType2 end_b)
    {
      //no data or different lengths
      checkIteratorsNotNULL(begin_a, end_a);

      SignedSize dist = std::distance(begin_a, end_a);
      double error = 0;
      IteratorType1 iter_a = begin_a;
      IteratorType2 iter_b = begin_b;
      for (; iter_a != end_a; ++iter_a, ++iter_b)
      {
        /* assure both ranges have the same number of elements */
        checkIteratorsAreValid(iter_b, end_b, iter_a, end_a);

        double tmp(*iter_a - *iter_b);
        error += tmp * tmp;
      }
      /* assure both ranges have the same number of elements */
      checkIteratorsEqual(iter_b, end_b);

      return error / dist;
    }

    /**
       @brief Calculates the classification rate for the values in [begin_a, end_a) and [begin_b, end_b)

       Calculates the classification rate for the data given by the two iterator ranges.

       @exception Exception::InvalidRange is thrown if the iterator ranges are not of the same length or empty.

       @ingroup MathFunctionsStatistics
    */
    template <typename IteratorType1, typename IteratorType2>
    static double classificationRate(IteratorType1 begin_a, IteratorType1 end_a,
                                     IteratorType2 begin_b, IteratorType2 end_b)
    {
      //no data or different lengths
      checkIteratorsNotNULL(begin_a, end_a);

      SignedSize dist = std::distance(begin_a, end_a);
      SignedSize correct = dist;
      IteratorType1 iter_a = begin_a;
      IteratorType2 iter_b = begin_b;
      for (; iter_a != end_a; ++iter_a, ++iter_b)
      {
        /* assure both ranges have the same number of elements */
        checkIteratorsAreValid(iter_b, end_b, iter_a, end_a);
        if ((*iter_a < 0 && *iter_b >= 0) || (*iter_a >= 0 && *iter_b < 0))
        {
          --correct;
        }

      }
      /* assure both ranges have the same number of elements */
      checkIteratorsEqual(iter_b, end_b);

      return double(correct) / dist;
    }

    /**
       @brief Calculates the Matthews correlation coefficient for the values in [begin_a, end_a) and [begin_b, end_b)

       Calculates the Matthews correlation coefficient for the data given by the
       two iterator ranges. The values in [begin_a, end_a) have to be the
       predicted labels and the values in [begin_b, end_b) have to be the real
       labels.

       @exception Exception::InvalidRange is thrown if the iterator ranges are not of the same length or empty.

       @ingroup MathFunctionsStatistics
    */
    template <typename IteratorType1, typename IteratorType2>
    static double matthewsCorrelationCoefficient(
      IteratorType1 begin_a, IteratorType1 end_a,
      IteratorType2 begin_b, IteratorType2 end_b)
    {
      //no data or different lengths
      checkIteratorsNotNULL(begin_a, end_b);

      double tp = 0;
      double fp = 0;
      double tn = 0;
      double fn = 0;
      IteratorType1 iter_a = begin_a;
      IteratorType2 iter_b = begin_b;
      for (; iter_a != end_a; ++iter_a, ++iter_b)
      {
        /* assure both ranges have the same number of elements */
        checkIteratorsAreValid(iter_b, end_b, iter_a, end_a);

        if (*iter_a < 0 && *iter_b >= 0)
        {
          ++fn;
        }
        else if (*iter_a < 0 && *iter_b < 0)
        {
          ++tn;
        }
        else if (*iter_a >= 0 && *iter_b >= 0)
        {
          ++tp;
        }
        else if (*iter_a >= 0 && *iter_b < 0)
        {
          ++fp;
        }
      }
      /* assure both ranges have the same number of elements */
      checkIteratorsEqual(iter_b, end_b);

      return (tp * tn - fp * fn) / std::sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));
    }

    /**
       @brief Calculates the Pearson correlation coefficient for the values in [begin_a, end_a) and [begin_b, end_b)

       Calculates the linear correlation coefficient for the data given by the two iterator ranges.

       If one of the ranges contains only the same values 'nan' is returned.

       @exception Exception::InvalidRange is thrown if the iterator ranges are not of the same length or empty.

       @ingroup MathFunctionsStatistics
    */
    template <typename IteratorType1, typename IteratorType2>
    static double pearsonCorrelationCoefficient(
      IteratorType1 begin_a, IteratorType1 end_a,
      IteratorType2 begin_b, IteratorType2 end_b)
    {
      //no data or different lengths
      checkIteratorsNotNULL(begin_a, end_a);

      //calculate average
      SignedSize dist = std::distance(begin_a, end_a);
      double avg_a = std::accumulate(begin_a, end_a, 0.0) / dist;
      double avg_b = std::accumulate(begin_b, end_b, 0.0) / dist;

      double numerator = 0;
      double denominator_a = 0;
      double denominator_b = 0;
      IteratorType1 iter_a = begin_a;
      IteratorType2 iter_b = begin_b;
      for (; iter_a != end_a; ++iter_a, ++iter_b)
      {
        /* assure both ranges have the same number of elements */
        checkIteratorsAreValid(iter_b, end_b, iter_a, end_a);
        double temp_a = *iter_a - avg_a;
        double temp_b = *iter_b - avg_b;
        numerator += (temp_a * temp_b);
        denominator_a += (temp_a * temp_a);
        denominator_b += (temp_b * temp_b);
      }
      /* assure both ranges have the same number of elements */
      checkIteratorsEqual(iter_b, end_b);
      return numerator / std::sqrt(denominator_a * denominator_b);
    }

    /// Replaces the elements in vector @p w by their ranks
    template <typename Value>
    static void computeRank(std::vector<Value> & w)
    {
      Size i = 0; // main index
      Size z  = 0;  // "secondary" index
      Value rank = 0;
      Size n = (w.size() - 1);
      //store original indices for later
      std::vector<std::pair<Size, Value> > w_idx;
      for (Size j = 0; j < w.size(); ++j)
      {
        w_idx.push_back(std::make_pair(j, w[j]));
      }
      //sort
      std::sort(w_idx.begin(), w_idx.end(),
                [](const auto& pair1, const auto& pair2) { return pair1.second < pair2.second; });
      //replace pairs <orig_index, value> in w_idx by pairs <orig_index, rank>
      while (i < n)
      {
        // test for equality with tolerance:
        if (fabs(w_idx[i + 1].second - w_idx[i].second) > 0.0000001 * fabs(w_idx[i + 1].second)) // no tie
        {
          w_idx[i].second = Value(i + 1);
          ++i;
        }
        else // tie, replace by mean rank
        {
          // count number of ties
          for (z = i + 1; (z <= n) && fabs(w_idx[z].second - w_idx[i].second) <= 0.0000001 * fabs(w_idx[z].second); ++z)
          {
          }
          // compute mean rank of tie
          rank = 0.5 * (i + z + 1);
          // replace intensities by rank
          for (Size v = i; v <= z - 1; ++v)
          {
            w_idx[v].second = rank;
          }
          i = z;
        }
      }
      if (i == n)
        w_idx[n].second = Value(n + 1);
      //restore original order and replace elements of w with their ranks
      for (Size j = 0; j < w.size(); ++j)
      {
        w[w_idx[j].first] = w_idx[j].second;
      }
    }

    /**
       @brief calculates the rank correlation coefficient for the values in [begin_a, end_a) and [begin_b, end_b)

       Calculates the rank correlation coefficient for the data given by the two iterator ranges.

       If one of the ranges contains only the same values 'nan' is returned.

       @exception Exception::InvalidRange is thrown if the iterator ranges are not of the same length or empty.

       @ingroup MathFunctionsStatistics
    */
    template <typename IteratorType1, typename IteratorType2>
    static double rankCorrelationCoefficient(
      IteratorType1 begin_a, IteratorType1 end_a,
      IteratorType2 begin_b, IteratorType2 end_b)
    {
      //no data or different lengths
      checkIteratorsNotNULL(begin_a, end_a);

      // store and sort intensities of model and data
      SignedSize dist = std::distance(begin_a, end_a);
      std::vector<double> ranks_data;
      ranks_data.reserve(dist);
      std::vector<double> ranks_model;
      ranks_model.reserve(dist);
      IteratorType1 iter_a = begin_a;
      IteratorType2 iter_b = begin_b;
      for (; iter_a != end_a; ++iter_a, ++iter_b)
      {
        /* assure both ranges have the same number of elements */
        checkIteratorsAreValid(iter_b, end_b, iter_a, end_a);

        ranks_model.push_back(*iter_a);
        ranks_data.push_back(*iter_b);
      }
      /* assure both ranges have the same number of elements */
      checkIteratorsEqual(iter_b, end_b);

      // replace entries by their ranks
      computeRank(ranks_data);
      computeRank(ranks_model);

      double mu = double(ranks_data.size() + 1) / 2.; // mean of ranks
      // Was the following, but I think the above is more correct ... (Clemens)
      // double mu = (ranks_data.size() + 1) / 2;

      double sum_model_data = 0;
      double sqsum_data = 0;
      double sqsum_model = 0;

      for (Int i = 0; i < dist; ++i)
      {
        sum_model_data += (ranks_data[i] - mu) * (ranks_model[i] - mu);
        sqsum_data += (ranks_data[i] - mu) * (ranks_data[i] - mu);
        sqsum_model += (ranks_model[i] - mu) * (ranks_model[i] - mu);
      }

      // check for division by zero
      if (!sqsum_data || !sqsum_model)
      {
        return 0;
      }

      return sum_model_data / (std::sqrt(sqsum_data) * std::sqrt(sqsum_model));
    }

    /// Helper class to gather (and dump) some statistics from a e.g. vector<double>.
    template<typename T>
    struct SummaryStatistics
    {
      SummaryStatistics() = default;

      // Ctor with data
      SummaryStatistics(T& data)
      {
        count = data.size();
        // Sanity check: avoid core dump if no data points present.
        if (data.empty())
        {
          mean = variance = min = lowerq = median = upperq = max = 0.0;
        }
        else
        {
          sort(data.begin(), data.end());
          mean = Math::mean(data.begin(), data.end());
          variance = Math::variance(data.begin(), data.end(), mean);
          min = data.front();
          lowerq = Math::quantile1st(data.begin(), data.end(), true);
          median = Math::median(data.begin(), data.end(), true);
          upperq = Math::quantile3rd(data.begin(), data.end(), true);
          max = data.back();
        }
      }

      double mean = 0, variance = 0 , lowerq = 0, median = 0, upperq = 0;
      typename T::value_type min = 0, max = 0;
      size_t count = 0;
    };

  }   // namespace Math
} // namespace OpenMS

