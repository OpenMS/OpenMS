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
// $Maintainer: Clemens Groepl $
// $Authors: Clemens Groepl, Johannes Junker, Mathias Walzer$
// --------------------------------------------------------------------------

#include <numeric>
#include <algorithm>
#include <OpenMS/CONCEPT/Types.h>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/casts.hpp>
#include <boost/function/function_base.hpp>

#ifndef OPENMS_MATH_STATISTICS_STATISTICFUNCTIONS_H
#define OPENMS_MATH_STATISTICS_STATISTICFUNCTIONS_H

namespace OpenMS
{

  namespace Math
  {
    /**
      @brief Calculates the sum of a range of values

      @ingroup MathFunctionsStatistics
    */
    template <typename IteratorType>
    static DoubleReal sum(IteratorType begin, IteratorType end)
    {
      return std::accumulate(begin, end, 0.0);
    }

    /**
      @brief Calculates the mean of a range of values

      @exception Exception::InvalidRange is thrown if the range is empty

      @ingroup MathFunctionsStatistics
    */
    template <typename IteratorType>
    static DoubleReal mean(IteratorType begin, IteratorType end)
    {
      SignedSize size = std::distance(begin, end);
      if (size <= 0)
      {
        throw Exception::InvalidRange(__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }
      return sum(begin, end) / size;
    }

    /**
      @brief Calculates the median of a range of values

      @param begin Start of range
      @param end End of range (past-the-end iterator)
      @param sorted Is the range already sorted? If not, it will be sorted.

      @exception Exception::InvalidRange is thrown if the range is empty

      @ingroup MathFunctionsStatistics
    */
    template <typename IteratorType>
    static DoubleReal median(IteratorType begin, IteratorType end, bool sorted = FALSE)
    {
      Size size = std::distance(begin, end);

      if (size == 0)
      {
        throw Exception::InvalidRange(__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }

      if (!sorted)
      {
        std::sort(begin, end);
      }

      if (size % 2 == 0)        // even size => average two middle values
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
      @brief Calculates the quantile of a range of values

      @param begin Start of range
      @param end End of range (past-the-end iterator)
      @param sorted Is the range already sorted? If not, it will be sorted.

      @exception Exception::InvalidRange is thrown if the range is empty or a quantile over 100 is given

      @ingroup MathFunctionsStatistics
    */
    template <typename IteratorType>
    static DoubleReal quantile(IteratorType begin, IteratorType end, UInt quantile, bool sorted = FALSE)
    {
      Size size = std::distance(begin, end);

      if (size == 0)
      {
        throw Exception::InvalidRange(__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }
      if (quantile > 100 || quantile < 1) //TODO is 0 quantile a valid request?
      {
        throw Exception::InvalidRange(__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }

      int l = floor( (double(quantile) * (double(size) / 100)) + 0.5); // will not be negative, so this is round nearest

      if (!sorted)
      {
        std::sort(begin, end);
      }

      IteratorType it = begin;
      std::advance(it, l - 1);
      return *it;

    }

    /**
      @brief Calculates the mean square error for the values in [begin_a, end_a) and [begin_b, end_b)

      Calculates the mean square error for the data given by the two iterator ranges.

      @exception Exception::InvalidRange is thrown if the iterator ranges are not of the same length or empty.

      @ingroup MathFunctionsStatistics
    */
    template <typename IteratorType1, typename IteratorType2>
    static DoubleReal meanSquareError(IteratorType1 begin_a, IteratorType1 end_a, IteratorType2 begin_b, IteratorType2 end_b)
    {
      //no data or different lengths
      SignedSize dist = std::distance(begin_a, end_a);
      if (dist == 0 || dist != std::distance(begin_b, end_b))
      {
        throw Exception::InvalidRange(__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }

      DoubleReal error = 0;
      while (begin_a != end_a)
      {
        DoubleReal tmp(*begin_a - *begin_b);
        error += tmp * tmp;
        ++begin_a;
        ++begin_b;
      }

      return error / dist;
    }

    /**
      @brief Calculates the classification rate for the values in [begin_a,	end_a) and [begin_b, end_b)

      Calculates the classification rate for the data given by the two iterator ranges.

      @exception Exception::InvalidRange is thrown if the iterator ranges are not of the same length or empty.

      @ingroup MathFunctionsStatistics
    */
    template <typename IteratorType1, typename IteratorType2>
    static DoubleReal classificationRate(IteratorType1 begin_a, IteratorType1 end_a, IteratorType2 begin_b, IteratorType2 end_b)
    {
      //no data or different lengths
      SignedSize dist = std::distance(begin_a, end_a);
      if (dist == 0 || dist != std::distance(begin_b, end_b))
      {
        throw Exception::InvalidRange(__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }

      DoubleReal correct = (DoubleReal) dist;
      while (begin_a != end_a)
      {
        if ((*begin_a < 0 && *begin_b >= 0) || (*begin_a >= 0 && *begin_b < 0))
        {
          --correct;
        }
        ++begin_a;
        ++begin_b;
      }

      return correct / dist;
    }

    /**
      @brief Calculates the Matthews correlation coefficient for the values in [begin_a, end_a) and [begin_b, end_b)

      Calculates the Matthews correlation coefficient for the data given by the two iterator ranges. The values in [begin_a, end_a) have to be the predicted labels and the values in [begin_b, end_b) have to be the real labels.

      @exception Exception::InvalidRange is thrown if the iterator ranges are not of the same length or empty.

      @ingroup MathFunctionsStatistics
    */
    template <typename IteratorType1, typename IteratorType2>
    static DoubleReal matthewsCorrelationCoefficient(IteratorType1 begin_a, IteratorType1 end_a, IteratorType2 begin_b, IteratorType2 end_b)
    {
      //no data or different lengths
      Int dist = std::distance(begin_a, end_a);
      if (dist == 0 || dist != std::distance(begin_b, end_b))
      {
        throw Exception::InvalidRange(__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }

      DoubleReal tp = 0;
      DoubleReal fp = 0;
      DoubleReal tn = 0;
      DoubleReal fn = 0;

      while (begin_a != end_a)
      {
        if (*begin_a < 0 && *begin_b >= 0)
        {
          ++fn;
        }
        else if (*begin_a < 0 && *begin_b < 0)
        {
          ++tn;
        }
        else if (*begin_a >= 0 && *begin_b >= 0)
        {
          ++tp;
        }
        else if (*begin_a >= 0 && *begin_b < 0)
        {
          ++fp;
        }

        ++begin_a;
        ++begin_b;
      }

      return (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));
    }

    /**
      @brief Calculates the Pearson correlation coefficient for the values in [begin_a, end_a) and [begin_b, end_b)

      Calculates the linear correlation coefficient for the data given by the two iterator ranges.

      If one of the ranges contains only the same values 'nan' is returned.

      @exception Exception::InvalidRange is thrown if the iterator ranges are not of the same length or empty.

      @ingroup MathFunctionsStatistics
    */
    template <typename IteratorType1, typename IteratorType2>
    static DoubleReal pearsonCorrelationCoefficient(IteratorType1 begin_a, IteratorType1 end_a, IteratorType2 begin_b, IteratorType2 end_b)
    {
      //no data or different lengths
      SignedSize dist = std::distance(begin_a, end_a);
      if (dist == 0 || dist != std::distance(begin_b, end_b))
      {
        throw Exception::InvalidRange(__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }

      //calculate average
      DoubleReal avg_a = std::accumulate(begin_a, end_a, 0.0) / dist;
      DoubleReal avg_b = std::accumulate(begin_b, end_b, 0.0) / dist;

      DoubleReal numerator = 0;
      DoubleReal denominator_a = 0;
      DoubleReal denominator_b = 0;
      while (begin_a != end_a)
      {
        DoubleReal temp_a = *begin_a - avg_a;
        DoubleReal temp_b = *begin_b - avg_b;
        numerator += (temp_a * temp_b);
        denominator_a += (temp_a * temp_a);
        denominator_b += (temp_b * temp_b);
        ++begin_a;
        ++begin_b;
      }

      return numerator / sqrt(denominator_a * denominator_b);
    }

    /// Replaces the elements in vector @p w by their ranks
    template <typename Value>
    static void computeRank(std::vector<Value> & w)
    {
      using namespace boost::lambda;
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
                ret<bool>((&_1->*& std::pair<Size, Value>::second) < (&_2->*& std::pair<Size, Value>::second)));
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
    static DoubleReal rankCorrelationCoefficient(IteratorType1 begin_a, IteratorType1 end_a, IteratorType2 begin_b, IteratorType2 end_b)
    {
      //no data or different lengths
      SignedSize dist = std::distance(begin_a, end_a);
      if (dist == 0 || dist != std::distance(begin_b, end_b))
      {
        throw Exception::InvalidRange(__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }

      // store and sort intensities of model and data
      std::vector<DoubleReal> ranks_data;
      ranks_data.reserve(dist);
      std::vector<DoubleReal> ranks_model;
      ranks_model.reserve(dist);

      while (begin_a != end_a)
      {
        ranks_model.push_back(*begin_a);
        ranks_data.push_back(*begin_b);
        ++begin_a;
        ++begin_b;
      }

      // replace entries by their ranks
      computeRank(ranks_data);
      computeRank(ranks_model);

      DoubleReal mu = DoubleReal(ranks_data.size() + 1) / 2.; // mean of ranks
      // Was the following, but I think the above is more correct ... (Clemens)
      // DoubleReal mu = (ranks_data.size() + 1) / 2;

      DoubleReal sum_model_data = 0;
      DoubleReal sqsum_data = 0;
      DoubleReal sqsum_model = 0;

      for (Int i = 0; i < dist; ++i)
      {
        sum_model_data += (ranks_data[i] - mu) * (ranks_model[i] - mu);
        sqsum_data += (ranks_data[i] - mu) * (ranks_data[i] - mu);
        sqsum_model += (ranks_model[i] - mu) * (ranks_model[i] - mu);
      }

      // check for division by zero
      if (!sqsum_data || !sqsum_model)
        return 0;

      return sum_model_data / (sqrt(sqsum_data) * sqrt(sqsum_model));
    }

  }   // namespace Math
} // namespace OpenMS

#endif // OPENMS_MATH_STATISTICS_STATISTICFUNCTIONS_H
