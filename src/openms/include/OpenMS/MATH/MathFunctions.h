// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/RangeManager.h>

#include <boost/random/mersenne_twister.hpp> // for mt19937_64
#include <boost/random/uniform_int.hpp>
#include <cmath>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/complement.hpp>
#include <algorithm>
#include <iterator>
#include <limits>
#include <numeric>
#include <utility> // for std::pair
#include <vector>

namespace OpenMS
{
/**
  @brief Math class.

  Contains mathematical auxiliary functions.

  @note Using a class with static methods instead of a namespace simplifies
        wrapping for other languages (e.g., Python bindings), as classes are
        more straightforward to expose than namespaces in many binding generators.

  @ingroup Concept
*/
class OPENMS_DLLAPI Math
{
public:

  /**
    @brief Given an interval/range and a new value, extend the range to include the new value if needed

    @param[in,out] min The current minimum of the range
    @param[in,out] max The current maximum of the range
    @param[in] value The new value which may extend the range
    @return true if the range was modified
  */
  template<typename T>
  static bool extendRange(T& min, T& max, const T& value)
  {
    if (value < min)
    {
      min = value;
      return true;
    }
    if (value > max)
    {
      max = value;
      return true;
    }
    return false;
  }

  /**
   * \brief Is a @p value contained in [min, max] ?
   * \tparam T Type, e.g. double
   * \return True if contained, false otherwise
   */
  template<typename T>
  static bool contains(T value, T min, T max)
  {
    return min <= value && value <= max;
  }

  /**
   * \brief Zoom into an interval [left, right], decreasing its width by @p factor (which must be in [0,inf]).
   *
   * To actually zoom in, the @p factor needs to be within [0,1]. Chosing a factor > 1 actually zooms out.
   * @p align (between [0,1]) determines where the zoom happens:
   *   i.e. align = 0 leaves @p left the same and reduces @p right (or extends if factor>1)
   *        whereas align = 0.5 zooms into the center of the range etc
   *
   * You can do round trips, i.e. undo a zoom in, by inverting the factor:
   * \code
   * [a2, b2] = zoomIn(a1, b1, 0.5, al); // zoom in
   * [a1, b1] === zoomIn(a2, b2, 2, al); // zoom out again (inverting)
   * \endcode
   *
   * \param left Start of interval
   * \param right End of interval
   * \param factor Number between [0,1] to shrink, or >1 to extend the span (=right-left)
   * \param align Where to position the smaller/shrunk interval (0 = left, 1 = right, 0.5=center etc)
   * \return [new_left, new_right] as pair
   */
  static std::pair<double, double> zoomIn(const double left, const double right, const float factor, const float align)
  {
    OPENMS_PRECONDITION(factor >= 0, "Factor must be >=0")
    OPENMS_PRECONDITION(align >= 0, "align must be >=0")
    OPENMS_PRECONDITION(align <= 1, "align must be <=1")
    std::pair<double, double> res;
    auto old_width = right - left;
    auto offset_left = (1.0f - factor) * old_width * align;
    res.first = left + offset_left;
    res.second = res.first + old_width * factor;
    return res;
  }

  using BinContainer = std::vector<RangeBase>;
  /**
    @brief Split a range [min,max] into @p number_of_bins (with optional overlap) and return the ranges of each bin.

    Optionally, bins can be made overlapping, by extending each bins' left and right margin by @p extend_margin. 
    The overlap between neighboring bins will thus be `2 x extend_margin`.
    The borders of the original interval will @em not be extended.
    
    @param[in] min The minimum of the range; must be smaller than @p max
    @param[in] max The maximum of the range
    @param[in] number_of_bins How many bins should the range be divided into? Must be 1 or larger
    @param[in] extend_margin Overlap of neighboring bins (=0 for no overlap). Negative values will shrink the range (feature).
    @return Vector with @p number_of_bins elements, each representing the margins of one bin

    @throws OpenMS::Precondition if `min >= max` or `number_of_bins == 0`
  */
  static BinContainer createBins(double min, double max, uint32_t number_of_bins, double extend_margin = 0)
  {
    OPENMS_PRECONDITION(number_of_bins >= 1, "Number of bins must be >= 1")
    OPENMS_PRECONDITION(min < max, "Require min < max");
    std::vector<RangeBase> res(number_of_bins);
    const double bin_width = (max - min) / number_of_bins;
    for (uint32_t i = 0; i < number_of_bins; ++i)
    {
      res[i] = RangeBase(min + i * bin_width, min + (i + 1) * bin_width);
      res[i].extendLeftRight(extend_margin);
    }
    res.front().setMin(min); // undo potential margin
    res.back().setMax(max);  // undo potential margin
    
    return res;
  }


  /**
    @brief rounds @p x up to the next decimal power 10 ^ @p decPow

    @verbatim
    e.g.: (123.0 , 1)  => 130
          (123.0 , 2)  => 200
              (0.123 ,-2)  => 0.13 ( 10^-2 = 0.01 )
    @endverbatim

    @ingroup MathFunctionsMisc
  */
  static double ceilDecimal(double x, int decPow)
  {
    return (ceil(x / pow(10.0, decPow))) * pow(10.0, decPow); // decimal shift right, ceiling, decimal shift left
  }

  /**
      @brief rounds @p x to the next decimal power 10 ^ @p decPow

      @verbatim
      e.g.: (123.0 , 1)  => 120
            (123.0 , 2)  => 100
      @endverbatim

      @ingroup MathFunctionsMisc
  */
  static double roundDecimal(double x, int decPow)
  {
    if (x > 0) return (floor(0.5 + x / pow(10.0, decPow))) * pow(10.0, decPow);

    return -((floor(0.5 + fabs(x) / pow(10.0, decPow))) * pow(10.0, decPow));
  }

  /**
      @brief transforms point @p x of interval [left1,right1] into interval [left2,right2]

      @ingroup MathFunctionsMisc
  */
  static double intervalTransformation(double x, double left1, double right1, double left2, double right2)
  {
    return left2 + (x - left1) * (right2 - left2) / (right1 - left1);
  }

  /**
      @brief Transforms a number from linear to log10 scale. Avoids negative logarithms by adding 1.

      @param[in] x The number to transform

      @ingroup MathFunctionsMisc
  */
  static double linear2log(double x)
  {
    return log10(x + 1); //+1 to avoid negative logarithms
  }

  /**
      @brief Transforms a number from log10 to to linear scale. Subtracts the 1 added by linear2log(double)

      @param[in] x The number to transform

      @ingroup MathFunctionsMisc
  */
  static double log2linear(double x)
  {
    return pow(10, x) - 1;
  }

  /**
      @brief Returns true if the given integer is odd

      @ingroup MathFunctionsMisc
  */
  static bool isOdd(UInt x)
  {
    return (x & 1) != 0;
  }

  /**
      @brief Rounds the value

      @ingroup MathFunctionsMisc
  */
  template<typename T>
  static T round(T x)
  {
    return std::round(x);
  }

  /**
  * 
    Rounds to the i'th digit after the decimal point, also works for negative digits.
    
    e.g.
    \code
     round_to(3.14159265,  2)  // 3.14
     round_to(1234.9    , -2)  // 1200
    \endcode

    @param[in] value The value to round
    @param[in] digits The number of digits to round to (can be negative)
    @return The rounded value
  */
  template<typename T>
  static T roundTo(const T value, int digits)
  {
    T factor = 1.0;
    if (digits > 0)
    {
      for (int i = 0; i < digits; ++i)
        factor *= 10.0;
    }
    else if (digits < 0)
    {
      for (int i = 0; i < -digits; ++i)
        factor /= 10.0;
    }

    return std::round(value * factor) / factor;
  }

  /**
    Computes the percentage of @p value in relation to @p total, rounded to @p digits.

    @note If @p total is zero, the function returns 0.0 to avoid division by zero.

    @param[in] value The value to compute the percentage for
    @param[in] total The total value to compute the percentage against
    @param[in] digits The number of digits to round the result to
    @return The percentage of @p value in relation to @p total, rounded to @p digits
    @throw OpenMS::Exception::InvalidValue if @p value or @p total is negative.


    \code
      auto one_third = 1.0/3;
      percentOf(one_third, 1.0, 2) // returns 33.33
    \endcode
  */
  template<typename T>
  static double percentOf(T value, T total, int digits)
  {
    if (value < 0) { throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value must be non-negative", String(value)); }
    if (total < 0) { throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Total must be non-negative", String(total)); }
    if (total <= 0) // avoid float equality compare
    {
      return 0.0; // avoid division by zero
    }
    return roundTo(value * 100.0 / total, digits);
  }

  /**
      @brief Returns if @p a is approximately equal @p b , allowing a tolerance of @p tol

      @ingroup MathFunctionsMisc
  */
  static bool approximatelyEqual(double a, double b, double tol)
  {
    return std::fabs(a - b) <= tol;
  }

  /**
    @brief Returns the greatest common divisor (gcd) of two numbers by applying the Euclidean algorithm.
    @param[in] a A number.
    @param[in] b A number.
    @return The greatest common divisor.
    @see gcd(T a, T b, T& a1, T& b1)
    @ingroup MathFunctionsMisc
   */
  template<typename T>
  static T gcd(T a, T b)
  {
    T c;
    while (b != 0)
    {
      c = a % b;
      a = b;
      b = c;
    }
    return a;
  }

  /**
   @brief Returns the greatest common divisor by applying the extended Euclidean algorithm (Knuth TAoCP vol. 2, p342).
   Calculates u1, u2 and u3 (which is returned) so that a * u1 + b * u2 = u3 = gcd(a, b, u1, u2)

   @param[in] a A number.
   @param[in] b A number.
   @param[out] u1 A reference to the number to be returned (see the above formula).
   @param[out] u2 A reference to the number to be returned (see the above formula).
   @return The greatest common divisor.
   @see gcd(T, T)
   @ingroup MathFunctionsMisc
   */
  template<typename T>
  static T gcd(T a, T b, T& u1, T& u2)
  {
    u1 = 1;
    u2 = 0;
    T u3 = a;

    T v1 = 0;
    T v2 = 1;
    T v3 = b;

    while (v3 != 0)
    {
      T q = u3 / v3;
      T t1 = u1 - v1 * q;
      T t2 = u2 - v2 * q;
      T t3 = u3 - v3 * q;

      u1 = v1;
      u2 = v2;
      u3 = v3;

      v1 = t1;
      v2 = t2;
      v3 = t3;
    }

    return u3;
  }

  /**
    @brief Compute parts-per-million of two @em m/z values.

    The returned ppm value can be either positive (mz_obs > mz_ref) or negative (mz_obs < mz_ref)!

    @param[in] mz_obs Observed (experimental) m/z
    @param[in] mz_ref Reference (theoretical) m/z
    @return The ppm value
  */
  template<typename T>
  static T getPPM(T mz_obs, T mz_ref)
  {
    return (mz_obs - mz_ref) / mz_ref * 1e6;
  }

  /**
    @brief Compute absolute parts-per-million of two @em m/z values.

    The returned ppm value is always >= 0.

    @param[in] mz_obs Observed (experimental) m/z
    @param[in] mz_ref Reference (theoretical) m/z
    @return The absolute ppm value
  */
  template<typename T>
  static T getPPMAbs(T mz_obs, T mz_ref)
  {
    return std::fabs(getPPM(mz_obs, mz_ref));
  }

  /**
    @brief Compute the mass diff in [Th], given a ppm value and a reference point.

    The returned mass diff can be either positive (ppm > 0) or negative (ppm < 0)!

    @param[in] ppm Parts-per-million error
    @param[in] mz_ref Reference m/z
    @return The mass diff in [Th]
  */
  template<typename T>
  static T ppmToMass(T ppm, T mz_ref)
  {
    return (ppm / T(1e6)) * mz_ref;
  }

  /*
    @brief Compute the absolute mass diff in [Th], given a ppm value and a reference point.

    The returned mass diff is always positive!

    @param[in] ppm Parts-per-million error
    @param[in] mz_ref Reference m/z
    @return The absolute mass diff in [Th]
  */
  template<typename T>
  static T ppmToMassAbs(T ppm, T mz_ref)
  {
    return std::fabs(ppmToMass(ppm, mz_ref));
  }

  /**
    @brief Return tolerance window around @p val given tolerance @p tol

    Note that when ppm is used, the window is not symmetric. In this case,
    (right - @p val) > (@p val - left), i.e., the tolerance window also
    includes the largest value x which still has @p val in *its* tolerance
    window for the given ppms, so the compatibility relation is symmetric.

    @param[in] val Value
    @param[in] tol Tolerance
    @param[in] ppm Whether @p tol is in ppm or absolute
    @return Tolerance window boundaries
  */
  static std::pair<double, double> getTolWindow(double val, double tol, bool ppm)
  {
    double left, right;

    if (ppm)
    {
      left = val - val * tol * 1e-6;
      right = val / (1.0 - tol * 1e-6);
    }
    else
    {
      left = val - tol;
      right = val + tol;
    }

    return std::make_pair(left, right);
  }

  /**
     @brief Returns the value of the @p q th quantile (0-1) in a sorted non-empty vector @p x
  */
  template<typename T1>
  static typename T1::value_type quantile(const T1& x, double q)
  {
    if (x.empty()) throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Quantile requested from empty container.");
    if (q < 0.0) q = 0.;
    if (q > 1.0) q = 1.;

    const auto n = x.size();
    const auto id = std::max(0., n * q - 1); // -1 for c++ index starting at 0
    const auto lo = floor(id);
    const auto hi = ceil(id);
    const auto qs = x[lo];
    const auto h = (id - lo);

    return (1.0 - h) * qs + h * x[hi];
  }

  // portable random shuffle
  class OPENMS_DLLAPI RandomShuffler
  {
  public:
    explicit RandomShuffler(int seed): rng_(boost::mt19937_64(seed))
    {
    }

    explicit RandomShuffler(const boost::mt19937_64& mt_rng): rng_(mt_rng)
    {
    }

    RandomShuffler() = default;
    ~RandomShuffler() = default;

    boost::mt19937_64 rng_;
    template<class RandomAccessIterator>
    void portable_random_shuffle(RandomAccessIterator first, RandomAccessIterator last)
    {
      for (auto i = (last - first) - 1; i > 0; --i) // OMS_CODING_TEST_EXCLUDE
      {
        boost::uniform_int<decltype(i)> d(0, i);
        std::swap(first[i], first[d(rng_)]);
      }
    }

    void seed(uint64_t val)
    {
      rng_.seed(val);
    }
  };

  /**
   * @brief Calculate logarithm of binomial coefficient C(n,k) using log-gamma function
   * 
   * @param[in] n Total number of items
   * @param[in] k Number of items to choose
   * @return Natural logarithm of binomial coefficient C(n,k)
   * @throws std::invalid_argument if k > n
   */
  static double log_binomial_coef(unsigned n, unsigned k) 
  {
    // Handle edge cases for improved numerical stability
    if (k > n) 
    {
      throw std::invalid_argument("k cannot be greater than n in binomial coefficient");
    }
    
    if (k == 0 || k == n) 
    {
      return 0.0;  // log(1) = 0
    }
    
    // Use symmetry to minimize computation for large k
    if (k > n / 2) 
    {
      k = n - k;
    }
    
    return boost::math::lgamma(n + 1.0) - boost::math::lgamma(k + 1.0) - boost::math::lgamma(n - k + 1.0);
  }

  /**
   * @brief Log-sum-exp operation for numerical stability
   * 
   * @param[in] x First logarithmic value
   * @param[in] y Second logarithmic value
   * @return Natural logarithm of (exp(x) + exp(y))
   */
  static double log_sum_exp(double x, double y) 
  {
    // Handle infinite cases
    if (std::isinf(x) && x < 0) return y;
    if (std::isinf(y) && y < 0) return x;
    
    // Use the maximum value for numerical stability
    double max_val = std::max(x, y);
    return max_val + std::log(std::exp(x - max_val) + std::exp(y - max_val));
  }

  /**
   * @brief Calculate binomial cumulative distribution function P(X ≥ n)
   * 
   * Calculates P(X ≥ n) for a binomial distribution with parameters N and p,
   * using numerically stable algorithms in the log domain to handle large values.
   * 
   * @param[in] N Total number of trials
   * @param[in] n Minimum number of successes
   * @param[in] p Probability of success in each trial
   * @return Probability P(X ≥ n) for binomial distribution B(N,p)
   * @throws std::invalid_argument if parameters are invalid
   */
  static double binomial_cdf_complement(unsigned N, unsigned n, double p)
  {
    if (p < 0.0 || p > 1.0)
    {
      throw std::invalid_argument("Probability p must be between 0 and 1");
    }
    if (n > N)
    {
      throw std::invalid_argument("n cannot be greater than N");
    }

    if (n == 0)   return 1.0;                // P(X ≥ 0) = 1
    if (p == 0.0) return (n == 0) ? 1.0 : 0.0;
    if (p == 1.0) return 1.0;               // all mass at N

    const boost::math::binomial_distribution<double> dist(N, p);
    return boost::math::cdf(boost::math::complement(dist, n - 1));
  }

  //
  // Statistical functions (formerly in StatisticFunctions.h)
  //

  /**
    @brief Result of adaptiveQuantile computation.

    Fields:
      - blended       : the final blended (adaptive) quantile
      - half_raw      : raw q-quantile of values
      - half_rob      : q-quantile after IQR-winsorization
      - upper_fence   : Tukey upper fence (Q3 + k*IQR), +inf if undefined
      - tail_fraction : fraction of values above upper_fence
      - weight        : blend weight w in [0,1] (0=robust, 1=raw)
  */
  struct AdaptiveQuantileResult
  {
    double blended{0.0};
    double half_raw{0.0};
    double half_rob{0.0};
    double upper_fence{std::numeric_limits<double>::infinity()};
    double tail_fraction{0.0};
    double weight{0.0};
  };

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
    if ((begin_b == end_b) ^ (begin_a == end_a))
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

     @param[in] begin Start of range
     @param[in] end End of range (past-the-end iterator)
     @param[in] sorted Is the range already sorted? If not, it will be sorted.
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

    @param[in] begin Start of range
    @param[in] end End of range (past-the-end iterator)
    @param[in] median_of_numbers The precomputed median of range @p begin - @p end.
    @return the MAD

    @ingroup MathFunctionsStatistics

  */
  template <typename IteratorType>
  static double MAD(IteratorType begin, IteratorType end, double median_of_numbers)
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

    @param[in] begin Start of range
    @param[in] end End of range (past-the-end iterator)
    @param[in] mean_of_numbers The precomputed mean of range @p begin - @p end.
    @return the MeanAbsoluteDeviation

    @ingroup MathFunctionsStatistics

  */
  template <typename IteratorType>
  static double MeanAbsoluteDeviation(IteratorType begin, IteratorType end, double mean_of_numbers)
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

     @param[in] begin Start of range
     @param[in] end End of range (past-the-end iterator)
     @param[in] sorted Is the range already sorted? If not, it will be sorted.

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

     @param[in] begin Start of range
     @param[in] end End of range (past-the-end iterator)
     @param[in] sorted Is the range already sorted? If not, it will be sorted.

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
    @brief Calculates the q-quantile (0 <= q <= 1) of a *sorted* range of values.

    Assumes the range [begin, end) is already sorted in non-decreasing order.
    Uses the common "Type 7" definition (linear interpolation):

      pos = q * (n - 1)
      idx = floor(pos),  frac = pos - idx
      quantile = (1 - frac) * x[idx] + frac * x[idx + 1]

    Exact endpoints:
      - q == 0 returns the first (minimum) element
      - q == 1 returns the last (maximum) element

    @param[in] begin  Start of range
    @param[in] end    End of range (past-the-end iterator)
    @param[in] q      Quantile in [0, 1]

    @pre Input range must be sorted ascending.

    @exception Exception::InvalidRange is thrown if the range is NULL or empty.
    @exception Exception::InvalidValue is thrown if q is outside [0, 1].

    @ingroup MathFunctionsStatistics
  */
  template <typename IteratorType>
  static double quantileFromSortedRange(IteratorType begin, IteratorType end, double q)
  {
    OPENMS_PRECONDITION(std::is_sorted(begin, end),
                        "Math::quantileFromSortedRange expects a sorted range. Sort before calling.");

    checkIteratorsNotNULL(begin, end);

    const Size n = std::distance(begin, end);
    if (n == 0)
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    if (q < 0.0 || q > 1.0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "q must be in [0,1]", String(q));
    }
    if (n == 1) return static_cast<double>(*begin);

    const double pos = q * static_cast<double>(n - 1);
    const Size i = static_cast<Size>(std::floor(pos));
    const double frac = pos - static_cast<double>(i);

    const auto it_i = begin + static_cast<typename std::iterator_traits<IteratorType>::difference_type>(i);
    if (frac == 0.0) return static_cast<double>(*it_i);

    const auto it_ip1 = it_i + 1;
    return (1.0 - frac) * static_cast<double>(*it_i) + frac * static_cast<double>(*it_ip1);
  }

  /**
    @brief Tukey upper fence (UF) for outlier detection.

    Computes Q3 + k * IQR on the (finite) values in [begin,end).
    If there are too few values or IQR ≤ 0, returns +infinity.

    References: J. W. Tukey (1977). Exploratory Data Analysis.

    @tparam IteratorType  input iterator over arithmetic values
    @param[in] begin          start iterator
    @param[in] end            past-the-end iterator
    @param[in] k              Tukey factor (default 1.5)
    @return               upper fence (Q3 + k*IQR) or +infinity if undefined
  */
  template <typename IteratorType>
  static double tukeyUpperFence(IteratorType begin, IteratorType end, double k = 1.5)
  {
      std::vector<double> v;
      v.reserve(std::distance(begin, end));
      for (auto it = begin; it != end; ++it)
      {
        if (std::isfinite(*it)) v.push_back(static_cast<double>(*it));
      }
      if (v.size() < 4) return std::numeric_limits<double>::infinity();

      std::sort(v.begin(), v.end());
      const double q1  = quantileFromSortedRange(v.begin(), v.end(), 0.25);
      const double q3  = quantileFromSortedRange(v.begin(), v.end(), 0.75);
      const double iqr = q3 - q1;
      if (!(iqr > 0.0)) return std::numeric_limits<double>::infinity();

      return q3 + k * iqr;
  }

  /**
    @brief Fraction of values above a threshold.

    @tparam IteratorType  input iterator over arithmetic values
    @param[in] begin          start iterator
    @param[in] end            past-the-end iterator
    @param[in] threshold      threshold T
    @return               (# { x > T } / N), ignoring non-finite x
  */
  template <typename IteratorType>
  static double tailFractionAbove(IteratorType begin, IteratorType end, double threshold)
  {
      size_t n = 0, n_tail = 0;
      for (auto it = begin; it != end; ++it)
      {
        const double x = static_cast<double>(*it);
        if (!std::isfinite(x)) continue;
        ++n;
        if (x > threshold) ++n_tail;
      }
      return (n == 0) ? 0.0 : static_cast<double>(n_tail) / static_cast<double>(n);
  }

  /**
    @brief Quantile after winsorizing at an upper fence.

    Copies the (finite) values in [begin,end), caps them at @p upper_fence
    (and at 0 on the lower side, which is convenient for absolute residuals),
    then returns the requested quantile.

    If @p upper_fence is not finite, this falls back to the raw quantile.

    References: J. W. Tukey (1962). The Future of Data Analysis.

    @tparam IteratorType  input iterator over arithmetic values
    @param[in] begin          start iterator
    @param[in] end            past-the-end iterator
    @param[in] q              quantile in [0,1]
    @param[in] upper_fence    winsorization cap (Q3+k*IQR), or +inf to disable
    @return               winsorized quantile
  */
  template <typename IteratorType>
  static double winsorizedQuantile(IteratorType begin, IteratorType end, double q, double upper_fence)
  {
      std::vector<double> v;
      v.reserve(std::distance(begin, end));
      for (auto it = begin; it != end; ++it)
      {
        const double x = static_cast<double>(*it);
        if (!std::isfinite(x)) continue;
        v.push_back(x);
      }
      if (v.empty()) return 0.0;

      if (std::isfinite(upper_fence))
      {
        for (double& x : v)
        {
          if (x > upper_fence) x = upper_fence;
          if (x < 0.0) x = 0.0; // defensive; useful when passing |residual|
        }
      }
      std::sort(v.begin(), v.end());
      return quantileFromSortedRange(v.begin(), v.end(), q);
  }

  /**
    @brief Adaptive quantile that blends RAW and IQR-winsorized quantiles
           based on tail density beyond the Tukey upper fence.

    Let UF = Q3 + k*IQR on the (finite) inputs. Compute:
     - half_raw = quantile(values, q)
     - half_rob = winsorizedQuantile(values, q, UF)
     - r       = fraction(values > UF)

    Blend with weight w(r):
      r ≤ r_sparse  -> w=0 (use robust)
      r ≥ r_dense   -> w=1 (use raw)
      otherwise     -> linear interpolation between 0 and 1

    Returned value = (1-w)*half_rob + w*half_raw.

    This keeps windows stable when outliers are sparse, while respecting
    genuinely broad tails (dense outliers) by leaning toward the raw quantile.

    References:
        - J. W. Tukey (1962). The Future of Data Analysis.
        - J. W. Tukey (1977). Exploratory Data Analysis.
        - R. J. Hyndman, Y. Fan (1996). Sample Quantiles in Statistical Packages

    @tparam IteratorType   input iterator over arithmetic values
    @param[in] begin           start iterator
    @param[in] end             past-the-end iterator
    @param[in] q               target quantile in [0,1] (e.g., 0.99 for 99% half-width)
    @param[in] k               Tukey factor (default 1.5)
    @param[in] r_sparse        tail density below which robust wins (default 0.01 = 1%)
    @param[in] r_dense         tail density above which raw wins (default 0.10 = 10%)

    @return                AdaptiveQuantileResult
  */
  template <typename IteratorType>
  static AdaptiveQuantileResult adaptiveQuantile(IteratorType begin, IteratorType end, double q,
                          double k = 1.5,
                          double r_sparse = 0.01,
                          double r_dense  = 0.10)
  {
      AdaptiveQuantileResult res;

      // Copy finite values
      std::vector<double> v;
      v.reserve(std::distance(begin, end));
      for (auto it = begin; it != end; ++it)
      {
        if (std::isfinite(*it)) v.push_back(static_cast<double>(*it));
      }
      if (v.empty())
      {
        return res;
      }

      std::sort(v.begin(), v.end());
      const double half_raw = quantileFromSortedRange(v.begin(), v.end(), q);

      // Robust path (winsorization at Tukey fence)
      const double uf       = tukeyUpperFence(v.begin(), v.end(), k);
      const double r        = std::isfinite(uf) ? tailFractionAbove(v.begin(), v.end(), uf) : 0.0;
      const double half_rob = winsorizedQuantile(v.begin(), v.end(), q, uf);

      // Blend weight w(r)
      double w = 0.0;
      if (r_dense <= r_sparse)
      {
        w = (r > r_sparse) ? 1.0 : 0.0;
      }
      else
      {
        const double t = (r - r_sparse) / (r_dense - r_sparse);
        w = std::max(0.0, std::min(1.0, t));
      }

      res.half_raw      = half_raw;
      res.half_rob      = half_rob;
      res.upper_fence   = uf;
      res.tail_fraction = r;
      res.weight        = w;
      res.blended       = (1.0 - w) * half_rob + w * half_raw;
      return res;
  }

  /**
     @brief Calculates the variance of a range of values

The @p mean_value can be provided explicitly to save computation time. If left at default, it will be computed internally.

     @exception Exception::InvalidRange is thrown if the range is empty

     @ingroup MathFunctionsStatistics
  */
  template <typename IteratorType>
  static double variance(IteratorType begin, IteratorType end,
                         double mean_value = std::numeric_limits<double>::max())
  {
    checkIteratorsNotNULL(begin, end);
    double sum_value = 0.0;
    if (mean_value == std::numeric_limits<double>::max())
    {
      mean_value = Math::mean(begin, end);
    }
    for (IteratorType iter=begin; iter!=end; ++iter)
    {
      double diff = *iter - mean_value;
      sum_value += diff * diff;
    }
    return sum_value / (std::distance(begin, end)-1);
  }

  /**
     @brief Calculates the standard deviation of a range of values.

     The @p mean_value can be provided explicitly to save computation time. If left at default, it will be computed internally.

     @exception Exception::InvalidRange is thrown if the range is empty

     @ingroup MathFunctionsStatistics
  */
  template <typename IteratorType>
  static double sd(IteratorType begin, IteratorType end,
                   double mean_value = std::numeric_limits<double>::max())
  {
    checkIteratorsNotNULL(begin, end);
    return std::sqrt( variance(begin, end, mean_value) );
  }

  /**
     @brief Calculates the absolute deviation of a range of values

     @exception Exception::InvalidRange is thrown if the range is empty

     @ingroup MathFunctionsStatistics
  */
  template <typename IteratorType>
  static double absdev(IteratorType begin, IteratorType end,
                       double mean_value = std::numeric_limits<double>::max())
  {
    checkIteratorsNotNULL(begin, end);
    double sum_value = 0.0;
    if (mean_value == std::numeric_limits<double>::max())
    {
      mean_value = Math::mean(begin, end);
    }
    for (IteratorType iter=begin; iter!=end; ++iter)
    {
      sum_value += *iter - mean_value;
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
     @brief Calculates the root mean square error (RMSE) for the values in
       [begin_a, end_a) and [begin_b, end_b)

     Computes the square root of the mean of the squared differences between the
      two iterator ranges (i.e., RMSE = sqrt(MSE)). .

     @exception Exception::InvalidRange is thrown if the iterator ranges are not
     of the same length or are empty.

     @ingroup MathFunctionsStatistics
  */
  template <typename IteratorType1, typename IteratorType2>
  static double rootMeanSquareError(IteratorType1 begin_a, IteratorType1 end_a,
                                    IteratorType2 begin_b, IteratorType2 end_b)
  {
    return std::sqrt(meanSquareError(begin_a, end_a, begin_b, end_b));
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

private:
  /// Private constructor to prevent instantiation (all methods are static)
  Math() = delete;
}; // class Math
} // namespace OpenMS
