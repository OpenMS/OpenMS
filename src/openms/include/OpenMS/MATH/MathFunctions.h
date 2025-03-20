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
#include <OpenMS/KERNEL/RangeManager.h>

#include <boost/random/mersenne_twister.hpp> // for mt19937_64
#include <boost/random/uniform_int.hpp>
#include <cmath>
#include <utility> // for std::pair
#include <vector>

namespace OpenMS
{
/**
  @brief %Math namespace.

  Contains mathematical auxiliary functions.

  @ingroup Concept
*/
namespace Math
{

  /**
    @brief Given an interval/range and a new value, extend the range to include the new value if needed

    @param min The current minimum of the range
    @param max The current maximum of the range
    @param value The new value which may extend the range
    @return true if the range was modified
  */
  template<typename T>
  bool extendRange(T& min, T& max, const T& value)
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
  bool contains(T value, T min, T max)
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
  inline std::pair<double, double> zoomIn(const double left, const double right, const float factor, const float align)
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
    
    @param min The minimum of the range; must be smaller than @p max
    @param max The maximum of the range
    @param number_of_bins How many bins should the range be divided into? Must be 1 or larger
    @param extend_margin Overlap of neighboring bins (=0 for no overlap). Negative values will shrink the range (feature).
    @return Vector with @p number_of_bins elements, each representing the margins of one bin

    @throws OpenMS::Precondition if `min >= max` or `number_of_bins == 0`
  */
  inline BinContainer createBins(double min, double max, uint32_t number_of_bins, double extend_margin = 0)
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
  inline double ceilDecimal(double x, int decPow)
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
  inline double roundDecimal(double x, int decPow)
  {
    if (x > 0) return (floor(0.5 + x / pow(10.0, decPow))) * pow(10.0, decPow);

    return -((floor(0.5 + fabs(x) / pow(10.0, decPow))) * pow(10.0, decPow));
  }

  /**
      @brief transforms point @p x of interval [left1,right1] into interval [left2,right2]

      @ingroup MathFunctionsMisc
  */
  inline double intervalTransformation(double x, double left1, double right1, double left2, double right2)
  {
    return left2 + (x - left1) * (right2 - left2) / (right1 - left1);
  }

  /**
      @brief Transforms a number from linear to log10 scale. Avoids negative logarithms by adding 1.

      @param x The number to transform

      @ingroup MathFunctionsMisc
  */
  inline double linear2log(double x)
  {
    return log10(x + 1); //+1 to avoid negative logarithms
  }

  /**
      @brief Transforms a number from log10 to to linear scale. Subtracts the 1 added by linear2log(double)

      @param x The number to transform

      @ingroup MathFunctionsMisc
  */
  inline double log2linear(double x)
  {
    return pow(10, x) - 1;
  }

  /**
      @brief Returns true if the given integer is odd

      @ingroup MathFunctionsMisc
  */
  inline bool isOdd(UInt x)
  {
    return (x & 1) != 0;
  }

  /**
      @brief Rounds the value

      @ingroup MathFunctionsMisc
  */
  template<typename T>
  T round(T x)
  {
    if (x >= T(0)) { return T(floor(x + T(0.5))); }
    else { return T(ceil(x - T(0.5))); }
  }

  /**
      @brief Returns if @p a is approximately equal @p b , allowing a tolerance of @p tol

      @ingroup MathFunctionsMisc
  */
  inline bool approximatelyEqual(double a, double b, double tol)
  {
    return std::fabs(a - b) <= tol;
  }

  /**
    @brief Returns the greatest common divisor (gcd) of two numbers by applying the Euclidean algorithm.
    @param a A number.
    @param b A number.
    @return The greatest common divisor.
    @see gcd(T a, T b, T& a1, T& b1)
    @ingroup MathFunctionsMisc
   */
  template<typename T>
  T gcd(T a, T b)
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

   @param a A number.
   @param b A number.
   @param u1 A reference to the number to be returned (see the above formula).
   @param u2 A reference to the number to be returned (see the above formula).
   @return The greatest common divisor.
   @see gcd(T, T)
   @ingroup MathFunctionsMisc
   */
  template<typename T>
  T gcd(T a, T b, T& u1, T& u2)
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

    @param mz_obs Observed (experimental) m/z
    @param mz_ref Reference (theoretical) m/z
    @return The ppm value
  */
  template<typename T>
  T getPPM(T mz_obs, T mz_ref)
  {
    return (mz_obs - mz_ref) / mz_ref * 1e6;
  }

  /**
    @brief Compute absolute parts-per-million of two @em m/z values.

    The returned ppm value is always >= 0.

    @param mz_obs Observed (experimental) m/z
    @param mz_ref Reference (theoretical) m/z
    @return The absolute ppm value
  */
  template<typename T>
  T getPPMAbs(T mz_obs, T mz_ref)
  {
    return std::fabs(getPPM(mz_obs, mz_ref));
  }

  /**
    @brief Compute the mass diff in [Th], given a ppm value and a reference point.

    The returned mass diff can be either positive (ppm > 0) or negative (ppm < 0)!

    @param ppm Parts-per-million error
    @param mz_ref Reference m/z
    @return The mass diff in [Th]
  */
  template<typename T>
  T ppmToMass(T ppm, T mz_ref)
  {
    return (ppm / T(1e6)) * mz_ref;
  }

  /*
    @brief Compute the absolute mass diff in [Th], given a ppm value and a reference point.

    The returned mass diff is always positive!

    @param ppm Parts-per-million error
    @param mz_ref Reference m/z
    @return The absolute mass diff in [Th]
  */
  template<typename T>
  T ppmToMassAbs(T ppm, T mz_ref)
  {
    return std::fabs(ppmToMass(ppm, mz_ref));
  }

  /**
    @brief Return tolerance window around @p val given tolerance @p tol

    Note that when ppm is used, the window is not symmetric. In this case,
    (right - @p val) > (@p val - left), i.e., the tolerance window also
    includes the largest value x which still has @p val in *its* tolerance
    window for the given ppms, so the compatibility relation is symmetric.

    @param val Value
    @param tol Tolerance
    @param ppm Whether @p tol is in ppm or absolute
    @return Tolerance window boundaries
  */
  inline std::pair<double, double> getTolWindow(double val, double tol, bool ppm)
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
  typename T1::value_type quantile(const T1& x, double q)
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
} // namespace Math
} // namespace OpenMS
