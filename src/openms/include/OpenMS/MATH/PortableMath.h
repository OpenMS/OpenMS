// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Kohlbacher $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------
//
// exp() and log() are ports of fdlibm 5.3 (e_exp.c, e_log.c; netlib.org/fdlibm), which carry this notice:
//
// ====================================================
// Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//
// Developed at SunPro, a Sun Microsystems, Inc. business.
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
// ====================================================

// exp and log with the same bits on every platform.
//
// The platform libm does not give them: Apple's libm, glibc and the Windows UCRT round exp and log differently in the
// last bit, and in the mass trace detection those bits move the +-3 sd trace windows, so one input gave different
// traces on macOS, Linux and Windows. These are fdlibm's functions (error < 1 ulp), built from IEEE + - * / and exact
// bit and exponent operations only, which every IEEE-754 double platform evaluates alike as long as nothing is
// contracted into an FMA. Contraction is switched off inside the functions for clang and GCC (OpenMS builds with
// -ffp-contract=off anyway); MSVC does not contract under /fp:precise.
#pragma once

#include <cstdint>
#include <cstring>
#include <cmath>
#include <limits>

#if defined(__GNUC__) && ! defined(__clang__)
  #pragma GCC push_options
  #pragma GCC optimize("fp-contract=off")
#endif
#if defined(__clang__)
  #define OPENMS_PORTABLEMATH_NO_FMA _Pragma("clang fp contract(off)")
#else
  #define OPENMS_PORTABLEMATH_NO_FMA
#endif

/// @brief exp and log with the same bits on every platform (fdlibm)
namespace OpenMS::PortableMath
{
/// @brief Bit-access helpers for IEEE-754 doubles used by the fdlibm-derived functions.
namespace detail
{
  /**
    @brief Returns the high 32 bits (sign, exponent, upper mantissa) of @p x.

    @param[in] x The double to inspect.
    @return The high word of the IEEE-754 representation of @p x, as a signed 32-bit integer.
  */
  inline int32_t hi(double x)
  {
    uint64_t b;
    std::memcpy(&b, &x, 8);
    return (int32_t)(b >> 32);
  }
  /**
    @brief Returns the low 32 bits (lower mantissa) of @p x.

    @param[in] x The double to inspect.
    @return The low word of the IEEE-754 representation of @p x.
  */
  inline uint32_t lo(double x)
  {
    uint64_t b;
    std::memcpy(&b, &x, 8);
    return (uint32_t)b;
  }
  /**
    @brief Returns @p x with its high 32 bits replaced by @p h; the low 32 bits are kept.

    @param[in] x The double whose low word is kept.
    @param[in] h The new high word.
    @return The double built from @p h and the low word of @p x.
  */
  inline double withHi(double x, int32_t h)
  {
    uint64_t b;
    std::memcpy(&b, &x, 8);
    b = ((uint64_t)(uint32_t)h << 32) | (b & 0xffffffffULL);
    std::memcpy(&x, &b, 8);
    return x;
  }
} // namespace detail

/**
  @brief e^x (fdlibm e_exp.c)

  Error below 1 ulp. Returns NaN for NaN, +inf above the overflow threshold and 0 below the underflow threshold.

  @param[in] x The exponent.
  @return e raised to @p x.
*/
inline double exp(double x)
{
  OPENMS_PORTABLEMATH_NO_FMA
  const double o_threshold = 7.09782712893383973096e+02, u_threshold = -7.45133219101941108420e+02;
  const double ln2HI[2] = {6.93147180369123816490e-01, -6.93147180369123816490e-01};
  const double ln2LO[2] = {1.90821492927058770002e-10, -1.90821492927058770002e-10};
  const double invln2 = 1.44269504088896338700e+00, halF[2] = {0.5, -0.5}, huge = 1.0e+300;
  const double P1 = 1.66666666666666019037e-01, P2 = -2.77777777770155933842e-03, P3 = 6.61375632143793436117e-05, P4 = -1.65339022054652515390e-06,
               P5 = 4.13813679705723846039e-08;
  double hi_ = 0.0, lo_ = 0.0;
  int k = 0;
  int32_t hx = detail::hi(x);
  const int xsb = (hx >> 31) & 1;
  hx &= 0x7fffffff;
  if (hx >= 0x40862E42) // |x| >= 709.78...
  {
    if (hx >= 0x7ff00000)
    {
      if (((hx & 0xfffff) | detail::lo(x)) != 0)
      {
        return x + x; // NaN
      }
      return xsb == 0 ? x : 0.0; // exp(+-inf) = inf, 0
    }
    if (x > o_threshold) { return std::numeric_limits<double>::infinity(); }
    if (x < u_threshold) { return 0.0; }
  }
  if (hx > 0x3fd62e42) // |x| > 0.5 ln2: reduce
  {
    if (hx < 0x3FF0A2B2)
    {
      hi_ = x - ln2HI[xsb];
      lo_ = ln2LO[xsb];
      k = 1 - xsb - xsb;
    }
    else
    {
      k = (int)(invln2 * x + halF[xsb]);
      const double t = k;
      hi_ = x - t * ln2HI[0];
      lo_ = t * ln2LO[0];
    }
    x = hi_ - lo_;
  }
  else if (hx < 0x3e300000) // |x| < 2^-28
  {
    if (huge + x > 1.0) { return 1.0 + x; }
  }
  const double t = x * x;
  const double c = x - t * (P1 + t * (P2 + t * (P3 + t * (P4 + t * P5))));
  if (k == 0) { return 1.0 - ((x * c) / (c - 2.0) - x); }
  const double y = 1.0 - ((lo_ - (x * c) / (2.0 - c)) - hi_);
  // fdlibm's scaling, not std::ldexp: a normal result takes k into its exponent bits (exact); a subnormal one is
  // scaled to a normal number first and rounded once by an IEEE multiplication. ldexp's rounding into the
  // subnormals is the platform libm's, and macOS 14's differs from glibc's and macOS 26's.
  if (k >= -1021) { return detail::withHi(y, detail::hi(y) + (k << 20)); }
  return detail::withHi(y, detail::hi(y) + ((k + 1000) << 20)) * 9.33263618503218878990e-302; // * 2^-1000
}

/**
  @brief natural logarithm (fdlibm e_log.c)

  Error below 1 ulp. Returns -inf for +-0, NaN for negative @p x and NaN, and +inf for +inf.

  @param[in] x The argument.
  @return The natural logarithm of @p x.
*/
inline double log(double x)
{
  OPENMS_PORTABLEMATH_NO_FMA
  const double ln2_hi = 6.93147180369123816490e-01, ln2_lo = 1.90821492927058770002e-10, two54 = 1.80143985094819840000e+16;
  const double Lg1 = 6.666666666666735130e-01, Lg2 = 3.999999999940941908e-01, Lg3 = 2.857142874366239149e-01, Lg4 = 2.222219843214978396e-01,
               Lg5 = 1.818357216161805012e-01, Lg6 = 1.531383769920937332e-01, Lg7 = 1.479819860511658591e-01;
  int32_t hx = detail::hi(x);
  const uint32_t lx = detail::lo(x);
  int32_t k = 0;
  if (hx < 0x00100000) // x < 2^-1022
  {
    if (((hx & 0x7fffffff) | lx) == 0) { return -std::numeric_limits<double>::infinity(); }
    if (hx < 0) { return std::numeric_limits<double>::quiet_NaN(); }
    k -= 54;
    x *= two54;
    hx = detail::hi(x); // subnormal: scale up
  }
  if (hx >= 0x7ff00000) { return x + x; }
  k += (hx >> 20) - 1023;
  hx &= 0x000fffff;
  int32_t i = (hx + 0x95f64) & 0x100000;
  x = detail::withHi(x, hx | (i ^ 0x3ff00000)); // normalise x or x/2
  k += (i >> 20);
  const double f = x - 1.0;
  if ((0x000fffff & (2 + hx)) < 3) // -2^-20 <= f < 2^-20
  {
    if (f == 0.0)
    {
      if (k == 0) { return 0.0; }
      const double dk = (double)k;
      return dk * ln2_hi + dk * ln2_lo;
    }
    const double R = f * f * (0.5 - 0.33333333333333333 * f);
    if (k == 0) { return f - R; }
    const double dk = (double)k;
    return dk * ln2_hi - ((R - dk * ln2_lo) - f);
  }
  const double s = f / (2.0 + f), dk = (double)k, z = s * s, w = z * z;
  i = hx - 0x6147a;
  const int32_t j = 0x6b851 - hx;
  const double t1 = w * (Lg2 + w * (Lg4 + w * Lg6));
  const double t2 = z * (Lg1 + w * (Lg3 + w * (Lg5 + w * Lg7)));
  i |= j;
  const double R = t2 + t1;
  if (i > 0)
  {
    const double hfsq = 0.5 * f * f;
    if (k == 0) { return f - (hfsq - s * (hfsq + R)); }
    return dk * ln2_hi - ((hfsq - (s * (hfsq + R) + dk * ln2_lo)) - f);
  }
  if (k == 0) { return f - s * (f - R); }
  return dk * ln2_hi - ((s * (f - R) - dk * ln2_lo) - f);
}

} // namespace OpenMS::PortableMath

#undef OPENMS_PORTABLEMATH_NO_FMA
#if defined(__GNUC__) && ! defined(__clang__)
  #pragma GCC pop_options
#endif
