// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <iostream>

namespace OpenMS
{

  /**
    @brief Wrapper class to implement output with appropriate precision.  See precisionWrapper().
  */
  template <typename FloatingPointType>
  struct PrecisionWrapper
  {
    /// Constructor.  Note: Normally you will prefer to use the "make"-function precisionWrapper(), which see.
    explicit PrecisionWrapper(const FloatingPointType rhs) :
      ref_(rhs) {}
    PrecisionWrapper(const PrecisionWrapper & rhs) :
      ref_(rhs.ref_) {}
    FloatingPointType const ref_;
private:
    PrecisionWrapper();     // intentionally not implemented
  };

  /**
    @brief Wrapper function that sets the appropriate precision for output
    temporarily.  The original precision is restored afterwards so that no side
    effects remain.  This is a "make"-function that deduces the typename
    FloatingPointType from its argument and returns a
    PrecisionWrapper<FloatingPointType>.

    Example:
    @code
    std::cout
    << 0.1234567890123456789f << ' ' << 0.1234567890123456789 << ' ' << 0.1234567890123456789l << '\n'
    << precisionWrapper(0.1234567890123456789f) << '\n' // float
    << 0.1234567890123456789f << ' ' << 0.1234567890123456789 << ' ' << 0.1234567890123456789l << '\n'
    << precisionWrapper(0.1234567890123456789) << '\n' // double
    << 0.1234567890123456789f << ' ' << 0.1234567890123456789 << ' ' << 0.1234567890123456789l << '\n'
    << precisionWrapper(0.1234567890123456789l) << '\n' // long double
    << 0.1234567890123456789f << ' ' << 0.1234567890123456789 << ' ' << 0.1234567890123456789l << '\n';
    @endcode
    Result:
    @code
    0.123457 0.123457 0.123457
    0.123457
    0.123457 0.123457 0.123457
    0.123456789012346
    0.123457 0.123457 0.123457
    0.123456789012345679
    0.123457 0.123457 0.123457
    @endcode

    Note: Unfortunately we cannot return a const& - this will change when rvalue
    references become part of the new C++ standard.  In the meantime, we need a
    copy constructor for PrecisionWrapper.
  */
  template <typename FloatingPointType>
  inline const PrecisionWrapper<FloatingPointType> precisionWrapper(const FloatingPointType rhs)
  {
    return PrecisionWrapper<FloatingPointType>(rhs);
  }

  /// Output operator for a PrecisionWrapper. Specializations are defined for float, double, long double.
  template <typename FloatingPointType>
  inline std::ostream & operator<<(std::ostream & os, const PrecisionWrapper<FloatingPointType> & rhs)
  {
    // manual conversion to String is much faster than ostreams internal conversion
    // (this calls the correct overload for extra precision)
    String s(rhs.ref_, true);
    os << s;
    return os;
  }
} // namespace OpenMS

