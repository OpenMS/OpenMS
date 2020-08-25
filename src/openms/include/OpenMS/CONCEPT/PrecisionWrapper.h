// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

