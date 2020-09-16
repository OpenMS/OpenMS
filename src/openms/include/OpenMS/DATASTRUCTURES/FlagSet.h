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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Macros.h> // for OPENMS_PRECONDITION

namespace OpenMS
{
  /**
    @brief Stores and handles combinations of enum values, e.g. a set of flags as bits flipped in an UInt64

    Conversion from the enum is computed as `pow(2, r)`.
    Thus make sure that 0 <= 'r' <=63 for all enum values.
    Multiple enum values can be computed by bitwise 'or' (operator|=)

    This class allows assignment and bit operations with itself and an object
    of type ENUM, i.e. not with any numeric types.

  **/
  template<class ENUM>
  class FlagSet
  {
  public:
    /// Constructors
    FlagSet() 
      : value_(0)
    {}

    /// C'tor from Enum
    explicit FlagSet(const ENUM& en)
      : value_(getPow_(en))
    {
    }

    FlagSet(const FlagSet& stat) = default;

    /// Assignment
    FlagSet& operator=(const FlagSet& stat) = default;

    /// no Assignment from Enum (would allow implicit conversion)
    FlagSet& operator=(const ENUM& en) = delete;


    /// Destructor (default)
    ~FlagSet() = default;

    /// Equality
    bool operator==(const FlagSet& stat) const
    {
      return (value_ == stat.value_);
    }

    /// bitwise AND
    FlagSet operator&(const ENUM& en) const
    {
      FlagSet s(*this) ;
      s &= en;
      return s;
    }

    /// bitwise AND
    FlagSet operator&(const FlagSet& rhs) const
    {
      FlagSet s(*this);
      s &= rhs;
      return s;
    }

    /// bitwise AND=
    FlagSet& operator&=(const ENUM& en)
    {
      value_ &= getPow_(en);
      return *this;
    }
  
    /// bitwise AND=
    FlagSet& operator&=(const FlagSet& rhs)
    {
      value_ &= rhs.value_;
      return *this;
    }

    /// bitwise OR
    FlagSet operator|(const ENUM& en) const
    {
      FlagSet s(*this);
      s.value_ |= getPow_(en);
      return s;
    }
    /// bitwise OR
    FlagSet operator|(const FlagSet& rhs) const
    {
      FlagSet s(*this);
      s.value_ |= rhs.value_;
      return s;
    }

    ///bitwise OR=
    FlagSet& operator|=(const ENUM& en)
    {
      value_ |= getPow_(en);
      return *this;
    }

    /// bitwise OR=
    FlagSet& operator|=(const FlagSet& rhs)
    {
      value_ |= rhs.value_;
      return *this;
    }

    /// bitwise OR (same as |)
    FlagSet operator+(const ENUM& en) const
    {
      return *this | en;
    }

    /// bitwise OR (same as |)
    FlagSet operator+(const FlagSet& en) const
    {
      return *this | en;
    }
    ///bitwise OR= (same as |=)
    FlagSet& operator+=(const ENUM& rhs)
    {
      return *this |= rhs;
    }

    /// bitwise OR= (same as |=)
    FlagSet& operator+=(const FlagSet& rhs)
    {
      return *this |= rhs;
    }

    /// remove all flags set in @p rhs from this
    FlagSet operator-(const FlagSet& rhs)
    {
      FlagSet r(*this);
      r -= rhs;
      return r;
    }

    /// remove all flags set in @p rhs from this
    FlagSet& operator-=(const FlagSet& rhs)
    {
      auto overlap = value_ & rhs.value_;
      value_ ^= overlap; // disable bits which overlap with rhs using XOR
      return *this;
    }

    /// remove flag in @p rhs from this
    FlagSet operator-(const ENUM& rhs)
    {
      FlagSet r(*this);
      r -= rhs;
      return r;
    }

    /// remove flag in @p rhs from this
    FlagSet& operator-=(const ENUM& rhs)
    {
      auto overlap = value_ & FlagSet(rhs).value_;
      value_ ^= overlap; // disable bits which overlap with rhs using XOR
      return *this;
    }

    /**
     * @brief Check if this FlagSet has at least the active bits of another @p required FlagSet
     */
    bool isSuperSetOf(const FlagSet& required) const
    {
      return ((*this | required) == *this);
    }

    /**
     * @brief Check if this FlagSet has the bit for @p required
     */
    bool isSuperSetOf(const ENUM& required) const
    {
      return ((*this | required) == *this);
    }

    /// internal representation (mostly for illustrative purposes)
    UInt64 value() const
    {
      return value_;
    }

  private:
    /// computes pow(2, r)
    UInt64 getPow_(const ENUM& en) const
    {
      OPENMS_PRECONDITION((int)en >= 0, "Enum value is too small!")
      OPENMS_PRECONDITION((int)en <= 63, "Enum value is too large!")
      return UInt64(1) << UInt64(en);
    }
    UInt64 value_;
  };
}
