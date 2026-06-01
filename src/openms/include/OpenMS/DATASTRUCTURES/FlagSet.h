// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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

    /// checks if any bit is set
    bool empty() const
    {
      return value_ == 0;
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
