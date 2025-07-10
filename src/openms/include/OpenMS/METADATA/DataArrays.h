// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/MetaInfoDescription.h>
#include <compare>

namespace OpenMS
{
  namespace DataArrays
  {

    /// Float data array class
    class FloatDataArray :
      public MetaInfoDescription,
      public std::vector<float>
    {
      using std::vector<float>::vector; // to allow for aggregate initialization of FloatDataArray

    public:
      /// Spaceship operator for three-way comparison
      auto operator<=>(const FloatDataArray& rhs) const
      {
        // Compare MetaInfoDescription base class first
        const MetaInfoDescription& lhs_meta = *this;
        const MetaInfoDescription& rhs_meta = rhs;
        if (auto cmp = lhs_meta <=> rhs_meta; cmp != 0)
          return std::partial_ordering(cmp);
        
        // Then compare the vector part (returns partial_ordering for floats)
        return static_cast<const std::vector<float>&>(*this) <=> static_cast<const std::vector<float>&>(rhs);
      }

      /// Equality operator
      bool operator==(const FloatDataArray& rhs) const
      {
        return static_cast<const MetaInfoDescription&>(*this) == static_cast<const MetaInfoDescription&>(rhs) && 
               static_cast<const std::vector<float>&>(*this) == static_cast<const std::vector<float>&>(rhs);
      }
    };

    /// Integer data array class
    class IntegerDataArray :
      public MetaInfoDescription,
      public std::vector<Int>
    {
      using std::vector<int>::vector; // to allow for aggregate initialization of IntegerDataArray

    public:
      /// Spaceship operator for three-way comparison
      auto operator<=>(const IntegerDataArray& rhs) const
      {
        // Compare MetaInfoDescription base class first
        const MetaInfoDescription& lhs_meta = *this;
        const MetaInfoDescription& rhs_meta = rhs;
        if (auto cmp = lhs_meta <=> rhs_meta; cmp != 0)
          return cmp;
        
        // Then compare the vector part
        return static_cast<const std::vector<Int>&>(*this) <=> static_cast<const std::vector<Int>&>(rhs);
      }

      /// Equality operator
      bool operator==(const IntegerDataArray& rhs) const
      {
        return static_cast<const MetaInfoDescription&>(*this) == static_cast<const MetaInfoDescription&>(rhs) && 
               static_cast<const std::vector<Int>&>(*this) == static_cast<const std::vector<Int>&>(rhs);
      }
    };

    /// String data array class
    class StringDataArray :
      public MetaInfoDescription,
      public std::vector<String>
    {
      using std::vector<String>::vector; // to allow for aggregate initialization of StringDataArray

    public:
      /// Spaceship operator for three-way comparison
      auto operator<=>(const StringDataArray& rhs) const
      {
        // Compare MetaInfoDescription base class first
        const MetaInfoDescription& lhs_meta = *this;
        const MetaInfoDescription& rhs_meta = rhs;
        if (auto cmp = lhs_meta <=> rhs_meta; cmp != 0)
          return cmp;
        
        // Then compare the vector part
        return static_cast<const std::vector<String>&>(*this) <=> static_cast<const std::vector<String>&>(rhs);
      }

      /// Equality operator
      bool operator==(const StringDataArray& rhs) const
      {
        return static_cast<const MetaInfoDescription&>(*this) == static_cast<const MetaInfoDescription&>(rhs) && 
               static_cast<const std::vector<String>&>(*this) == static_cast<const std::vector<String>&>(rhs);
      }
    };

  }
} // namespace OpenMS

