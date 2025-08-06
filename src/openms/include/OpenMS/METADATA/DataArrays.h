// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/MetaInfoDescription.h>
#include <OpenMS/CONCEPT/Helpers.h> // for ADL version of <=> for STL containers

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
      /// Less than operator
      bool operator<(const FloatDataArray& rhs) const
      {
        const MetaInfoDescription& lhs_meta = *this;
        const MetaInfoDescription& rhs_meta = rhs;
        if (lhs_meta != rhs_meta)
          return lhs_meta < rhs_meta;
        return static_cast<const std::vector<float>&>(*this) < static_cast<const std::vector<float>&>(rhs);
      }

      /// Less than or equal operator
      bool operator<=(const FloatDataArray& rhs) const
      {
        return *this < rhs || *this == rhs;
      }

      /// Greater than operator
      bool operator>(const FloatDataArray& rhs) const
      {
        return !(*this <= rhs);
      }

      /// Greater than or equal operator
      bool operator>=(const FloatDataArray& rhs) const
      {
        return !(*this < rhs);
      }

      /// Not equal operator
      bool operator!=(const FloatDataArray& rhs) const
      {
        return !(*this == rhs);
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
      /// Less than operator
      bool operator<(const IntegerDataArray& rhs) const
      {
        const MetaInfoDescription& lhs_meta = *this;
        const MetaInfoDescription& rhs_meta = rhs;
        if (lhs_meta != rhs_meta)
          return lhs_meta < rhs_meta;
        return static_cast<const std::vector<Int>&>(*this) < static_cast<const std::vector<Int>&>(rhs);
      }

      /// Less than or equal operator
      bool operator<=(const IntegerDataArray& rhs) const
      {
        return *this < rhs || *this == rhs;
      }

      /// Greater than operator
      bool operator>(const IntegerDataArray& rhs) const
      {
        return !(*this <= rhs);
      }

      /// Greater than or equal operator
      bool operator>=(const IntegerDataArray& rhs) const
      {
        return !(*this < rhs);
      }

      /// Not equal operator
      bool operator!=(const IntegerDataArray& rhs) const
      {
        return !(*this == rhs);
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
      /// Less than operator
      bool operator<(const StringDataArray& rhs) const
      {
        const MetaInfoDescription& lhs_meta = *this;
        const MetaInfoDescription& rhs_meta = rhs;
        if (lhs_meta != rhs_meta)
          return lhs_meta < rhs_meta;
        return static_cast<const std::vector<String>&>(*this) < static_cast<const std::vector<String>&>(rhs);
      }

      /// Less than or equal operator
      bool operator<=(const StringDataArray& rhs) const
      {
        return *this < rhs || *this == rhs;
      }

      /// Greater than operator
      bool operator>(const StringDataArray& rhs) const
      {
        return !(*this <= rhs);
      }

      /// Greater than or equal operator
      bool operator>=(const StringDataArray& rhs) const
      {
        return !(*this < rhs);
      }

      /// Not equal operator
      bool operator!=(const StringDataArray& rhs) const
      {
        return !(*this == rhs);
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

