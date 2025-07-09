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
        if (lhs_meta != rhs_meta)
        {
          // Since MetaInfoDescription doesn't have spaceship, we provide a consistent ordering
          // by comparing names first, then other fields deterministically
          if (auto cmp = lhs_meta.getName() <=> rhs_meta.getName(); cmp != 0)
            return std::partial_ordering(cmp);
          
          // If names are equal, compare data processing vectors by size first
          const auto& lhs_dp = lhs_meta.getDataProcessing();
          const auto& rhs_dp = rhs_meta.getDataProcessing();
          if (auto cmp = lhs_dp.size() <=> rhs_dp.size(); cmp != 0)
            return std::partial_ordering(cmp);
          
          // Compare MetaInfoInterface by getting keys and comparing them
          const MetaInfoInterface& lhs_iface = static_cast<const MetaInfoInterface&>(lhs_meta);
          const MetaInfoInterface& rhs_iface = static_cast<const MetaInfoInterface&>(rhs_meta);
          std::vector<UInt> lhs_keys, rhs_keys;
          lhs_iface.getKeys(lhs_keys);
          rhs_iface.getKeys(rhs_keys);
          
          // Compare the keys first
          if (auto cmp = lhs_keys <=> rhs_keys; cmp != 0)
            return std::partial_ordering(cmp);
          
          // If keys are same, compare values for each key
          for (UInt key : lhs_keys)
          {
            const DataValue& lhs_val = lhs_iface.getMetaValue(key);
            const DataValue& rhs_val = rhs_iface.getMetaValue(key);
            if (lhs_val != rhs_val)
            {
              // DataValue doesn't have spaceship, so we compare string representations
              String lhs_str = lhs_val.toString();
              String rhs_str = rhs_val.toString();
              if (auto cmp = lhs_str <=> rhs_str; cmp != 0)
                return std::partial_ordering(cmp);
            }
          }
          
          // If all comparable fields are equal, consider them equivalent for ordering
          return std::partial_ordering::equivalent;
        }
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
        if (lhs_meta != rhs_meta)
        {
          // Since MetaInfoDescription doesn't have spaceship, we provide a consistent ordering
          // by comparing names first, then other fields deterministically
          if (auto cmp = lhs_meta.getName() <=> rhs_meta.getName(); cmp != 0)
            return cmp;
          
          // If names are equal, compare data processing vectors by size first
          const auto& lhs_dp = lhs_meta.getDataProcessing();
          const auto& rhs_dp = rhs_meta.getDataProcessing();
          if (auto cmp = lhs_dp.size() <=> rhs_dp.size(); cmp != 0)
            return cmp;
          
          // Compare MetaInfoInterface by getting keys and comparing them
          const MetaInfoInterface& lhs_iface = static_cast<const MetaInfoInterface&>(lhs_meta);
          const MetaInfoInterface& rhs_iface = static_cast<const MetaInfoInterface&>(rhs_meta);
          std::vector<UInt> lhs_keys, rhs_keys;
          lhs_iface.getKeys(lhs_keys);
          rhs_iface.getKeys(rhs_keys);
          
          // Compare the keys first
          if (auto cmp = lhs_keys <=> rhs_keys; cmp != 0)
            return cmp;
          
          // If keys are same, compare values for each key
          for (UInt key : lhs_keys)
          {
            const DataValue& lhs_val = lhs_iface.getMetaValue(key);
            const DataValue& rhs_val = rhs_iface.getMetaValue(key);
            if (lhs_val != rhs_val)
            {
              // DataValue doesn't have spaceship, so we compare string representations
              String lhs_str = lhs_val.toString();
              String rhs_str = rhs_val.toString();
              if (auto cmp = lhs_str <=> rhs_str; cmp != 0)
                return cmp;
            }
          }
          
          // If all comparable fields are equal, consider them equivalent for ordering
          return std::strong_ordering::equivalent;
        }
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
        if (lhs_meta != rhs_meta)
        {
          // Since MetaInfoDescription doesn't have spaceship, we provide a consistent ordering
          // by comparing names first, then other fields deterministically
          if (auto cmp = lhs_meta.getName() <=> rhs_meta.getName(); cmp != 0)
            return cmp;
          
          // If names are equal, compare data processing vectors by size first
          const auto& lhs_dp = lhs_meta.getDataProcessing();
          const auto& rhs_dp = rhs_meta.getDataProcessing();
          if (auto cmp = lhs_dp.size() <=> rhs_dp.size(); cmp != 0)
            return cmp;
          
          // Compare MetaInfoInterface by getting keys and comparing them
          const MetaInfoInterface& lhs_iface = static_cast<const MetaInfoInterface&>(lhs_meta);
          const MetaInfoInterface& rhs_iface = static_cast<const MetaInfoInterface&>(rhs_meta);
          std::vector<UInt> lhs_keys, rhs_keys;
          lhs_iface.getKeys(lhs_keys);
          rhs_iface.getKeys(rhs_keys);
          
          // Compare the keys first
          if (auto cmp = lhs_keys <=> rhs_keys; cmp != 0)
            return cmp;
          
          // If keys are same, compare values for each key
          for (UInt key : lhs_keys)
          {
            const DataValue& lhs_val = lhs_iface.getMetaValue(key);
            const DataValue& rhs_val = rhs_iface.getMetaValue(key);
            if (lhs_val != rhs_val)
            {
              // DataValue doesn't have spaceship, so we compare string representations
              String lhs_str = lhs_val.toString();
              String rhs_str = rhs_val.toString();
              if (auto cmp = lhs_str <=> rhs_str; cmp != 0)
                return cmp;
            }
          }
          
          // If all comparable fields are equal, consider them equivalent for ordering
          return std::strong_ordering::equivalent;
        }
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

