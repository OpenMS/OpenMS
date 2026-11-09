// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <algorithm>
#include <string>
#include <vector>

namespace OpenMS
{
  namespace Helpers
  {
    /// return the index of an element in a container
    /// useful for matching 'names_of_...' arrays to their enum value
    template <class ContainerType>
    Size indexOf(const ContainerType& cont, const typename ContainerType::value_type& val)
    {
      auto it = std::find(cont.begin(), cont.end(), val);
      if (it == cont.end())
      {
        throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, val);
      }
      return std::distance(cont.begin(), it);
    }

    /**
      @brief Returns the names of all enum values, i.e. the first @p size_value entries of @p names.

      Implements the <tt>getAllNamesOfX()</tt> functions of enums with a <tt>NamesOfX[]</tt> array and a
      trailing <tt>SIZE_OF_X</tt> value.

      @param[in] names Array of names, indexed by enum value
      @param[in] size_value The enum's <tt>SIZE_OF_X</tt> value (number of valid enum values)
    */
    template <class EnumType>
    std::vector<std::string> enumNames(const std::string* names, EnumType size_value)
    {
      std::vector<std::string> result;
      result.reserve(static_cast<size_t>(size_value));
      for (size_t i = 0; i < static_cast<size_t>(size_value); ++i)
      {
        result.push_back(names[i]);
      }
      return result;
    }

    /**
      @brief Converts an enum value to its entry in @p names.

      Implements the <tt>xToString(X)</tt> functions of enums with a <tt>NamesOfX[]</tt> array and a
      trailing <tt>SIZE_OF_X</tt> value.

      @param[in] names Array of names, indexed by enum value
      @param[in] value The enum value to convert
      @param[in] size_value The enum's <tt>SIZE_OF_X</tt> value, which is rejected
      @param[in] size_name Name of @p size_value, reported as the invalid value (e.g. "SIZE_OF_X")
      @param[in] function Name of the calling function, reported by the exception

      @exception Exception::InvalidValue ("Value not allowed") if @p value is @p size_value
    */
    template <class EnumType>
    const std::string& enumToName(const std::string* names, EnumType value, EnumType size_value, const char* size_name, const char* function)
    {
      if (value == size_value)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, function, "Value not allowed", size_name);
      }
      return names[static_cast<size_t>(value)];
    }

    /**
      @brief Converts a name to the enum value at its index in @p names.

      Implements the <tt>toX(const std::string&)</tt> functions of enums with a <tt>NamesOfX[]</tt> array
      and a trailing <tt>SIZE_OF_X</tt> value. Unlike indexOf(), this throws Exception::InvalidValue.

      @param[in] names Array of names, indexed by enum value
      @param[in] size_value The enum's <tt>SIZE_OF_X</tt> value (number of valid enum values)
      @param[in] name The name to convert
      @param[in] message Message of the exception thrown if @p name is unknown (e.g. "Value unknown")
      @param[in] function Name of the calling function, reported by the exception

      @exception Exception::InvalidValue if @p name is not among the first @p size_value entries of @p names
    */
    template <class EnumType>
    EnumType enumFromName(const std::string* names, EnumType size_value, const std::string& name, const char* message, const char* function)
    {
      const std::string* first = names;
      const std::string* last = names + static_cast<size_t>(size_value);
      const auto it = std::find(first, last, name);
      if (it == last)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, function, message, name);
      }
      return static_cast<EnumType>(it - first);
    }

  }
}


