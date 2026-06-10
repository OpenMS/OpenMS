// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/OpenMSConfig.h>
#include <OpenMS/config.h>

#include <cmath>
#include <iterator>
#include <vector>
#include <algorithm>

namespace OpenMS
{

  /**
    @brief Vector of signed integers.

    @ingroup Datastructures
  */
  typedef std::vector<Int> IntList;

  /**
   @brief Vector of double precision real types.

   @ingroup Datastructures
   */
  typedef std::vector<double> DoubleList;


  /**
   @brief Vector of String.

   @ingroup Datastructures
  */
  typedef std::vector<std::string> StringList;

  /**
    @brief Collection of utility functions for management of vectors.

    @ingroup Datastructures
  */

  class OPENMS_DLLAPI ListUtils
  {
private:
    /**
      @brief Predicate to check double equality with a given tolerance.
    */
    struct DoubleTolerancePredicate_
    {
      DoubleTolerancePredicate_(const double& target, const double& tolerance) :
        tolerance_(tolerance),
        target_(target)
      {}

      /**
        @brief Returns true if \| @p value - @p target \| \< @p tolerance.

        @param[in] value The value to test.
        @return true if \| @p value - @p target \| \< @p tolerance, false otherwise.
      */
      inline bool operator()(const double& value) const
      {
        return std::fabs(value - target_) < tolerance_;
      }

private:
      /// The allowed tolerance.
      double tolerance_;
      /// The target value that should be found.
      double target_;
    };

public:
    /**
      @brief Returns a list that is created by splitting the given comma-separated string.
      @note If converted to vector<std::string> the strings are not trimmed.
      @note The values get converted by boost::lexical_cast so a valid conversion from std::string to T needs to be available.

      @param[in] str The string that should be split and converted to a list.
      @param[in] splitter The separator to look for in @p str
      @return A vector containing the elements of the string converted into type T.
    */
    template <typename T>
    static std::vector<T> create(const std::string& str, const char splitter = ',')
    {
      // temporary storage for the individual elements of the string
      std::vector<std::string> temp_string_vec;
      StringUtils::split(str, splitter, temp_string_vec);
      return create<T>(temp_string_vec);
    }

    /**
      @brief Converts a vector of strings to a vector of the target type T.
      @note The strings are trimmed before conversion.
      @note The values get converted by boost::lexical_cast so a valid conversion from std::string to T needs to be available.

      @param[in] s The vector of strings that should be converted.
      @return A vector containing the elements of input vector converted into type T.
    */
    template <typename T>
    static std::vector<T> create(const std::vector<std::string>& s);


    /**
  @brief Converts a vector of T's to a vector of Strings.

  @param[in] s The vector of T's that should be converted.
  @return A vector containing the elements of input vector converted into Strings.
*/
    template <typename T>
    static std::vector<std::string> toStringList(const std::vector<T>& s)
    {
      StringList out;
      out.reserve(s.size());
      for (const auto& elem : s) out.push_back(StringUtils::toStr(elem));
      return out;
    }

    /**
      @brief Checks whether the element @p elem is contained in the given container.

      @param[in] container The container to check.
      @param[in] elem The element to check whether it is in the container or not.

      @return True if @p elem is contained in @p container, false otherwise.
    */
    template <typename T, typename E>
    static bool contains(const std::vector<T>& container, const E& elem)
    {
      return find(container.begin(), container.end(), elem) != container.end();
    }

    /**
      @brief Checks whether the element @p elem is contained in the given container of floating point numbers.

      @param[in] container The container of doubles to check.
      @param[in] elem The element to check whether it is in the container or not.
      @param[in] tolerance The allowed tolerance for the double.

      @return True if @p elem is contained in @p container, false otherwise.
    */
    static bool contains(const std::vector<double>& container, const double& elem, double tolerance = 0.00001)
    {
      return find_if(container.begin(), container.end(), DoubleTolerancePredicate_(elem, tolerance)) != container.end();
    }


    enum class CASE { SENSITIVE, INSENSITIVE};
    /**
    @brief Checks whether the String @p elem is contained in the given container (potentially case insensitive)

    @param[in] container The container of std::string to check.
    @param[in] elem The element to check whether it is in the container or not.
    @param[in] case_sensitive Do the comparison case sensitive or insensitive

    @return True if @p elem is contained in @p container, false otherwise.
    */
    static bool contains(const std::vector<std::string>& container, std::string elem, const CASE case_sensitive)
    {
      if (case_sensitive == CASE::SENSITIVE) return contains(container, elem);
      // case insensitive ...
      StringUtils::toLower(elem);
      return find_if(container.begin(), container.end(), [&elem](std::string ce) {
        return elem == StringUtils::toLower(ce);
      }) != container.end();
    }

    /**
      @brief Concatenates all elements of the @p container and puts the @p glue string between elements.

      @param[in] container The container to concatenate;
      @param[in] glue The string to add in between elements.
    */
    template <typename T>
    static std::string concatenate(const std::vector<T>& container, const std::string& glue = "")
    {
      return concatenate< std::vector<T> >(container, glue);
    }

    /**
      @brief Concatenates all elements of the @p container and puts the @p glue string between elements.

      @param[in] container The container to concatenate; must have begin() and end() iterator.
      @param[in] glue The string to add in between elements.
    */
    template <typename T>
    static std::string concatenate(const T& container, const std::string& glue = "")
    {
      // handle empty containers
      if (container.empty()) return "";

      typename T::const_iterator it = container.begin();
      std::string ret =StringUtils::toStr(*it);
      // we have handled the first element
      ++it;
      // add the rest
      for (; it != container.end(); ++it)
      {
        ret += (glue + StringUtils::toStr(*it));
      }

      return ret;
    }

    /**
       @brief Get the index of the first occurrence of an element in the vector (or -1 if not found)
    */
    template <typename T, typename E>
    static Int getIndex(const std::vector<T>& container, const E& elem)
    {
      typename std::vector<T>::const_iterator pos =
        std::find(container.begin(), container.end(), elem);
      if (pos == container.end()) return -1;

      return static_cast<Int>(std::distance(container.begin(), pos));
    }

  };

  namespace detail
  {
    template <typename T>
    T convert(const std::string& s);
  
    template<>
    inline Int32 convert(const std::string& s)
    {
      return StringUtils::toInt32(s);
    }
    template<>
    inline double convert(const std::string& s)
    {
      return StringUtils::toDouble(s);
    }
    template<>
    inline float convert(const std::string& s)
    {
      return StringUtils::toFloat(s);
    }
    template<>
    inline std::string convert(const std::string& s)
    {
        return static_cast<std::string>(s);
    }
  }

  template <typename T>
  inline std::vector<T> ListUtils::create(const std::vector<std::string>& s)
  {
    std::vector<T> c;
    c.reserve(s.size());
    for (std::vector<std::string>::const_iterator it = s.begin(); it != s.end(); ++it)
    {
      try
      {
        c.push_back(detail::convert<T>(StringUtils::trimmed(*it))); // succeeds only if the whole output can be explained, i.e. "1.3 3" will fail (which is good)
      }
      catch (...)
      {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,std::string("Could not convert string '") + *it + "'");
      }
    }

    return c;
  }

  /// create specialization for std::string since we do not need to cast here
  template <>
  inline std::vector<std::string> ListUtils::create(const std::vector<std::string>& s)
  {
    return s;
  }

} // namespace OpenMS
