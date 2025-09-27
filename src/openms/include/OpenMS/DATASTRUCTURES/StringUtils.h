// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg, Chris Bielow $
// $Authors: Marc Sturm, Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/StringUtilsSimple.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/CONCEPT/PrecisionWrapper.h>

#include <QtCore/QString>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/type_traits.hpp>

#include <string>
#include <vector>
#include <charconv>
#include <cctype>
#include <cstring>


namespace OpenMS
{
  class String;

  class OPENMS_DLLAPI StringUtilsHelper
  {

public:

    //
    /// Functions
    //
    static Int toInt32(const String& this_s)
    {
      // std::from_chars provides better performance and locale independence compared to boost::spirit::qi
      return parseInteger_<Int>(this_s);
    }

    static Int64 toInt64(const String& this_s)
    {
      // std::from_chars provides better performance and locale independence compared to boost::spirit::qi
      return parseInteger_<Int64>(this_s);
    }

    static float toFloat(const String& this_s)
    {
      // std::from_chars provides better performance and locale independence compared to boost::spirit::qi
      return parseFloatingPoint_<float>(this_s);
    }

    /**
      @brief convert String (leading and trailing whitespace allowed) to double

      @p s Input string which represents a double, e.g. " 12.3 "
      @return A double representation of @p s
      @throws Exception::ConversionError if the string is not completely explained by the double (whitespaces are allowed)
    */
    static double toDouble(const String& s)
    {
      // std::from_chars provides better performance and locale independence compared to boost::spirit::qi
      return parseFloatingPoint_<double>(s);
    }

    /// Reads a double from an iterator position.
    /// The begin iterator is modified (advanced) if parsing was successful.
    /// The @p target only contains a valid result if the functions returns true (i.e. parsing succeeded).
    /// Whitespaces before and after the double are NOT consumed!
    template <typename IteratorT>
    static bool extractDouble(IteratorT& begin, const IteratorT& end, double& target)
    {
      // std::from_chars provides better performance and locale independence compared to boost::spirit::qi
      
      // For std::from_chars we need contiguous char data
      // This assumes the iterators are pointing to contiguous char data (like std::string iterators)
      if (begin == end) {
        return false;
      }
      
      const char* char_begin = &(*begin);
      const char* char_end = char_begin + std::distance(begin, end);
      
      // Check for NaN first (to maintain compatibility with boost::spirit::qi behavior)
      if (isNaNString_(char_begin, char_end)) {
        target = std::numeric_limits<double>::quiet_NaN();
        begin = end; // consume entire string for NaN
        return true;
      }
      
      auto [ptr, ec] = std::from_chars(char_begin, char_end, target);
      
      if (ec == std::errc{}) {
        // Advance the iterator by the number of characters consumed
        std::advance(begin, ptr - char_begin);
        return true;
      }
      
      return false;
    }

    /// Reads an int from an iterator position.
    /// The begin iterator is modified (advanced) if parsing was successful.
    /// The @p target only contains a valid result if the functions returns true (i.e. parsing succeeded).
    /// Whitespaces before and after the int are NOT consumed!
    template <typename IteratorT>
    static bool extractInt(IteratorT& begin, const IteratorT& end, int& target)
    {
      // std::from_chars provides better performance and locale independence compared to boost::spirit::qi
      
      // For std::from_chars we need contiguous char data
      // This assumes the iterators are pointing to contiguous char data (like std::string iterators)
      if (begin == end) {
        return false;
      }
      
      const char* char_begin = &(*begin);
      const char* char_end = char_begin + std::distance(begin, end);
      
      auto [ptr, ec] = std::from_chars(char_begin, char_end, target);
      
      if (ec == std::errc{}) {
        // Advance the iterator by the number of characters consumed
        std::advance(begin, ptr - char_begin);
        return true;
      }
      
      return false;
    }

  private:
  
    // Helper functions for std::from_chars-based conversion
    // These provide higher performance and locale independence compared to boost::spirit::qi
    
    /// Skip whitespace from the beginning of a range
    static const char* skipWhitespace_(const char* begin, const char* end)
    {
      while (begin != end && std::isspace(*begin)) {
        ++begin;
      }
      return begin;
    }
    
    /// Skip whitespace from the end of a range (reverse direction)
    static const char* skipTrailingWhitespace_(const char* begin, const char* end)
    {
      while (end != begin && std::isspace(*(end - 1))) {
        --end;
      }
      return end;
    }
    
    /// Check if string (after trimming whitespace) represents NaN
    static bool isNaNString_(const char* begin, const char* end)
    {
      // Trim whitespace
      begin = skipWhitespace_(begin, end);
      end = skipTrailingWhitespace_(begin, end);
      
      const size_t len = end - begin;
      if (len == 3) {
        return (strncmp(begin, "nan", 3) == 0 || strncmp(begin, "NaN", 3) == 0);
      }
      return false;
    }
    
    /// Parse floating point using std::from_chars with proper NaN and whitespace handling
    template<typename T>
    static T parseFloatingPoint_(const String& s)
    {
      const char* begin = s.data();
      const char* end = s.data() + s.size();
      
      // Check for NaN first (to maintain compatibility with boost::spirit::qi behavior)
      if (isNaNString_(begin, end)) {
        return std::numeric_limits<T>::quiet_NaN();
      }
      
      // Trim leading whitespace
      begin = skipWhitespace_(begin, end);
      if (begin == end) {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
          String("Could not convert string '") + s + "' to a floating point value");
      }
      
      T result;
      auto [ptr, ec] = std::from_chars(begin, end, result);
      
      if (ec != std::errc{}) {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
          String("Could not convert string '") + s + "' to a floating point value");
      }
      
      // Skip trailing whitespace
      ptr = skipWhitespace_(ptr, end);
      
      // Check if we consumed the entire (trimmed) string
      if (ptr != end) {
        String error_msg = String("Prefix of string '") + s + "' successfully converted to a floating point value. Additional characters found at position " + String((int)(ptr - s.data() + 1));
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, error_msg);
      }
      
      return result;
    }
    
    /// Parse integer using std::from_chars with proper whitespace handling
    template<typename T>
    static T parseInteger_(const String& s)
    {
      const char* begin = s.data();
      const char* end = s.data() + s.size();
      
      // Trim leading whitespace
      begin = skipWhitespace_(begin, end);
      if (begin == end) {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
          String("Could not convert string '") + s + "' to an integer value");
      }
      
      T result;
      auto [ptr, ec] = std::from_chars(begin, end, result);
      
      if (ec != std::errc{}) {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
          String("Could not convert string '") + s + "' to an integer value");
      }
      
      // Skip trailing whitespace
      ptr = skipWhitespace_(ptr, end);
      
      // Check if we consumed the entire (trimmed) string
      if (ptr != end) {
        String error_msg = String("Prefix of string '") + s + "' successfully converted to an integer value. Additional characters found at position " + String((int)(ptr - s.data() + 1));
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, error_msg);
      }
      
      return result;
    }
  
    /*
      @brief A fixed Boost:pi real parser policy, capable of dealing with 'nan' without crashing
  
      The original Boost implementation has a bug, see https://svn.boost.org/trac/boost/ticket/6955.
      Can be removed if Boost 1.60 or above is required
      
    */
    template <typename T>
    struct real_policies_NANfixed_ : boost::spirit::qi::real_policies<T>
    {
      template <typename Iterator, typename Attribute>
      static bool
        parse_nan(Iterator& first, Iterator const& last, Attribute& attr_)
      {
        if (first == last)
          return false;   // end of input reached
  
        if (*first != 'n' && *first != 'N')
          return false;   // not "nan"
  
        // nan[(...)] ?
        if (boost::spirit::qi::detail::string_parse("nan", "NAN", first, last, boost::spirit::qi::unused))
        {
          if (first != last && *first == '(')  /* this check is broken in boost 1.49 - (at least) 1.54; fixed in 1.60 */
          {
            // skip trailing (...) part
            Iterator i = first;
  
            while (++i != last && *i != ')')
              ;
            if (i == last)
              return false;     // no trailing ')' found, give up
  
            first = ++i;
          }
          attr_ = std::numeric_limits<T>::quiet_NaN();
          return true;
        }
        return false;
      }
    };
    
    // Qi parsers using the 'real_policies_NANfixed_' template which allows for 'nan'
    // (the original Boost implementation has a bug, see https://svn.boost.org/trac/boost/ticket/6955)
    static boost::spirit::qi::real_parser<double, real_policies_NANfixed_<double> > parse_double_;
    static boost::spirit::qi::real_parser<float, real_policies_NANfixed_<float> > parse_float_;
    static boost::spirit::qi::int_parser<> parse_int_;

  };

  namespace StringUtils 
  {

    [[maybe_unused]] static String number(double d, UInt n)
    {
      return QString::number(d, 'f', n);
    }

    [[maybe_unused]] static QString toQString(const String & this_s) 
    {
      return QString(this_s.c_str());
    }

    [[maybe_unused]] static Int32 toInt32(const String & this_s)
    {
      return StringUtilsHelper::toInt32(this_s);
    }

    [[maybe_unused]] static Int64 toInt64(const String& this_s)
    {
      return StringUtilsHelper::toInt64(this_s);
    }

    [[maybe_unused]] static float toFloat(const String & this_s)
    {
      return StringUtilsHelper::toFloat(this_s);
    }

    [[maybe_unused]] static double toDouble(const String & this_s)
    {
      return StringUtilsHelper::toDouble(this_s);
    }

    template <typename IteratorT>
    static bool extractDouble(IteratorT& begin, const IteratorT& end, double& target)
    {
      return StringUtilsHelper::extractDouble(begin, end, target);
    }

    template <typename IteratorT>
    static bool extractInt(IteratorT& begin, const IteratorT& end, int& target)
    {
      return StringUtilsHelper::extractInt(begin, end, target);
    }
  }
} // namespace OPENMS

