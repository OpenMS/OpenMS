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
      Int ret;

      // boost::spirit::qi was found to be vastly superior to boost::lexical_cast or stringstream extraction (especially for VisualStudio),
      // so don't change this unless you have benchmarks for all platforms!
      String::ConstIterator it = this_s.begin();
      if (!boost::spirit::qi::phrase_parse(it, this_s.end(), boost::spirit::qi::int_, boost::spirit::ascii::space, ret))
      {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Could not convert string '") + this_s + "' to an integer value");
      }
      // was the string parsed (white spaces are skipped automatically!) completely? If not, we have a problem because a previous split might have used the wrong split char
      if (it != this_s.end())
      {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Prefix of string '") + this_s + "' successfully converted to an int32 value. Additional characters found at position " + (int)(distance(this_s.begin(), it) + 1));
      }
      return ret;
    }

    static Int64 toInt64(const String& this_s)
    {
      Int64 ret;

      // boost::spirit::qi was found to be vastly superior to boost::lexical_cast or stringstream extraction (especially for VisualStudio),
      // so don't change this unless you have benchmarks for all platforms!
      String::ConstIterator it = this_s.begin();
      if (!boost::spirit::qi::phrase_parse(it, this_s.end(), boost::spirit::qi::long_long, boost::spirit::ascii::space, ret))
      {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Could not convert string '") + this_s + "' to an int64 value");
      }
      // was the string parsed (white spaces are skipped automatically!) completely? If not, we have a problem because a previous split might have used the wrong split char
      if (it != this_s.end())
      {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         String("Prefix of string '") + this_s + "' successfully converted to an integer value. Additional characters found at position " +
                                           (int)(distance(this_s.begin(), it) + 1));
      }
      return ret;
    }

    static float toFloat(const String& this_s)
    {
      float ret;

      // boost::spirit::qi was found to be vastly superior to boost::lexical_cast or stringstream extraction (especially for VisualStudio),
      // so don't change this unless you have benchmarks for all platforms!
      String::ConstIterator it = this_s.begin();
      if (!boost::spirit::qi::phrase_parse(it, this_s.end(), parse_float_, boost::spirit::ascii::space, ret))
      {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Could not convert string '") + this_s + "' to a float value");
      }
      // was the string parsed (white spaces are skipped automatically!) completely? If not, we have a problem because a previous split might have used the wrong split char
      if (it != this_s.end())
      {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Prefix of string '") + this_s + "' successfully converted to a float value. Additional characters found at position " + (int)(distance(this_s.begin(), it) + 1));
      }
      return ret;
    }

    /**
      @brief convert String (leading and trailing whitespace allowed) to double

      @p s Input string which represents a double, e.g. " 12.3 "
      @return A double representation of @p s
      @throws Exception::ConversionError if the string is not completely explained by the double (whitespaces are allowed)
    */
    static double toDouble(const String& s)
    {
      double ret;
      // boost::spirit::qi was found to be vastly superior to boost::lexical_cast or stringstream extraction (especially for VisualStudio),
      // so don't change this unless you have benchmarks for all platforms!
      String::ConstIterator it = s.begin();
      if (!boost::spirit::qi::phrase_parse(it, s.end(), parse_double_, boost::spirit::ascii::space, ret))
      {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Could not convert string '") + s + "' to a double value");
      }
      // was the string parsed (white spaces are skipped automatically!) completely? If not, we have a problem because a previous split might have used the wrong split char
      if (it != s.end())
      {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Prefix of string '") + s + "' successfully converted to a double value. Additional characters found at position " + (int)(distance(s.begin(), it) + 1));
      }
      return ret;
    }

    /// Reads a double from an iterator position.
    /// The begin iterator is modified (advanced) if parsing was successful.
    /// The @p target only contains a valid result if the functions returns true (i.e. parsing succeeded).
    /// Whitespaces before and after the double are NOT consumed!
    template <typename IteratorT>
    static bool extractDouble(IteratorT& begin, const IteratorT& end, double& target)
    {
      // boost::spirit::qi was found to be vastly superior to boost::lexical_cast or stringstream extraction (especially for VisualStudio),
      // so don't change this unless you have benchmarks for all platforms!

      // qi::parse() does not consume whitespace before or after the double (qi::parse_phrase() would).
      return boost::spirit::qi::parse(begin, end, parse_double_, target);
    }

  private:
  
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

  }
} // namespace OPENMS

