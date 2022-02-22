// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
    static Int toInt(const String & this_s)
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
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Prefix of string '") + this_s + "' successfully converted to an integer value. Additional characters found at position " + (int)(distance(this_s.begin(), it) + 1));
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

    [[maybe_unused]] static Int toInt(const String & this_s)
    {
      return StringUtilsHelper::toInt(this_s);
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

