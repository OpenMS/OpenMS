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
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/CONCEPT/PrecisionWrapper.h>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/type_traits.hpp>

#include <boost/spirit/home/karma/numeric/detail/real_utils.hpp>

#include <string>
#include <vector>


namespace OpenMS
{
  class String;

  namespace StringConversions
  {

    namespace Detail
    {
      // Karma full precision float policy
      template<typename T>
      class BK_PrecPolicyFull : public boost::spirit::karma::real_policies<T>
      {
        typedef boost::spirit::karma::real_policies<T> base_policy_type;

      public:
        static unsigned precision(T /*n*/)
        {
          /* The following would be the only way for a lossless double-string-double
          * roundtrip but:
          * a) We only care about speed
          * b) Many tests have to be changed
          * c) In the end boost::karma is bugged and hard limits the fractional digits
          *    even though you have leading zeros (basically forcing scientific notation)
          *    for full precision https://github.com/boostorg/spirit/issues/585
          if (BK_PrecPolicyFull::floatfield(n))
          {
              T abs_n = boost::spirit::traits::get_absolute_value(n);
              if (abs_n >= 1)
              {
                  return std::numeric_limits<T>::max_digits10 - (floor(log10(abs_n)) + 1);
              }
              else
              {
                  return std::numeric_limits<T>::max_digits10 - (floor(log10(abs_n)));
              }
          }
          else
          {
              return std::numeric_limits<T>::max_digits10 - 1;
          }
          */
          return writtenDigits<T>();
        }

        //  we want the numbers always to be in scientific format
        static unsigned floatfield(T n)
        {
          if (boost::spirit::traits::test_zero(n)) return base_policy_type::fmtflags::fixed;

          T abs_n = boost::spirit::traits::get_absolute_value(n);
          // this is due to a bug in downstream thirdparty tools that only can read
          // up to 19 digits. https://github.com/OpenMS/OpenMS/issues/4627
          return (abs_n >= 1e4 || abs_n < 1e-2) ? base_policy_type::fmtflags::scientific : base_policy_type::fmtflags::fixed;
        }

        // we need this special 'NaN' since the default 'nan' is not recognized by downstream tools such as any Java-based tool (e.g. KNIME) trying to
        // parse our output files
        template<typename CharEncoding, typename Tag, typename OutputIterator>
        static bool nan(OutputIterator& sink, T n, bool force_sign)
        {
          return boost::spirit::karma::sign_inserter::call(sink, false, boost::spirit::traits::test_negative(n), force_sign)
                 && boost::spirit::karma::string_inserter<CharEncoding, Tag>::call(sink, "NaN");
        }
      };

    // Karma default (3-digits) precision float policy
      template<typename T>
      class BK_PrecPolicyShort : public boost::spirit::karma::real_policies<T>
      {
        typedef boost::spirit::karma::real_policies<T> base_policy_type;

      public:
        // we need this special 'NaN' since the default 'nan' is not recognized by downstream tools such as any Java-based tool (e.g. KNIME) trying to
        // parse our output files
        template<typename CharEncoding, typename Tag, typename OutputIterator>
        static bool nan(OutputIterator& sink, T n, bool force_sign)
        {
          return boost::spirit::karma::sign_inserter::call(sink, false, boost::spirit::traits::test_negative(n), force_sign)
                 && boost::spirit::karma::string_inserter<CharEncoding, Tag>::call(sink, "NaN");
        }
      };

      using BK_PrecPolicyFloatFull_type = boost::spirit::karma::real_generator<float, BK_PrecPolicyFull<float>>;
      const BK_PrecPolicyFloatFull_type BK_PrecPolicyFloatFull;
      using BK_PrecPolicyDoubleFull_type = boost::spirit::karma::real_generator<double, BK_PrecPolicyFull<double>>;
      const BK_PrecPolicyDoubleFull_type BK_PrecPolicyDoubleFull;
      using BK_PrecPolicyLongDoubleFull_type = boost::spirit::karma::real_generator<long double, BK_PrecPolicyFull<long double>>;
      const BK_PrecPolicyLongDoubleFull_type BK_PrecPolicyLongDoubleFull;

      using BK_PrecPolicyFloatShort_type = boost::spirit::karma::real_generator<float, BK_PrecPolicyShort<float>>;
      const BK_PrecPolicyFloatShort_type BK_PrecPolicyFloatShort;
      using BK_PrecPolicyDoubleShort_type = boost::spirit::karma::real_generator<double, BK_PrecPolicyShort<double>>;
      const BK_PrecPolicyDoubleShort_type BK_PrecPolicyDoubleShort;
      using BK_PrecPolicyLongDoubleShort_type = boost::spirit::karma::real_generator<long double, BK_PrecPolicyShort<long double>>;
      const BK_PrecPolicyLongDoubleShort_type BK_PrecPolicyLongDoubleShort;

    } // namespace Detail



    // toString functions (single argument)

    /// fallback template for general purpose using Boost::Karma; more specializations below
    /// does NOT clear the input string @p target, so appending is as efficient as possible
    template <typename T>
    inline void append(const T& i, String& target)
    {
      std::back_insert_iterator<std::string> sink(target);
      boost::spirit::karma::generate(sink, i);
    }

    /// fallback template for general purpose using Boost::Karma; more specializations below
    template <typename T>
    inline String toString(const T& i)
    {
      //std::stringstream s;
      //s << i;
      //return s.str();
      String str;
      append(i, str);
      return str;
    }
    

    /// low precision (3 fractional digits) conversion to string (Karma default)
    /// does NOT clear the input string @p target, so appending is as efficient as possible
    inline void appendLowP(float f, String& target)
    {
      std::back_insert_iterator<std::string> sink(target);
      boost::spirit::karma::generate(sink, Detail::BK_PrecPolicyFloatShort, f);
    }
    /// low precision (3 fractional digits) conversion to string (Karma default)
    inline String toStringLowP(float f)
    {
      String str;
      appendLowP(f, str);
      return str;
    }


    /// low precision (3 fractional digits) conversion to string (Karma default)
    /// does NOT clear the input string @p target, so appending is as efficient as possible
    inline void appendLowP(double d, String& target)
    {
      std::back_insert_iterator<std::string> sink(target);
      boost::spirit::karma::generate(sink, Detail::BK_PrecPolicyDoubleShort, d);
    }
    /// low precision (3 fractional digits) conversion to string (Karma default)
    inline String toStringLowP(double d)
    {
      String str;
      appendLowP(d, str);
      return str;
    }


    /// low precision (3 fractional digits) conversion to string (Karma default)
    inline void appendLowP(long double ld, String& target)
    {
      std::back_insert_iterator<std::string> sink(target);
      boost::spirit::karma::generate(sink, Detail::BK_PrecPolicyLongDoubleShort, ld);
    }
    /// low precision (3 fractional digits) conversion to string (Karma default)
    inline String toStringLowP(long double ld)
    {
      String str;
      appendLowP(ld, str);
      return str;
    }



    /// high precision (6 fractional digits) conversion to String
    inline void append(float f, String& target)
    {
      std::back_insert_iterator<std::string> sink(target);
      boost::spirit::karma::generate(sink, Detail::BK_PrecPolicyFloatFull, f);
    }
    /// high precision (6 fractional digits) conversion to String
    inline String toString(float f)
    {
      String str;
      append(f, str);
      return str;
    }



    /// high precision (15 fractional digits) conversion to String
    inline void append(double d, String& target)
    {
      std::back_insert_iterator<std::string> sink(target);
      boost::spirit::karma::generate(sink, Detail::BK_PrecPolicyDoubleFull, d);
    }
    /// high precision (15 fractional digits) conversion to String
    inline String toString(double d)
    {
      String str;
      append(d, str);
      return str;
    }


    /// high precision (15 fractional digits) conversion to String
    inline void append(long double ld, String& target)
    {
      std::back_insert_iterator<std::string> sink(target);
      boost::spirit::karma::generate(sink, Detail::BK_PrecPolicyLongDoubleFull, ld);
    }
    /// high precision (15 fractional digits) conversion to String
    inline String toString(long double ld)
    {
      String str;
      append(ld, str);
      return str;
    }

    
    inline void append(const DataValue& d, bool full_precision, String& target)
    {
      target += d.toString(full_precision);
    }
    inline String toString(const DataValue& d, bool full_precision)
    {
      return d.toString(full_precision);
    }



    inline String toString(const char c)
    {
      return std::string(1, c);
    }

    inline String toString(const std::string& s)
    {
      return s;
    }

    inline String toString(const char* s)
    {
      return std::string(s);
    }

    /// Other toString functions (different number of arguments)
    inline String toString()
    {
      return String();
    }

  }

} // namespace OPENMS

