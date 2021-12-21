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
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/CONCEPT/PrecisionWrapper.h>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/type_traits.hpp>

#include <string>
#include <vector>


namespace OpenMS
{
  class String;

  namespace StringConversions
  {

    // Karma full precision float policy
    template <typename T> 
    class BK_PrecPolicy : public boost::spirit::karma::real_policies<T>
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
            if (BK_PrecPolicy::floatfield(n))
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
            if (boost::spirit::traits::test_zero(n))
                return base_policy_type::fmtflags::fixed;

            T abs_n = boost::spirit::traits::get_absolute_value(n);
            // this is due to a bug in downstream thirdparty tools that only can read
            // up to 19 digits. https://github.com/OpenMS/OpenMS/issues/4627
            return (abs_n >= 1e4 || abs_n < 1e-2) 
                ? base_policy_type::fmtflags::scientific : base_policy_type::fmtflags::fixed;
        }
    };
    typedef boost::spirit::karma::real_generator<float, BK_PrecPolicy<float> > BK_PrecPolicyFloat_type;
    const BK_PrecPolicyFloat_type BK_PrecPolicyFloat;
    typedef boost::spirit::karma::real_generator<double, BK_PrecPolicy<double> > BK_PrecPolicyDouble_type;
    const BK_PrecPolicyDouble_type BK_PrecPolicyDouble;
    typedef boost::spirit::karma::real_generator<long double, BK_PrecPolicy<long double> > BK_PrecPolicyLongDouble_type;
    const BK_PrecPolicyLongDouble_type BK_PrecPolicyLongDouble;
    
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
      boost::spirit::karma::generate(sink, f);
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
      boost::spirit::karma::generate(sink, d);
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
      boost::spirit::karma::generate(sink, ld);
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
      boost::spirit::karma::generate(sink, BK_PrecPolicyFloat, f);
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
      boost::spirit::karma::generate(sink, BK_PrecPolicyDouble, d);
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
      boost::spirit::karma::generate(sink, BK_PrecPolicyLongDouble, ld);
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

    inline String toString(const char* s, size_t length)
    {
      String res;
      size_t count = 0;
      while (count < length)
      {
        res += *(s + count);
        ++count;
      }
      return res;
    }
  }

} // namespace OPENMS

