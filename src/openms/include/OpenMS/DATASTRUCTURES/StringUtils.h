// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#ifndef OPENMS_DATASTRUCTURES_STRINGUTILS_H
#define OPENMS_DATASTRUCTURES_STRINGUTILS_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>

#include <QtCore/QString>
#include <boost/spirit/include/qi.hpp>

#include <string>
#include <vector>


namespace OpenMS
{
  class String;

  namespace StringConversions
  {

    /// toString functions (for floating point types)
    template <typename T>
    inline String floatToString(T f)
    {
      std::stringstream s;
      s.precision(writtenDigits(f));
      s << f;
      return s.str();
    }

    /// toString functions (single argument)
    template <typename T>
    inline String toString(T i)
    {
      std::stringstream s;
      s << i;
      return s.str();
    }

    template <>
    inline String toString<const char>(const char c)
    {
      return std::string(1,c);
    }

    template <>
    inline String toString(const std::string& s)
    {
      return s;
    }

    template <>
    inline String toString(const char* s)
    {
      return std::string(s);
    }

    template <>
    inline String toString(float f)
    {
      return floatToString(f);
    }

    template <>
    inline String toString(double f)
    {
      return floatToString(f);
    }

    template <>
    inline String toString(long double f)
    {
      return floatToString(f);
    }

    template <>
    inline String toString(const DataValue& d)
    {
      return d.toString();
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
      while (count < length && *(s + count) != 0)
      {
        res += *(s + count);
        ++count;
      }
      return res;
    }
  }

  class OPENMS_DLLAPI StringUtils
  {

public:

    //
    /// Functions
    //
    static String numberLength(double d, UInt n)
    {
      std::stringstream s;
      //reserve one space for the minus sign
      Int sign = 0;
      if (d < 0)
        sign = 1;
      d = fabs(d);

      if (d < pow(10.0, Int(n - sign - 2)))
      {
        s.precision(writtenDigits(d));
        if (sign == 1)
          s << "-";
        s << d;
      }
      else
      {
        UInt exp = 0;
        while (d > pow(10.0, Int(n - sign - 4)))
        {
          d /= 10;
          ++exp;
        }
        d = Int(d) / 10.0;
        exp += 1;
        if (sign == 1)
          s << "-";
        s << d << "e";
        if (exp < 10)
          s << "0";
        s << exp;
      }
      return s.str().substr(0, n);
    }

    static String number(double d, UInt n)
    {
      return QString::number(d, 'f', n);
    }

    static String& fillLeft(String & this_s, char c, UInt size)
    {
      if (this_s.size() < size)
      {
        this_s.std::string::operator=(String(size - this_s.size(), c) + this_s);
      }
      return this_s;
    }

    static String& fillRight(String & this_s, char c, UInt size)
    {
      if (this_s.size() < size)
      {
        this_s.std::string::operator=(this_s + String(size - this_s.size(), c));
      }
      return this_s;
    }

    /*
    /// Constructor from char* (only @p length characters)
    static String toString(const char * s, size_t length)
    {
      std::string res;
      size_t count = 0;
      while (count < length && *(s + count) != 0)
      {
        res += *(s + count);
        ++count;
      }
      return String(res);
    }

    /// Constructor from char (repeats the char @p len times)
    static String toString(size_t len, char c)
    {
      return String(std::string(len, c));
    }
    */

    static bool hasPrefix(const String & this_s, const String & string)
    {
      if (string.size() > this_s.size())
      {
        return false;
      }
      if (string.empty())
      {
        return true;
      }
      return this_s.compare(0, string.size(), string) == 0;
    }

    static bool hasSuffix(const String & this_s, const String& string)
    {
      if (string.size() > this_s.size())
      {
        return false;
      }
      if (string.empty())
      {
        return true;
      }
      return this_s.compare(this_s.size() - string.size(), string.size(), string) == 0;
    }

    static bool hasSubstring(const String & this_s, const String& string)
    {
      return this_s.find(string) != std::string::npos;
    }

    static bool has(const String & this_s, Byte byte)
    {
      return this_s.find(char(byte)) != std::string::npos;
    }

    static String prefix(const String & this_s, size_t length)
    {
      if (length > this_s.size())
      {
        throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, length, this_s.size());
      }
      return this_s.substr(0, length);
    }

    static String suffix(const String & this_s, size_t length)
    {
      if (length > this_s.size())
      {
        throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, length, this_s.size());
      }
      return this_s.substr(this_s.size() - length, length);
    }

    static String prefix(const String & this_s, Int length)
    {
      if (length < 0)
      {
        throw Exception::IndexUnderflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, length, 0);
      }
      if (length > Int(this_s.size()))
      {
        throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, length, this_s.size());
      }
      return this_s.substr(0, length);
    }

    static String suffix(const String & this_s, Int length)
    {
      if (length < 0)
      {
        throw Exception::IndexUnderflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, length, 0);
      }
      if (length > Int(this_s.size()))
      {
        throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, length, this_s.size());
      }
      return this_s.substr(this_s.size() - length, length);
    }

    static String prefix(const String & this_s, char delim)
    {
      Size pos = this_s.find(delim);
      if (pos == std::string::npos) //char not found
      {
        throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         String(delim));
      }
      return this_s.substr(0, pos);
    }

    static String suffix(const String & this_s, char delim)
    {
      Size pos = this_s.rfind(delim);
      if (pos == std::string::npos) //char not found
      {
        throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         String(delim));
      }
      return this_s.substr(++pos);
    }

    static String substr(const String & this_s, size_t pos, size_t n)
    {
      Size begin = std::min(pos, this_s.size());
      return static_cast<String>(this_s.std::string::substr(begin, n));
    }

    static String chop(const String & this_s, Size n)
    {
      Size end = 0;
      if (n < this_s.size())
      {
        end = this_s.size() - n;
      }
      return String(this_s.begin(), this_s.begin() + end);
    }

    static String& trim(String & this_s)
    {
      //search for the begin of truncated string
      std::string::iterator begin = this_s.begin();
      while (begin != this_s.end() && (*begin == ' ' || *begin == '\t' || *begin == '\n'  || *begin == '\r'))
      {
        ++begin;
      }

      //all characters are whitespaces
      if (begin == this_s.end())
      {
        this_s.clear();
        return this_s;
      }

      //search for the end of truncated string
      std::string::iterator end = this_s.end();
      end--;
      while (end != begin && (*end == ' ' || *end == '\n' || *end == '\t' || *end == '\r'))
      {
        --end;
      }
      ++end;

      //no characters are whitespaces
      if (begin == this_s.begin() && end == this_s.end())
      {
        return this_s;
      }

      // TODO:
      // string::operator=(std::string(begin, end));
      this_s.std::string::operator=(std::string(begin, end));

      return this_s;
    }

    static String& quote(String & this_s, char q, String::QuotingMethod method)
    {
      if (method == String::ESCAPE)
      {
        this_s.substitute(String("\\"), String("\\\\"));
        this_s.substitute(String(q), "\\" + String(q));
      }
      else if (method == String::DOUBLE)
        this_s.substitute(String(q), String(q) + String(q));
      this_s.std::string::operator=(q + this_s + q);
      return this_s;
    }

    static String& unquote(String & this_s, char q, String::QuotingMethod method)
    {
      // check if input string matches output format of the "quote" method:
      if ((this_s.size() < 2) || (this_s[0] != q) || (this_s[this_s.size() - 1] != q))
      {
        throw Exception::ConversionError(
                __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                "'" + this_s + "' does not have the expected format of a quoted string");
      }
      this_s.std::string::operator=(this_s.substr(1, this_s.size() - 2)); // remove quotation marks
      if (method == String::ESCAPE)
      {
        this_s.substitute("\\" + String(q), String(q));
        this_s.substitute(String("\\\\"), String("\\"));
      }
      else if (method == String::DOUBLE)
        this_s.substitute(String(q) + String(q), String(q));
      return this_s;
    }

    static String& simplify(String & this_s)
    {
      String simple;

      bool last_was_whitespace = false;
      for (std::string::iterator it = this_s.begin(); it != this_s.end(); ++it)
      {
        if (*it == ' ' || *it == '\n' || *it == '\t' || *it == '\r')
        {
          if (!last_was_whitespace)
          {
            simple += ' ';
          }
          last_was_whitespace = true;
        }
        else
        {
          simple += *it;
          last_was_whitespace = false;
        }
      }

      this_s.swap(simple);
      return this_s;
    }

    static String random(UInt length)
    {
      srand(time(nullptr));
      String tmp(length, '.');
      size_t random;
      for (Size i = 0; i < length; ++i)
      {
        random = static_cast<size_t>(floor((static_cast<double>(rand()) / (double(RAND_MAX) + 1)) * 62.0));
        if (random < 10)
        {
          tmp[i] = static_cast<char>(random + 48);
        }
        else if (random < 36)
        {
          tmp[i] = static_cast<char>(random + 55);
        }
        else
        {
          tmp[i] = static_cast<char>(random + 61);
        }
      }
      return tmp;
    }

    static String& reverse(String & this_s)
    {
      String tmp = this_s;
      for (Size i = 0; i != this_s.size(); ++i)
      {
        this_s[i] = tmp[this_s.size() - 1 - i];
      }
      return this_s;
    }

    static bool split(const String & this_s, const char splitter, std::vector<String>& substrings,
                       bool quote_protect)
    {
      substrings.clear();
      if (this_s.empty())
        return false;

      Size nsplits = count(this_s.begin(), this_s.end(), splitter);

      if (!quote_protect && (nsplits == 0))
      {
        substrings.push_back(this_s);
        return false;
      }

      // splitter(s) found
      substrings.reserve(nsplits + 1);

      // why is "this_s." needed here?
      std::string::const_iterator begin = this_s.begin();
      std::string::const_iterator end = this_s.begin();

      if (quote_protect)
      {
        Int quote_count(0);
        for (; end != this_s.end(); ++end)
        {
          if (*end == '"')
          {
            ++quote_count;
          }
          if ((quote_count % 2 == 0) && (*end == splitter))
          {
            String block = String(begin, end);
            block.trim();
            if ((block.size() >= 2) && ((block.prefix(1) == String("\"")) ^
                                        (block.suffix(1) == String("\""))))
            { // block has start or end quote, but not both
              // (one quote is somewhere in the middle)
              throw Exception::ConversionError(
                      __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                      String("Could not dequote string '") + block +
                      "' due to wrongly placed '\"'.");
            }
            else if ((block.size() >= 2) && (block.prefix(1) == String("\"")) &&
                     (block.suffix(1) == String("\"")))
            { // block has start and end quotes --> remove them
              block = block.substr(1, block.size() - 2);
            }
            substrings.push_back(block);
            begin = end + 1;
          }
        }
        // no valid splitter found - return empty list
        if (substrings.empty())
        {
          substrings.push_back(this_s);
          return false;
        }

        String block = String(begin, end);
        block.trim();
        if ((block.size() >= 2) && ((block.prefix(1) == String("\"")) ^
                                    (block.suffix(1) == String("\""))))
        { // block has start or end quote but not both
          // (one quote is somewhere in the middle)
          throw Exception::ConversionError(
                  __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                  String("Could not dequote string '") + block +
                  "' due to wrongly placed '\"'.");
        }
        else if ((block.size() >= 2) && (block.prefix(1) == String("\"")) &&
                 (block.suffix(1) == String("\"")))
        { // block has start and end quotes --> remove them
          block = block.substr(1, block.size() - 2);
        }
        substrings.push_back(block);
      }
      else // do not honor quotes
      {
        for (; end != this_s.end(); ++end)
        {
          if (*end == splitter)
          {
            substrings.push_back(String(begin, end));
            begin = end + 1;
          }
        }
        substrings.push_back(String(begin, end));
      }

      // at this point we are sure that there are at least two components
      return true; 
    }

    static bool split(const String & this_s, const String& splitter, std::vector<String>& substrings)
    {
      substrings.clear();
      if (this_s.empty())
        return false;

      if (splitter.empty()) // split after every character:
      {
        substrings.resize(this_s.size());
        for (Size i = 0; i < this_s.size(); ++i)
          substrings[i] = this_s[i];
        return true;
      }

      Size len = splitter.size(), start = 0, pos = this_s.find(splitter);
      if (len == 0)
        len = 1;
      while (pos != std::string::npos)
      {
        substrings.push_back(this_s.substr(start, pos - start));
        start = pos + len;
        pos = this_s.find(splitter, start);
      }
      substrings.push_back(this_s.substr(start, this_s.size() - start));
      return substrings.size() > 1;
    }

    static bool split_quoted(const String & this_s, const String& splitter, std::vector<String>& substrings,
                              char q, String::QuotingMethod method)
    {
      substrings.clear();
      if (this_s.empty() || splitter.empty())
        return false;

      bool in_quote = false;
      char targets[2] = {q, splitter[0]}; // targets for "find_first_of"
      std::string rest = splitter.substr(1, splitter.size() - 1);
      Size start = 0;
      for (Size i = 0; i < this_s.size(); ++i)
      {
        if (in_quote) // skip to closing quotation mark
        {
          bool embedded = false;
          if (method == String::ESCAPE)
          {
            for (; i < this_s.size(); ++i)
            {
              if (this_s[i] == '\\')
                embedded = !embedded;
              else if ((this_s[i] == q) && !embedded)
                break;
              else
                embedded = false;
            }
          }
          else // method: NONE or DOUBLE
          {
            for (; i < this_s.size(); ++i)
            {
              if (this_s[i] == q)
              {
                if (method == String::NONE)
                  break; // found
                // next character is also closing quotation mark:
                if ((i < this_s.size() - 1) && (this_s[i + 1] == q))
                  embedded = !embedded;
                // even number of subsequent quotes (doubled) => found
                else if (!embedded)
                  break;
                // odd number of subsequent quotes => belongs to a pair
                else
                  embedded = false;
              }
            }
          }
          in_quote = false; // end of quote reached
        }
        else
        {
          i = this_s.find_first_of(targets, i, 2);
          if (i == std::string::npos)
            break; // nothing found
          if (this_s[i] == q)
            in_quote = true;
          else if (this_s.compare(i + 1, rest.size(), rest) == 0) // splitter found
          {
            substrings.push_back(this_s.substr(start, i - start));
            start = i + splitter.size();
            i = start - 1; // increased by loop
          }
        }
      }
      if (in_quote) // reached end without finding closing quotation mark
      {
        throw Exception::ConversionError(
                __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                "unbalanced quotation marks in string '" + this_s + "'");
      }
      substrings.push_back(this_s.substr(start, this_s.size() - start));
      return substrings.size() > 1;
    }

    static QString toQString(const String & this_s) 
    {
      return QString(this_s.c_str());
    }

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

    static double toDouble(const String& this_s)
    {
      double ret;
      // boost::spirit::qi was found to be vastly superior to boost::lexical_cast or stringstream extraction (especially for VisualStudio),
      // so don't change this unless you have benchmarks for all platforms!
      String::ConstIterator it = this_s.begin();
      if (!boost::spirit::qi::phrase_parse(it, this_s.end(), parse_double_, boost::spirit::ascii::space, ret))
      {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Could not convert string '") + this_s + "' to a double value");
      }
      // was the string parsed (white spaces are skipped automatically!) completely? If not, we have a problem because a previous split might have used the wrong split char
      if (it != this_s.end())
      {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Prefix of string '") + this_s + "' successfully converted to a double value. Additional characters found at position " + (int)(distance(this_s.begin(), it) + 1));
      }
      return ret;
    }


    static String& toUpper(String & this_s)
    {
      std::transform(this_s.begin(), this_s.end(), this_s.begin(), (int (*)(int))toupper);
      return this_s;
    }

    static String& firstToUpper(String & this_s)
    {
      if (this_s.size() != 0)
      {
        this_s[0] = toupper(this_s[0]);
      }
      return this_s;
    }

    static String& toLower(String & this_s)
    {
      std::transform(this_s.begin(), this_s.end(), this_s.begin(), (int (*)(int))tolower);
      return this_s;
    }

    static String& substitute(String & this_s, char from, char to)
    {
      std::replace(this_s.begin(), this_s.end(), from, to);
      return this_s;
    }

  static String& substitute(String & this_s, const String& from, const String& to)
  {
    if (!from.empty())
    {
      std::vector<String> parts;
      this_s.split(from, parts);
      this_s.concatenate(parts.begin(), parts.end(), to);
    }
    return this_s;
  }

  static String& remove(String & this_s, char what)
  {
    this_s.erase(std::remove(this_s.begin(), this_s.end(), what), this_s.end());
    return this_s;
  }

  static String& ensureLastChar(String & this_s, char end)
  {
    if (!this_s.hasSuffix(end))
      this_s.append(1, end);
    return this_s;
  }

  static String& removeWhitespaces(String & this_s)
  {
    bool contains_ws = false;
    for (std::string::const_iterator it = this_s.begin(); it != this_s.end(); ++it)
    {
      char c = *it;
      if (c == ' ' || c == '\t' || c == '\n' || c == '\r')
      {
        contains_ws = true;
        break;
      }
    }

    if (contains_ws)
    {
      std::string tmp;
      tmp.reserve(this_s.size());
      for (std::string::const_iterator it = this_s.begin(); it != this_s.end(); ++it)
      {
        char c = *it;
        if (c != ' ' && c != '\t' && c != '\n' && c != '\r')
          tmp += c;
      }
      this_s.swap(tmp);
    }

    return this_s;
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

} // namespace OPENMS

#endif // OPENMS_DATASTRUCTURES_STRINGUTILS_H
