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

#include <algorithm>
#include <cmath>
#include <sstream>
#include <string_view>
#include <string>
#include <vector>

namespace OpenMS
{
  class String;

  namespace StringUtils
  {

    //
    /// Functions
    //
    static inline String numberLength(double d, UInt n)
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

    static inline String& fillLeft(String & this_s, char c, UInt size)
    {
      if (this_s.size() < size)
      {
        this_s.std::string::operator=(String(size - this_s.size(), c) + this_s);
      }
      return this_s;
    }

    static inline String& fillRight(String & this_s, char c, UInt size)
    {
      if (this_s.size() < size)
      {
        this_s.std::string::operator=(this_s + String(size - this_s.size(), c));
      }
      return this_s;
    }

    static inline bool hasPrefix(const String & this_s, const String & string)
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

    static inline bool hasSuffix(const String & this_s, const String& string)
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

    static inline bool hasSubstring(const String & this_s, const String& string)
    {
      return this_s.find(string) != std::string::npos;
    }

    static inline bool has(const String & this_s, Byte byte)
    {
      return this_s.find(char(byte)) != std::string::npos;
    }

    static inline String prefix(const String & this_s, size_t length)
    {
      if (length > this_s.size())
      {
        throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, length, this_s.size());
      }
      return this_s.substr(0, length);
    }

    static inline String suffix(const String & this_s, size_t length)
    {
      if (length > this_s.size())
      {
        throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, length, this_s.size());
      }
      return this_s.substr(this_s.size() - length, length);
    }

    static inline String prefix(const String & this_s, Int length)
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

    static inline String suffix(const String & this_s, Int length)
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

    static inline String prefix(const String & this_s, char delim)
    {
      Size pos = this_s.find(delim);
      if (pos == std::string::npos) //char not found
      {
        throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         String(delim));
      }
      return this_s.substr(0, pos);
    }

    static inline String suffix(const String & this_s, char delim)
    {
      Size pos = this_s.rfind(delim);
      if (pos == std::string::npos) //char not found
      {
        throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         String(delim));
      }
      return this_s.substr(++pos);
    }

    static inline String substr(const String & this_s, size_t pos, size_t n)
    {
      Size begin = std::min(pos, this_s.size());
      return static_cast<String>(this_s.std::string::substr(begin, n));
    }

    static inline String chop(const String & this_s, Size n)
    {
      Size end = 0;
      if (n < this_s.size())
      {
        end = this_s.size() - n;
      }
      return String(this_s.begin(), this_s.begin() + end);
    }

    static inline String& trim(String & this_s)
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

    static inline bool isQuoted(const String & this_s, char q)
    {
      return (this_s.size() < 2) || (this_s[0] != q) || (this_s[this_s.size() - 1] != q);
    }

    static inline String& quote(String & this_s, char q, String::QuotingMethod method)
    {
      if (method == String::ESCAPE)
      {
        this_s.substitute(String(R"(\)"), String(R"(\\)"));
        this_s.substitute(String(q), R"(\)" + String(q));
      }
      else if (method == String::DOUBLE)
        this_s.substitute(String(q), String(q) + String(q));
      this_s.std::string::operator=(q + this_s + q);
      return this_s;
    }

    static inline String& unquote(String & this_s, char q, String::QuotingMethod method)
    {
      // check if input string matches output format of the "quote" method:
      if (isQuoted(this_s, q))
      {
        throw Exception::ConversionError(
                __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                "'" + this_s + "' does not have the expected format of a quoted string");
      }
      this_s.std::string::operator=(this_s.substr(1, this_s.size() - 2)); // remove quotation marks
      if (method == String::ESCAPE)
      {
        this_s.substitute(R"(\)" + String(q), String(q));
        this_s.substitute(String(R"(\\)"), String(R"(\)"));
      }
      else if (method == String::DOUBLE)
        this_s.substitute(String(q) + String(q), String(q));
      return this_s;
    }

    static inline String& simplify(String & this_s)
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

    static inline String random(UInt length)
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

    static inline String& reverse(String & this_s)
    {
      String tmp = this_s;
      for (Size i = 0; i != this_s.size(); ++i)
      {
        this_s[i] = tmp[this_s.size() - 1 - i];
      }
      return this_s;
    }

    static inline bool split(const String & this_s, const char splitter, std::vector<String>& substrings,
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

    static inline bool split(const String & this_s, const String& splitter, std::vector<String>& substrings)
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

    static inline bool split_quoted(const String & this_s, const String& splitter, std::vector<String>& substrings,
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

    static inline String& toUpper(String & this_s)
    {
      std::transform(this_s.begin(), this_s.end(), this_s.begin(), (int (*)(int))toupper);
      return this_s;
    }

    static inline String& firstToUpper(String & this_s)
    {
      if (!this_s.empty())
      {
        this_s[0] = toupper(this_s[0]);
      }
      return this_s;
    }

    static inline String& toLower(String & this_s)
    {
      std::transform(this_s.begin(), this_s.end(), this_s.begin(), (int (*)(int))tolower);
      return this_s;
    }

    static inline String& substitute(String & this_s, char from, char to)
    {
      std::replace(this_s.begin(), this_s.end(), from, to);
      return this_s;
    }

    static inline String& substitute(String & this_s, const String& from, const String& to)
    {
      if (!from.empty())
      {
        std::vector<String> parts;
        this_s.split(from, parts);
        this_s.concatenate(parts.begin(), parts.end(), to);
      }
      return this_s;
    }
  
    static inline String& remove(String & this_s, char what)
    {
      this_s.erase(std::remove(this_s.begin(), this_s.end(), what), this_s.end());
      return this_s;
    }
  
    static inline String& ensureLastChar(String & this_s, char end)
    {
      if (!this_s.hasSuffix(end))
        this_s.append(1, end);
      return this_s;
    }

    /**
     @brief Get the first non-whitespace character (anything but \\n, \\t, \\r, ' ') in the string pointed to by @p p (where @p p_end is past the end of the string).

     If only whitespaces are contained, then @p p_end is returned.
    */
    OPENMS_DLLAPI const char* skipWhitespace(const char* p, const char* p_end);

    /**
     @brief Get the number of whitespace characters (\\n, \\t, \\r, ' ') in the prefix of @p data
    */
    inline int skipWhitespace(const std::string_view& data)
    {
      auto pos = skipWhitespace(data.data(), data.data() + data.size());
      return pos - data.data();
    }

    /**
     @brief Get the first whitespace character (\\n, \\t, \\r, ' ') in the string pointed to by @p p (where @p p_end is past the end of the string).

     If only non-whitespaces are contained, then @p p_end is returned.
    */
    OPENMS_DLLAPI const char* skipNonWhitespace(const char* p, const char* p_end);

    /**
     @brief return the number of non-whitespace characters (anything but \\n, \\t, \\r, ' ') in the prefix of @p data
    */
    inline int skipNonWhitespace(const std::string_view& data)
    {
      auto pos = skipNonWhitespace(data.data(), data.data() + data.size());
      return pos - data.data();
    }
    
    static inline String& removeWhitespaces(String& this_s)
    {
      auto start = skipNonWhitespace(std::string_view(this_s.data()));
      std::string::const_iterator it = this_s.begin() + start;
      std::string::iterator dest = this_s.begin() + start;
      std::string::const_iterator it_end = this_s.end();
      bool has_spaces(false);
      while (it != it_end)
      {
        const char c = *it;
        if (c == ' ' || c == '\t' || c == '\n' || c == '\r')
        {
          ++it;
          has_spaces = true;
          continue; // no need to copy a whitespace
        }
        // copy to the left, if we had a whitespace before
        if (has_spaces) *dest = *it;
        // advance both
        ++dest;
        ++it;
      }
  
      // shorten result
      if (has_spaces) this_s.resize(dest - this_s.begin());
  
      return this_s;
    }
  
  }

} // namespace OPENMS

