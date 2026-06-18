// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg, Chris Bielow $
// $Authors: Marc Sturm, Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

// NOTE: Do NOT include String.h here — String.h includes us. Use std::string directly.

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/OpenMSConfig.h>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/home/karma/numeric/detail/real_utils.hpp>
#include <boost/type_traits.hpp>
#include <boost/functional/hash.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <ranges>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace OpenMS
{
  class DataValue; // forward — DataValue.h includes String.h, so we cannot include it here
  class ParamValue; // forward — same reason

  // -------------------------------------------------------------------------
  // QuotingMethod — formerly String::QuotingMethod
  // -------------------------------------------------------------------------
  /// How to handle embedded quotes when quoting strings
  enum class QuotingMethod { NONE, ESCAPE, DOUBLE };


  // =========================================================================
  // Internal Boost.Karma generator policies (precision, NaN handling)
  // =========================================================================
  namespace StringConversions
  {
    namespace Detail
    {
      /// Karma full-precision float/double/long double policy with NaN→"NaN"
      template<typename T>
      class BK_PrecPolicyFull : public boost::spirit::karma::real_policies<T>
      {
        typedef boost::spirit::karma::real_policies<T> base_policy_type;
      public:
        static unsigned precision(T /*n*/) { return writtenDigits<T>(); }

        static unsigned floatfield(T n)
        {
          if (boost::spirit::traits::test_zero(n)) return base_policy_type::fmtflags::fixed;
          T abs_n = boost::spirit::traits::get_absolute_value(n);
          // downstream tools (Java-based) read only up to 19 digits — keep scientific for large/small values
          return (abs_n >= 1e4 || abs_n < 1e-2) ? base_policy_type::fmtflags::scientific
                                                 : base_policy_type::fmtflags::fixed;
        }

        template<typename CharEncoding, typename Tag, typename OutputIterator>
        static bool nan(OutputIterator& sink, T n, bool force_sign)
        {
          return boost::spirit::karma::sign_inserter::call(
                   sink, false, boost::spirit::traits::test_negative(n), force_sign)
                 && boost::spirit::karma::string_inserter<CharEncoding, Tag>::call(sink, "NaN");
        }
      };

      /// Karma low-precision (≤3 fractional digits) policy with NaN→"NaN"
      template<typename T>
      class BK_PrecPolicyShort : public boost::spirit::karma::real_policies<T>
      {
        typedef boost::spirit::karma::real_policies<T> base_policy_type;
      public:
        template<typename CharEncoding, typename Tag, typename OutputIterator>
        static bool nan(OutputIterator& sink, T n, bool force_sign)
        {
          return boost::spirit::karma::sign_inserter::call(
                   sink, false, boost::spirit::traits::test_negative(n), force_sign)
                 && boost::spirit::karma::string_inserter<CharEncoding, Tag>::call(sink, "NaN");
        }
      };

      using FloatFullGen  = boost::spirit::karma::real_generator<float,       BK_PrecPolicyFull<float>>;
      using DoubleFullGen = boost::spirit::karma::real_generator<double,      BK_PrecPolicyFull<double>>;
      using LDFullGen     = boost::spirit::karma::real_generator<long double, BK_PrecPolicyFull<long double>>;

      using FloatShortGen  = boost::spirit::karma::real_generator<float,       BK_PrecPolicyShort<float>>;
      using DoubleShortGen = boost::spirit::karma::real_generator<double,      BK_PrecPolicyShort<double>>;
      using LDShortGen     = boost::spirit::karma::real_generator<long double, BK_PrecPolicyShort<long double>>;

      inline const FloatFullGen   floatFull{};
      inline const DoubleFullGen  doubleFull{};
      inline const LDFullGen      ldFull{};
      inline const FloatShortGen  floatShort{};
      inline const DoubleShortGen doubleShort{};
      inline const LDShortGen     ldShort{};
    } // namespace Detail
  } // namespace StringConversions


  // =========================================================================
  // Internal Boost.Spirit.Qi parser helper (NaN-safe)
  // =========================================================================
  class OPENMS_DLLAPI StringUtilsHelper
  {
  public:
    /// Parse int32 from string (leading/trailing whitespace allowed)
    static Int32 toInt32(const std::string& s);
    /// Parse int64 from string (leading/trailing whitespace allowed)
    static Int64 toInt64(const std::string& s);
    /// Parse float from string (leading/trailing whitespace allowed)
    static float toFloat(const std::string& s);
    /// Parse double from string (leading/trailing whitespace allowed)
    static double toDouble(const std::string& s);

    /// Parse a double starting at @p begin; advances @p begin on success.
    /// Whitespace is NOT consumed. Returns false if no double found.
    template <typename IteratorT>
    static bool extractDouble(IteratorT& begin, const IteratorT& end, double& target)
    {
      return boost::spirit::qi::parse(begin, end, parse_double_, target);
    }

    /// Parse an int starting at @p begin; advances @p begin on success.
    /// Whitespace is NOT consumed. Returns false if no int found.
    template <typename IteratorT>
    static bool extractInt(IteratorT& begin, const IteratorT& end, int& target)
    {
      return boost::spirit::qi::parse(begin, end, parse_int_, target);
    }

  private:
    template <typename T>
    struct real_policies_NANfixed_ : boost::spirit::qi::real_policies<T>
    {
      template <typename Iterator, typename Attribute>
      static bool parse_nan(Iterator& first, Iterator const& last, Attribute& attr_)
      {
        if (first == last) return false;
        if (*first != 'n' && *first != 'N') return false;
        if (boost::spirit::qi::detail::string_parse("nan", "NAN", first, last, boost::spirit::qi::unused))
        {
          if (first != last && *first == '(')
          {
            Iterator i = first;
            while (++i != last && *i != ')') ;
            if (i == last) return false;
            first = ++i;
          }
          attr_ = std::numeric_limits<T>::quiet_NaN();
          return true;
        }
        return false;
      }
    };

    static boost::spirit::qi::real_parser<double, real_policies_NANfixed_<double>> parse_double_;
    static boost::spirit::qi::real_parser<float,  real_policies_NANfixed_<float>>  parse_float_;
    static boost::spirit::qi::int_parser<>                                          parse_int_;
  };


  // =========================================================================
  // StringUtils — unified string utility namespace
  // =========================================================================
  namespace StringUtils
  {

    // -----------------------------------------------------------------------
    // Append numeric value to a string (no copy; fastest path)
    // -----------------------------------------------------------------------

    /// Append @p i to @p target (Boost.Karma, fast)
    template <typename T>
    inline void appendToStr(const T& i, std::string& target)
    {
      std::back_insert_iterator<std::string> sink(target);
      boost::spirit::karma::generate(sink, i);
    }

    /// Append integer @p i to @p target
    inline void appendToStr(int i, std::string& target)
    { appendToStr<int>(i, target); }
    inline void appendToStr(unsigned int i, std::string& target)
    { appendToStr<unsigned int>(i, target); }
    inline void appendToStr(short int i, std::string& target)
    { appendToStr<short int>(i, target); }
    inline void appendToStr(short unsigned int i, std::string& target)
    { appendToStr<short unsigned int>(i, target); }
    inline void appendToStr(long int i, std::string& target)
    { appendToStr<long int>(i, target); }
    inline void appendToStr(long unsigned int i, std::string& target)
    { appendToStr<long unsigned int>(i, target); }
    inline void appendToStr(long long unsigned int i, std::string& target)
    { appendToStr<long long unsigned int>(i, target); }
    inline void appendToStr(long long signed int i, std::string& target)
    { appendToStr<long long signed int>(i, target); }

    /// Append float (high precision, 6 fractional digits) to @p target
    inline void appendToStr(float f, std::string& target)
    {
      std::back_insert_iterator<std::string> sink(target);
      boost::spirit::karma::generate(sink, StringConversions::Detail::floatFull, f);
    }
    /// Append float (low precision, 3 fractional digits) to @p target
    inline void appendToStrLowP(float f, std::string& target)
    {
      std::back_insert_iterator<std::string> sink(target);
      boost::spirit::karma::generate(sink, StringConversions::Detail::floatShort, f);
    }

    /// Append double (high precision, 15 fractional digits) to @p target
    inline void appendToStr(double d, std::string& target)
    {
      std::back_insert_iterator<std::string> sink(target);
      boost::spirit::karma::generate(sink, StringConversions::Detail::doubleFull, d);
    }
    /// Append double (low precision, 3 fractional digits) to @p target
    inline void appendToStrLowP(double d, std::string& target)
    {
      std::back_insert_iterator<std::string> sink(target);
      boost::spirit::karma::generate(sink, StringConversions::Detail::doubleShort, d);
    }

    /// Append long double (high precision) to @p target
    inline void appendToStr(long double ld, std::string& target)
    {
      std::back_insert_iterator<std::string> sink(target);
      boost::spirit::karma::generate(sink, StringConversions::Detail::ldFull, ld);
    }
    /// Append long double (low precision) to @p target
    inline void appendToStrLowP(long double ld, std::string& target)
    {
      std::back_insert_iterator<std::string> sink(target);
      boost::spirit::karma::generate(sink, StringConversions::Detail::ldShort, ld);
    }

    /// Append DataValue string representation to @p target (non-inline; defined in StringUtils.cpp)
    OPENMS_DLLAPI void appendToStr(const DataValue& d, bool full_precision, std::string& target);


    // -----------------------------------------------------------------------
    // toStr() — convert a value to a new std::string
    // -----------------------------------------------------------------------

    inline std::string toStr(int i)
    { std::string s; appendToStr(i, s); return s; }
    inline std::string toStr(unsigned int i)
    { std::string s; appendToStr(i, s); return s; }
    inline std::string toStr(short int i)
    { std::string s; appendToStr(i, s); return s; }
    inline std::string toStr(short unsigned int i)
    { std::string s; appendToStr(i, s); return s; }
    inline std::string toStr(long int i)
    { std::string s; appendToStr(i, s); return s; }
    inline std::string toStr(long unsigned int i)
    { std::string s; appendToStr(i, s); return s; }
    inline std::string toStr(long long unsigned int i)
    { std::string s; appendToStr(i, s); return s; }
    inline std::string toStr(long long signed int i)
    { std::string s; appendToStr(i, s); return s; }

    /// Float to string; @p full_precision selects 6-digit (true) or 3-digit (false) output
    inline std::string toStr(float f, bool full_precision = true)
    {
      std::string s;
      full_precision ? appendToStr(f, s) : appendToStrLowP(f, s);
      return s;
    }
    /// Double to string; @p full_precision selects 15-digit (true) or 3-digit (false) output
    inline std::string toStr(double d, bool full_precision = true)
    {
      std::string s;
      full_precision ? appendToStr(d, s) : appendToStrLowP(d, s);
      return s;
    }
    /// Long double to string; @p full_precision selects full-precision or 3-digit output
    inline std::string toStr(long double ld, bool full_precision = true)
    {
      std::string s;
      full_precision ? appendToStr(ld, s) : appendToStrLowP(ld, s);
      return s;
    }

    inline std::string toStr(char c)
    { return std::string(1, c); }

    inline std::string toStr(const std::string& s)
    { return s; }

    inline std::string toStr(const char* s)
    { return std::string(s); }

    /**
      @brief DataValue to string — the LENIENT stringification (non-inline; defined in StringUtils.cpp).

      "Class B" gotcha of the OpenMS::String -> std::string migration (PR 9450):
      this is the drop-in replacement for the removed lenient @c String(const DataValue&)
      constructor. It calls @c DataValue::toString() and therefore NEVER throws: numbers
      become their decimal text, lists are joined, EMPTY becomes "", and a genuine STRING
      is returned verbatim (a true no-op).

      Prefer this over the bare cast whenever you stringify a DataValue or a getMetaValue()
      result, because @c DataValue::operator std::string() is STRICT and throws
      Exception::ConversionError for non-string values:
      @code
        // WRONG: throws at runtime if "foo" is ever Int/Double/List/Empty
        std::string s = mi.getMetaValue("foo");
        std::string s = (std::string)dv;
        // RIGHT: lenient, matches the pre-migration String(DataValue) behaviour
        std::string s = StringUtils::toStr(mi.getMetaValue("foo"));
        std::string s = StringUtils::toStr(dv);
      @endcode
    */
    OPENMS_DLLAPI std::string toStr(const DataValue& d, bool full_precision = true);

    /// ParamValue to string — lenient, like toStr(const DataValue&) (non-inline; defined in StringUtils.cpp)
    OPENMS_DLLAPI std::string toStr(const ParamValue& p, bool full_precision = true);


    // -----------------------------------------------------------------------
    // Parse string → numeric
    // -----------------------------------------------------------------------

    inline Int32  toInt32(const std::string& s) { return StringUtilsHelper::toInt32(s); }
    inline Int64  toInt64(const std::string& s) { return StringUtilsHelper::toInt64(s); }
    inline float  toFloat(const std::string& s) { return StringUtilsHelper::toFloat(s); }
    inline double toDouble(const std::string& s) { return StringUtilsHelper::toDouble(s); }

    template <typename IteratorT>
    inline bool extractDouble(IteratorT& begin, const IteratorT& end, double& target)
    { return StringUtilsHelper::extractDouble(begin, end, target); }

    template <typename IteratorT>
    inline bool extractInt(IteratorT& begin, const IteratorT& end, int& target)
    { return StringUtilsHelper::extractInt(begin, end, target); }


    // -----------------------------------------------------------------------
    // Whitespace scanning (SIMD-accelerated; non-inline)
    // -----------------------------------------------------------------------

    /// Returns pointer to first non-whitespace character in [p, p_end), or p_end
    OPENMS_DLLAPI const char* skipWhitespace(const char* p, const char* p_end);

    /// Returns count of leading whitespace characters in @p data
    inline int skipWhitespace(const std::string_view& data)
    {
      auto pos = skipWhitespace(data.data(), data.data() + data.size());
      return static_cast<int>(pos - data.data());
    }

    /// Returns pointer to first whitespace character in [p, p_end), or p_end
    OPENMS_DLLAPI const char* skipNonWhitespace(const char* p, const char* p_end);

    /// Returns count of leading non-whitespace characters in @p data
    inline int skipNonWhitespace(const std::string_view& data)
    {
      auto pos = skipNonWhitespace(data.data(), data.data() + data.size());
      return static_cast<int>(pos - data.data());
    }


    // -----------------------------------------------------------------------
    // Static factory methods
    // -----------------------------------------------------------------------

    /// Returns a string with exactly @p n decimal places for @p d
    [[maybe_unused]] inline std::string number(double d, UInt n)
    {
      char buf[64];
      std::snprintf(buf, sizeof(buf), "%.*f", static_cast<int>(n), d);
      return std::string(buf);
    }

    /// Returns a string for @p d with at most @p n characters total (scientific notation if needed)
    [[maybe_unused]] inline std::string numberLength(double d, UInt n)
    {
      std::stringstream s;
      Int sign = (d < 0) ? 1 : 0;
      d = std::fabs(d);
      if (d < std::pow(10.0, Int(n - sign - 2)))
      {
        s.precision(writtenDigits(d));
        if (sign == 1) s << '-';
        s << d;
      }
      else
      {
        UInt exp = 0;
        while (d > std::pow(10.0, Int(n - sign - 4)))
        {
          d /= 10;
          ++exp;
        }
        d = static_cast<int>(d) / 10.0;
        exp += 1;
        if (sign == 1) s << '-';
        s << d << 'e';
        if (exp < 10) s << '0';
        s << exp;
      }
      return s.str().substr(0, n);
    }

    /// Returns a random string of @p length characters from [0-9a-zA-Z]
    inline std::string random(UInt length)
    {
      srand(time(nullptr));
      std::string tmp(length, '.');
      for (Size i = 0; i < length; ++i)
      {
        size_t r = static_cast<size_t>(
          std::floor((static_cast<double>(rand()) / (double(RAND_MAX) + 1)) * 62.0));
        if      (r < 10) tmp[i] = static_cast<char>(r + 48);
        else if (r < 36) tmp[i] = static_cast<char>(r + 55);
        else             tmp[i] = static_cast<char>(r + 61);
      }
      return tmp;
    }


    // -----------------------------------------------------------------------
    // Predicates
    // -----------------------------------------------------------------------

    inline bool hasPrefix(const std::string& s, const std::string& prefix)
    {
      if (prefix.size() > s.size()) return false;
      if (prefix.empty())           return true;
      return s.compare(0, prefix.size(), prefix) == 0;
    }
    inline bool hasPrefix(const std::string& s, char c)
    { return !s.empty() && s.front() == c; }

    inline bool hasSuffix(const std::string& s, const std::string& sfx)
    {
      if (sfx.size() > s.size()) return false;
      if (sfx.empty())           return true;
      return s.compare(s.size() - sfx.size(), sfx.size(), sfx) == 0;
    }
    inline bool hasSuffix(const std::string& s, char c)
    { return !s.empty() && s.back() == c; }

    inline bool hasSubstring(const std::string& s, const std::string& sub)
    {
      return s.find(sub) != std::string::npos;
    }
    inline bool hasSubstring(const std::string& s, char c)
    { return s.find(c) != std::string::npos; }

    inline bool has(const std::string& s, Byte byte)
    {
      return s.find(static_cast<char>(byte)) != std::string::npos;
    }
    inline bool has(const std::string& s, char c)
    { return s.find(c) != std::string::npos; }


    // -----------------------------------------------------------------------
    // Accessors (return new string)
    // -----------------------------------------------------------------------

    inline std::string prefix(const std::string& s, size_t length)
    {
      if (length > s.size())
        throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, length, s.size());
      return s.substr(0, length);
    }

    inline std::string suffix(const std::string& s, size_t length)
    {
      if (length > s.size())
        throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, length, s.size());
      return s.substr(s.size() - length, length);
    }

    inline std::string prefix(const std::string& s, Int length)
    {
      if (length < 0)
        throw Exception::IndexUnderflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, length, 0);
      return prefix(s, static_cast<size_t>(length));
    }

    inline std::string suffix(const std::string& s, Int length)
    {
      if (length < 0)
        throw Exception::IndexUnderflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, length, 0);
      return suffix(s, static_cast<size_t>(length));
    }

    /**
      @brief Returns the part of @p s before the first occurrence of @p delim (excluding it).

      If @p delim is not contained in @p s, the whole string is returned.
    */
    inline std::string prefix(const std::string& s, char delim)
    {
      return s.substr(0, s.find(delim)); // find() returns npos if absent -> substr(0, npos) -> whole string
    }

    /**
      @brief Returns the part of @p s after the last occurrence of @p delim (excluding it).

      If @p delim is not contained in @p s, the whole string is returned.
    */
    inline std::string suffix(const std::string& s, char delim)
    {
      size_t pos = s.rfind(delim);
      return pos == std::string::npos ? s : s.substr(pos + 1);
    }

    /// Wrapper around std::string::substr; clamps @p pos to [0, size]
    inline std::string substr(const std::string& s, size_t pos = 0, size_t n = std::string::npos)
    {
      size_t begin = std::min(pos, s.size());
      return s.substr(begin, n);
    }
    /// string_view overload — implicit conversions from string_view to const string& don't work
    inline std::string substr(std::string_view s, size_t pos = 0, size_t n = std::string::npos)
    {
      size_t begin = std::min(pos, s.size());
      return std::string(s.substr(begin, n));
    }

    /// toStr overload for string_view (returns a copy as std::string)
    inline std::string toStr(std::string_view sv) { return std::string(sv); }

    /// hasPrefix / hasSuffix / hasSubstring for string_view
    inline bool hasPrefix(std::string_view s, const std::string& prefix)
    { return s.size() >= prefix.size() && s.substr(0, prefix.size()) == prefix; }
    inline bool hasSuffix(std::string_view s, const std::string& sfx)
    { return s.size() >= sfx.size() && s.substr(s.size() - sfx.size()) == sfx; }
    inline bool hasSubstring(std::string_view s, const std::string& sub)
    { return s.find(sub) != std::string_view::npos; }

    /// Remove @p n characters from the end; returns empty string if n >= size
    inline std::string chop(const std::string& s, Size n)
    {
      size_t end = (n < s.size()) ? s.size() - n : 0;
      return std::string(s.begin(), s.begin() + end);
    }


    // -----------------------------------------------------------------------
    // In-place mutators (all return std::string& for chaining via free fns)
    // -----------------------------------------------------------------------

    inline std::string& trim(std::string& s)
    {
      auto begin = s.begin();
      while (begin != s.end() && (*begin == ' ' || *begin == '\t' || *begin == '\n' || *begin == '\r'))
        ++begin;
      if (begin == s.end()) { s.clear(); return s; }
      auto end = s.end() - 1;
      while (end != begin && (*end == ' ' || *end == '\n' || *end == '\t' || *end == '\r'))
        --end;
      ++end;
      if (begin != s.begin() || end != s.end())
        s.assign(begin, end);
      return s;
    }

    /// Returns a trimmed copy of @p s (for use in chained/rvalue expressions)
    inline std::string trimmed(std::string s) { return std::move(trim(s)); }

    inline std::string& toUpper(std::string& s)
    {
      std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::toupper(c); });
      return s;
    }

    inline std::string& toLower(std::string& s)
    {
      std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::tolower(c); });
      return s;
    }
    /// Returns an upper-cased copy of @p s (for use in chained/rvalue expressions)
    inline std::string toUppered(std::string s) { return std::move(toUpper(s)); }
    /// Returns a lower-cased copy of @p s (for use in chained/rvalue expressions)
    inline std::string toLowered(std::string s) { return std::move(toLower(s)); }

    inline std::string& firstToUpper(std::string& s)
    {
      if (!s.empty()) s[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(s[0])));
      return s;
    }

    inline std::string& reverse(std::string& s)
    {
      std::reverse(s.begin(), s.end());
      return s;
    }

    inline std::string& simplify(std::string& s)
    {
      std::string result;
      result.reserve(s.size());
      bool last_ws = false;
      for (char c : s)
      {
        if (c == ' ' || c == '\n' || c == '\t' || c == '\r')
        {
          if (!last_ws) result += ' ';
          last_ws = true;
        }
        else
        {
          result += c;
          last_ws = false;
        }
      }
      s.swap(result);
      return s;
    }

    inline std::string& fillLeft(std::string& s, char c, UInt size)
    {
      if (s.size() < size)
        s.insert(s.begin(), size - s.size(), c);
      return s;
    }

    inline std::string& fillRight(std::string& s, char c, UInt size)
    {
      if (s.size() < size)
        s.append(size - s.size(), c);
      return s;
    }

    inline std::string& substitute(std::string& s, char from, char to)
    {
      std::replace(s.begin(), s.end(), from, to);
      return s;
    }

    /// Replace all occurrences of @p from with @p to in @p s
    inline std::string& substitute(std::string& s, const std::string& from, const std::string& to)
    {
      if (from.empty()) return s;
      std::string result;
      result.reserve(s.size());
      size_t start = 0;
      size_t pos;
      while ((pos = s.find(from, start)) != std::string::npos)
      {
        result.append(s, start, pos - start);
        result += to;
        start = pos + from.size();
      }
      result.append(s, start, std::string::npos);
      s.swap(result);
      return s;
    }
    /// substitute on a copy (for chained/rvalue expressions)
    inline std::string substituted(std::string s, char from, char to) { return std::move(substitute(s, from, to)); }
    inline std::string substituted(std::string s, const std::string& from, const std::string& to) { return std::move(substitute(s, from, to)); }

    inline std::string& remove(std::string& s, char what)
    {
      s.erase(std::remove(s.begin(), s.end(), what), s.end());
      return s;
    }

    // rvalue-ref overloads so chained/temporary strings can be padded in place
    inline std::string fillLeft(std::string&& s, char c, UInt size) { fillLeft(s, c, size); return std::move(s); }
    inline std::string fillRight(std::string&& s, char c, UInt size) { fillRight(s, c, size); return std::move(s); }

    inline std::string& ensureLastChar(std::string& s, char end)
    {
      if (s.empty() || s.back() != end) s.push_back(end);
      return s;
    }
    /// rvalue-ref overload (defined after the lvalue version it delegates to)
    inline std::string ensureLastChar(std::string&& s, char end) { ensureLastChar(s, end); return std::move(s); }

    inline std::string& removeWhitespaces(std::string& s)
    {
      // skip unmodified prefix
      int start = skipNonWhitespace(std::string_view(s.data(), s.size()));
      auto it     = s.cbegin() + start;
      auto dest   = s.begin()  + start;
      auto it_end = s.cend();
      bool has_spaces = false;
      while (it != it_end)
      {
        const char c = *it;
        if (c == ' ' || c == '\t' || c == '\n' || c == '\r')
        {
          ++it;
          has_spaces = true;
          continue;
        }
        if (has_spaces) *dest = *it;
        ++dest;
        ++it;
      }
      if (has_spaces) s.resize(static_cast<size_t>(dest - s.begin()));
      return s;
    }

    inline bool isQuoted(const std::string& s, char q)
    {
      return s.size() >= 2 && s.front() == q && s.back() == q;
    }

    inline std::string& quote(std::string& s, char q = '"', QuotingMethod method = QuotingMethod::ESCAPE)
    {
      if (method == QuotingMethod::ESCAPE)
      {
        substitute(s, std::string("\\"), std::string("\\\\"));
        substitute(s, std::string(1, q), std::string("\\") + q);
      }
      else if (method == QuotingMethod::DOUBLE)
      {
        substitute(s, std::string(1, q), std::string(2, q));
      }
      s.insert(s.begin(), q);
      s.push_back(q);
      return s;
    }

    inline std::string& unquote(std::string& s, char q = '"', QuotingMethod method = QuotingMethod::ESCAPE)
    {
      if (!isQuoted(s, q))
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "'" + s + "' does not have the expected format of a quoted string");
      s.erase(s.begin());
      s.pop_back();
      if (method == QuotingMethod::ESCAPE)
      {
        substitute(s, std::string("\\") + q, std::string(1, q));
        substitute(s, std::string("\\\\"), std::string("\\"));
      }
      else if (method == QuotingMethod::DOUBLE)
      {
        substitute(s, std::string(2, q), std::string(1, q));
      }
      return s;
    }


    // -----------------------------------------------------------------------
    // Split / Join
    // -----------------------------------------------------------------------

    /// Split @p s on single character @p splitter into @p substrings.
    /// If @p quote_protect is true, splits only outside double-quoted blocks (quotes stripped).
    /// Returns true if at least one split occurred.
    inline bool split(const std::string& s, char splitter,
                      std::vector<std::string>& substrings,
                      bool quote_protect = false)
    {
      substrings.clear();
      if (s.empty()) return false;

      size_t nsplits = static_cast<size_t>(std::count(s.begin(), s.end(), splitter));

      if (!quote_protect && nsplits == 0)
      {
        substrings.push_back(s);
        return false;
      }

      substrings.reserve(nsplits + 1);

      if (quote_protect)
      {
        int quote_count = 0;
        auto begin = s.cbegin();
        auto end   = s.cbegin();
        for (; end != s.cend(); ++end)
        {
          if (*end == '"') ++quote_count;
          if ((quote_count % 2 == 0) && (*end == splitter))
          {
            std::string block(begin, end);
            trim(block);
            bool has_start = block.size() >= 1 && block.front() == '"';
            bool has_end   = block.size() >= 1 && block.back()  == '"';
            if (block.size() >= 2 && (has_start ^ has_end))
              throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                "Could not dequote string '" + block + "' due to wrongly placed '\"'.");
            if (block.size() >= 2 && has_start && has_end)
              block = block.substr(1, block.size() - 2);
            substrings.push_back(std::move(block));
            begin = end + 1;
          }
        }
        if (substrings.empty()) { substrings.push_back(s); return false; }
        std::string block(begin, end);
        trim(block);
        bool has_start = block.size() >= 1 && block.front() == '"';
        bool has_end   = block.size() >= 1 && block.back()  == '"';
        if (block.size() >= 2 && (has_start ^ has_end))
          throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Could not dequote string '" + block + "' due to wrongly placed '\"'.");
        if (block.size() >= 2 && has_start && has_end)
          block = block.substr(1, block.size() - 2);
        substrings.push_back(std::move(block));
      }
      else
      {
        auto begin = s.cbegin();
        for (auto it = s.cbegin(); it != s.cend(); ++it)
        {
          if (*it == splitter)
          {
            substrings.emplace_back(begin, it);
            begin = it + 1;
          }
        }
        substrings.emplace_back(begin, s.cend());
      }
      return true;
    }

    /// Split @p s on multi-character @p splitter string.
    /// Returns true if at least one split occurred.
    inline bool split(const std::string& s, const std::string& splitter,
                      std::vector<std::string>& substrings)
    {
      substrings.clear();
      if (s.empty()) return false;
      if (splitter.empty())
      {
        substrings.resize(s.size());
        for (size_t i = 0; i < s.size(); ++i) substrings[i] = std::string(1, s[i]);
        return true;
      }
      size_t len   = splitter.size();
      size_t start = 0;
      size_t pos   = s.find(splitter);
      while (pos != std::string::npos)
      {
        substrings.push_back(s.substr(start, pos - start));
        start = pos + len;
        pos   = s.find(splitter, start);
      }
      substrings.push_back(s.substr(start));
      return substrings.size() > 1;
    }

    /// Split @p s on @p splitter, but do not split within quoted substrings.
    /// Returns true if at least one split occurred.
    inline bool split_quoted(const std::string& s, const std::string& splitter,
                             std::vector<std::string>& substrings,
                             char q = '"', QuotingMethod method = QuotingMethod::ESCAPE)
    {
      substrings.clear();
      if (s.empty() || splitter.empty()) return false;

      bool   in_quote = false;
      char   targets[2] = {q, splitter[0]};
      std::string rest = splitter.substr(1);
      size_t start = 0;

      for (size_t i = 0; i < s.size(); ++i)
      {
        if (in_quote)
        {
          bool embedded = false;
          if (method == QuotingMethod::ESCAPE)
          {
            for (; i < s.size(); ++i)
            {
              if      (s[i] == '\\') embedded = !embedded;
              else if (s[i] == q && !embedded) break;
              else embedded = false;
            }
          }
          else
          {
            for (; i < s.size(); ++i)
            {
              if (s[i] == q)
              {
                if (method == QuotingMethod::NONE) break;
                if (i + 1 < s.size() && s[i + 1] == q) embedded = !embedded;
                else if (!embedded) break;
                else embedded = false;
              }
            }
          }
          in_quote = false;
        }
        else
        {
          i = s.find_first_of(targets, i, 2);
          if (i == std::string::npos) break;
          if (s[i] == q) { in_quote = true; }
          else if (s.compare(i + 1, rest.size(), rest) == 0)
          {
            substrings.push_back(s.substr(start, i - start));
            start = i + splitter.size();
            i     = start - 1;
          }
        }
      }
      if (in_quote)
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "unbalanced quotation marks in string '" + s + "'");
      substrings.push_back(s.substr(start));
      return substrings.size() > 1;
    }

    /// Join elements [first, last) with @p glue between them, storing result in @p target
    template <class StringIterator>
    inline void concatenate(std::string& target,
                            StringIterator first, StringIterator last,
                            const std::string& glue = "")
    {
      if (first == last) { target.clear(); return; }
      target = *first;
      for (auto it = ++first; it != last; ++it)
      {
        target += glue;
        target += *it;
      }
    }

    /// Join elements [first, last) with @p glue and return the result
    template <class StringIterator>
    inline std::string concatenate(StringIterator first, StringIterator last,
                                   const std::string& glue = "")
    {
      std::string result;
      concatenate(result, first, last, glue);
      return result;
    }

    /**
      @brief Join the string elements of a @p container with @p glue between them and return the result.

      Convenience overload of concatenate(first, last, glue) for any container exposing
      begin()/end() (e.g. std::vector<std::string>, std::set<std::string>), so callers can write
      @code std::string s = StringUtils::concatenate(container, ", "); @endcode
      instead of the iterator form
      @code std::string s; s = StringUtils::concatenate(container, ", "); @endcode

      Constraining @p Container to a range keeps this from competing with the iterator overloads.
    */
    template <std::ranges::range Container>
    inline std::string concatenate(const Container& container, const std::string& glue = "")
    {
      return concatenate(container.begin(), container.end(), glue);
    }

  } // namespace StringUtils


  // =========================================================================
  // Legacy StringConversions namespace — thin forwarding layer
  // Prefer StringUtils:: for new code.
  // =========================================================================
  namespace StringConversions
  {
    template <typename T>
    inline void append(const T& i, std::string& target)
    { StringUtils::appendToStr(i, target); }

    inline void append(float f, std::string& target)
    { StringUtils::appendToStr(f, target); }
    inline void append(double d, std::string& target)
    { StringUtils::appendToStr(d, target); }
    inline void append(long double ld, std::string& target)
    { StringUtils::appendToStr(ld, target); }
    inline void appendLowP(float f, std::string& target)
    { StringUtils::appendToStrLowP(f, target); }
    inline void appendLowP(double d, std::string& target)
    { StringUtils::appendToStrLowP(d, target); }
    inline void appendLowP(long double ld, std::string& target)
    { StringUtils::appendToStrLowP(ld, target); }
    inline void append(const DataValue& d, bool full_precision, std::string& target)
    { StringUtils::appendToStr(d, full_precision, target); }

    template <typename T>
    inline std::string toString(const T& i)
    { return StringUtils::toStr(i); }
    inline std::string toString(float f, bool fp = true)
    { return StringUtils::toStr(f, fp); }
    inline std::string toString(double d, bool fp = true)
    { return StringUtils::toStr(d, fp); }
    inline std::string toString(long double ld, bool fp = true)
    { return StringUtils::toStr(ld, fp); }
    inline std::string toString(const DataValue& d, bool fp = true)
    { return StringUtils::toStr(d, fp); }
    inline std::string toString()
    { return {}; }

    inline std::string toStringLowP(float f)   { return StringUtils::toStr(f, false); }
    inline std::string toStringLowP(double d)   { return StringUtils::toStr(d, false); }
    inline std::string toStringLowP(long double ld) { return StringUtils::toStr(ld, false); }

  } // namespace StringConversions

} // namespace OpenMS


// =============================================================================
// Global operator+ / operator+= for std::string + numeric types
//
// Placing these in the global namespace ensures they are found by unqualified
// lookup regardless of which OpenMS namespace the caller is in.
// Using pass-by-value for operator+ enables move semantics for rvalue strings.
// =============================================================================

inline std::string& operator+=(std::string& s, int i)
{ OpenMS::StringUtils::appendToStr(i, s); return s; }
inline std::string& operator+=(std::string& s, unsigned int i)
{ OpenMS::StringUtils::appendToStr(i, s); return s; }
inline std::string& operator+=(std::string& s, short int i)
{ OpenMS::StringUtils::appendToStr(i, s); return s; }
inline std::string& operator+=(std::string& s, short unsigned int i)
{ OpenMS::StringUtils::appendToStr(i, s); return s; }
inline std::string& operator+=(std::string& s, long int i)
{ OpenMS::StringUtils::appendToStr(i, s); return s; }
inline std::string& operator+=(std::string& s, long unsigned int i)
{ OpenMS::StringUtils::appendToStr(i, s); return s; }
inline std::string& operator+=(std::string& s, long long unsigned int i)
{ OpenMS::StringUtils::appendToStr(i, s); return s; }
inline std::string& operator+=(std::string& s, long long signed int i)
{ OpenMS::StringUtils::appendToStr(i, s); return s; }
inline std::string& operator+=(std::string& s, float f)
{ OpenMS::StringUtils::appendToStr(f, s); return s; }
inline std::string& operator+=(std::string& s, double d)
{ OpenMS::StringUtils::appendToStr(d, s); return s; }
inline std::string& operator+=(std::string& s, long double ld)
{ OpenMS::StringUtils::appendToStr(ld, s); return s; }

// operator+ reuses operator+= on a value-copy, enabling move-from-rvalue
inline std::string operator+(std::string s, int i)
{ return std::move(s += i); }
inline std::string operator+(std::string s, unsigned int i)
{ return std::move(s += i); }
inline std::string operator+(std::string s, short int i)
{ return std::move(s += i); }
inline std::string operator+(std::string s, short unsigned int i)
{ return std::move(s += i); }
inline std::string operator+(std::string s, long int i)
{ return std::move(s += i); }
inline std::string operator+(std::string s, long unsigned int i)
{ return std::move(s += i); }
inline std::string operator+(std::string s, long long unsigned int i)
{ return std::move(s += i); }
inline std::string operator+(std::string s, long long signed int i)
{ return std::move(s += i); }
inline std::string operator+(std::string s, float f)
{ return std::move(s += f); }
inline std::string operator+(std::string s, double d)
{ return std::move(s += d); }
inline std::string operator+(std::string s, long double ld)
{ return std::move(s += ld); }
