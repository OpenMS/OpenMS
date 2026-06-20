// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg, Chris Bielow $
// $Authors: Marc Sturm, Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/home/karma/numeric/detail/real_utils.hpp>
#include <boost/type_traits.hpp>

#include <charconv>
#include <cstdlib>
#include <ctime>
#include <iterator>

// DataValue can now be included here; it depends on String.h which includes StringUtils.h,
// but since this is a .cpp there is no circularity at compile time.
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/DATASTRUCTURES/ParamValue.h>
#include <OpenMS/SYSTEM/SIMDe.h>

namespace OpenMS
{
  // =========================================================================
  // Boost.Karma generator policies (precision, NaN handling) — internal
  // =========================================================================
  namespace
  {
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

    template<typename T>
    class BK_PrecPolicyShort : public boost::spirit::karma::real_policies<T>
    {
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

    const FloatFullGen   floatFull{};
    const DoubleFullGen  doubleFull{};
    const LDFullGen      ldFull{};
    const FloatShortGen  floatShort{};
    const DoubleShortGen doubleShort{};
    const LDShortGen     ldShort{};
  } // anonymous namespace


  // =========================================================================
  // StringUtilsHelper — parsing implementations (std::from_chars-based)
  // =========================================================================

  namespace
  {
    /// Try to parse a NaN variant ("nan", "NAN", "NaN(...)") starting at @p p.
    template <typename T>
    bool tryParseNaN(const char*& p, const char* end, T& target)
    {
      if (p == end) return false;
      if (*p != 'n' && *p != 'N') return false;

      const char* start = p;
      if ((end - start) >= 3 &&
          (start[0] == 'n' || start[0] == 'N') &&
          (start[1] == 'a' || start[1] == 'A') &&
          (start[2] == 'n' || start[2] == 'N'))
      {
        p += 3;
        if (p != end && *p == '(')
        {
          ++p;
          while (p != end && *p != ')') ++p;
          if (p == end) { p = start; return false; }
          ++p;
        }
        target = std::numeric_limits<T>::quiet_NaN();
        return true;
      }
      return false;
    }

    inline const char* skipWS(const char* p, const char* end)
    {
      while (p != end && (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r'))
        ++p;
      return p;
    }
  } // anonymous namespace


  Int32 StringUtilsHelper::toInt32(const std::string& s)
  {
    const char* f = s.data();
    const char* l = s.data() + s.size();

    f = skipWS(f, l);
    if (f == l)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Could not convert string '" + s + "' to an integer value");
    }

    Int32 ret{};
    auto fc = std::from_chars(f, l, ret);
    if (fc.ec != std::errc{})
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Could not convert string '" + s + "' to an integer value");
    }
    const char* after = skipWS(fc.ptr, l);
    if (after != l)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Prefix of string '" + s + "' successfully converted to an int32 value. "
        "Additional characters found at position " +
        std::to_string(static_cast<int>(std::distance(s.data(), after) + 1)));
    }
    return ret;
  }

  Int64 StringUtilsHelper::toInt64(const std::string& s)
  {
    const char* f = s.data();
    const char* l = s.data() + s.size();

    f = skipWS(f, l);
    if (f == l)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Could not convert string '" + s + "' to an integer value");
    }

    Int64 ret{};
    auto fc = std::from_chars(f, l, ret);
    if (fc.ec != std::errc{})
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Could not convert string '" + s + "' to an integer value");
    }
    const char* after = skipWS(fc.ptr, l);
    if (after != l)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Prefix of string '" + s + "' successfully converted to an integer value. "
        "Additional characters found at position " +
        std::to_string(static_cast<int>(std::distance(s.data(), after) + 1)));
    }
    return ret;
  }

  float StringUtilsHelper::toFloat(const std::string& s)
  {
    const char* f = s.data();
    const char* l = s.data() + s.size();

    f = skipWS(f, l);
    if (f == l)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Could not convert string '" + s + "' to a float value");
    }

    float nan_val{};
    if (tryParseNaN(f, l, nan_val))
    {
      const char* after = skipWS(f, l);
      if (after == l) return nan_val;
    }

    float ret{};
    const char* start = skipWS(s.data(), l);
    auto fc = std::from_chars(start, l, ret);
    if (fc.ec != std::errc{})
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Could not convert string '" + s + "' to a float value");
    }
    const char* after = skipWS(fc.ptr, l);
    if (after != l)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Prefix of string '" + s + "' successfully converted to a float value. "
        "Additional characters found at position " +
        std::to_string(static_cast<int>(std::distance(s.data(), after) + 1)));
    }
    return ret;
  }

  double StringUtilsHelper::toDouble(const std::string& s)
  {
    const char* f = s.data();
    const char* l = s.data() + s.size();

    f = skipWS(f, l);
    if (f == l)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Could not convert string '" + s + "' to a double value");
    }

    double nan_val{};
    if (tryParseNaN(f, l, nan_val))
    {
      const char* after = skipWS(f, l);
      if (after == l) return nan_val;
    }

    double ret{};
    const char* start = skipWS(s.data(), l);
    auto fc = std::from_chars(start, l, ret);
    if (fc.ec != std::errc{})
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Could not convert string '" + s + "' to a double value");
    }
    const char* after = skipWS(fc.ptr, l);
    if (after != l)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Prefix of string '" + s + "' successfully converted to a double value. "
        "Additional characters found at position " +
        std::to_string(static_cast<int>(std::distance(s.data(), after) + 1)));
    }
    return ret;
  }

  bool StringUtilsHelper::extractDouble(const char*& begin, const char* end, double& target)
  {
    if (tryParseNaN(begin, end, target))
      return true;

    // std::from_chars doesn't handle leading '+', but boost::spirit did
    const char* start = begin;
    if (start != end && *start == '+') ++start;

    auto fc = std::from_chars(start, end, target);
    if (fc.ec == std::errc{})
    {
      begin = fc.ptr;
      return true;
    }
    return false;
  }

  bool StringUtilsHelper::extractInt(const char*& begin, const char* end, int& target)
  {
    // std::from_chars doesn't handle leading '+', but boost::spirit did
    const char* start = begin;
    if (start != end && *start == '+') ++start;

    auto fc = std::from_chars(start, end, target);
    if (fc.ec == std::errc{})
    {
      begin = fc.ptr;
      return true;
    }
    return false;
  }


  // =========================================================================
  // StringUtils — appendToStr implementations (Boost.Karma)
  // =========================================================================

  namespace StringUtils
  {
    void appendToStr(int i, std::string& target)
    {
      std::back_insert_iterator<std::string> sink(target);
      boost::spirit::karma::generate(sink, i);
    }
    void appendToStr(unsigned int i, std::string& target)
    {
      std::back_insert_iterator<std::string> sink(target);
      boost::spirit::karma::generate(sink, i);
    }
    void appendToStr(short int i, std::string& target)
    {
      std::back_insert_iterator<std::string> sink(target);
      boost::spirit::karma::generate(sink, i);
    }
    void appendToStr(short unsigned int i, std::string& target)
    {
      std::back_insert_iterator<std::string> sink(target);
      boost::spirit::karma::generate(sink, i);
    }
    void appendToStr(long int i, std::string& target)
    {
      std::back_insert_iterator<std::string> sink(target);
      boost::spirit::karma::generate(sink, i);
    }
    void appendToStr(long unsigned int i, std::string& target)
    {
      std::back_insert_iterator<std::string> sink(target);
      boost::spirit::karma::generate(sink, i);
    }
    void appendToStr(long long unsigned int i, std::string& target)
    {
      std::back_insert_iterator<std::string> sink(target);
      boost::spirit::karma::generate(sink, i);
    }
    void appendToStr(long long signed int i, std::string& target)
    {
      std::back_insert_iterator<std::string> sink(target);
      boost::spirit::karma::generate(sink, i);
    }

    void appendToStr(float f, std::string& target)
    {
      std::back_insert_iterator<std::string> sink(target);
      boost::spirit::karma::generate(sink, floatFull, f);
    }
    void appendToStrLowP(float f, std::string& target)
    {
      std::back_insert_iterator<std::string> sink(target);
      boost::spirit::karma::generate(sink, floatShort, f);
    }

    void appendToStr(double d, std::string& target)
    {
      std::back_insert_iterator<std::string> sink(target);
      boost::spirit::karma::generate(sink, doubleFull, d);
    }
    void appendToStrLowP(double d, std::string& target)
    {
      std::back_insert_iterator<std::string> sink(target);
      boost::spirit::karma::generate(sink, doubleShort, d);
    }

    void appendToStr(long double ld, std::string& target)
    {
      std::back_insert_iterator<std::string> sink(target);
      boost::spirit::karma::generate(sink, ldFull, ld);
    }
    void appendToStrLowP(long double ld, std::string& target)
    {
      std::back_insert_iterator<std::string> sink(target);
      boost::spirit::karma::generate(sink, ldShort, ld);
    }

    void appendToStr(const DataValue& d, bool full_precision, std::string& target)
    {
      target += d.toString(full_precision);
    }

    std::string toStr(float f, bool full_precision)
    {
      std::string s;
      full_precision ? appendToStr(f, s) : appendToStrLowP(f, s);
      return s;
    }

    std::string toStr(double d, bool full_precision)
    {
      std::string s;
      full_precision ? appendToStr(d, s) : appendToStrLowP(d, s);
      return s;
    }

    std::string toStr(long double ld, bool full_precision)
    {
      std::string s;
      full_precision ? appendToStr(ld, s) : appendToStrLowP(ld, s);
      return s;
    }

    std::string toStr(const DataValue& d, bool full_precision)
    {
      return d.toString(full_precision);
    }

    std::string toStr(const ParamValue& p, bool full_precision)
    {
      return p.toString(full_precision);
    }
  } // namespace StringUtils


  // -------------------------------------------------------------------------
  // SIMD-accelerated whitespace scanners
  // -------------------------------------------------------------------------

  namespace StringUtils
  {
    const char* skipWhitespace(const char* p, const char* p_end)
    {
      const simde__m128i w0 = simde_mm_set1_epi8(' ');
      const simde__m128i w1 = simde_mm_set1_epi8('\t');
      const simde__m128i w2 = simde_mm_set1_epi8('\n');
      const simde__m128i w3 = simde_mm_set1_epi8('\r');

      for (; p <= p_end - 16; p += 16)
      {
        const simde__m128i s = simde_mm_loadu_si128(reinterpret_cast<const simde__m128i*>(p));
        simde__m128i x = simde_mm_cmpeq_epi8(s, w0);
        x = simde_mm_or_si128(x, simde_mm_cmpeq_epi8(s, w1));
        simde__m128i y = simde_mm_cmpeq_epi8(s, w2);
        y = simde_mm_or_si128(y, simde_mm_cmpeq_epi8(s, w3));
        x = simde_mm_or_si128(x, y);
        auto non_space = static_cast<uint16_t>(~simde_mm_movemask_epi8(x));
        if (non_space != 0)
        {
#ifdef _MSC_VER
          unsigned long offset;
          _BitScanForward(&offset, non_space);
          return p + offset;
#else
          return p + __builtin_ffs(non_space) - 1;
#endif
        }
      }
      while (p != p_end)
      {
        if (*p == ' ' || *p == '\n' || *p == '\r' || *p == '\t') ++p;
        else return p;
      }
      return p_end;
    }

    const char* skipNonWhitespace(const char* p, const char* p_end)
    {
      const simde__m128i w0 = simde_mm_set1_epi8(' ');
      const simde__m128i w1 = simde_mm_set1_epi8('\t');
      const simde__m128i w2 = simde_mm_set1_epi8('\n');
      const simde__m128i w3 = simde_mm_set1_epi8('\r');

      for (; p <= p_end - 16; p += 16)
      {
        const simde__m128i s = simde_mm_loadu_si128(reinterpret_cast<const simde__m128i*>(p));
        simde__m128i x = simde_mm_cmpeq_epi8(s, w0);
        x = simde_mm_or_si128(x, simde_mm_cmpeq_epi8(s, w1));
        simde__m128i y = simde_mm_cmpeq_epi8(s, w2);
        y = simde_mm_or_si128(y, simde_mm_cmpeq_epi8(s, w3));
        x = simde_mm_or_si128(x, y);
        auto spaces = static_cast<uint16_t>(simde_mm_movemask_epi8(x));
        if (spaces != 0)
        {
#ifdef _MSC_VER
          unsigned long offset;
          _BitScanForward(&offset, spaces);
          return p + offset;
#else
          return p + __builtin_ffs(spaces) - 1;
#endif
        }
      }
      while (p != p_end)
      {
        if (*p == ' ' || *p == '\n' || *p == '\r' || *p == '\t') return p;
        else ++p;
      }
      return p_end;
    }
  } // namespace StringUtils


  // -------------------------------------------------------------------------
  // StringUtils — factory methods
  // -------------------------------------------------------------------------

  namespace StringUtils
  {
    std::string number(double d, UInt n)
    {
      char buf[64];
      std::snprintf(buf, sizeof(buf), "%.*f", static_cast<int>(n), d);
      return std::string(buf);
    }

    std::string numberLength(double d, UInt n)
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

    std::string random(UInt length)
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
  } // namespace StringUtils

} // namespace OpenMS
