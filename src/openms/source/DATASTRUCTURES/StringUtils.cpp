// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg, Chris Bielow $
// $Authors: Marc Sturm, Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <cerrno>
#include <charconv>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iterator>
#include <memory>
#include <type_traits>

// DataValue can now be included here; it depends on String.h which includes StringUtils.h,
// but since this is a .cpp there is no circularity at compile time.
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/DATASTRUCTURES/ParamValue.h>
#include <OpenMS/SYSTEM/SIMDe.h>

namespace OpenMS
{
  // =========================================================================
  // StringUtilsHelper — parsing implementations (std::from_chars-based)
  // =========================================================================

  namespace
  {
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

    /// Append a double/float/long double via std::to_chars.
    /// NaN is output as "NaN" (uppercase) for backward compatibility, infinities as "inf"/"-inf".
    /// Trailing zeros are trimmed but at least one digit after '.' is kept (matches old karma behavior).
    template <typename T>
    inline void appendNumeric(T value, std::string& target, int precision, bool fixed_format)
    {
      if (std::isnan(value)) { target += "NaN"; return; }
      // Without this, std::to_chars writes "inf", which carries neither '.' nor 'e', so the
      // "keep at least one digit after the decimal point" rule below (there so that 5 prints as 5.0)
      // appended ".0" and produced "inf.0" - a token nothing can read back, which turned any file
      // containing an infinity into one that fails to load. "inf"/"-inf" round-trip through
      // toDouble()/toFloat() as they stand, so only the writing side was ever wrong.
      if (std::isinf(value)) { target += (value < T(0)) ? "-inf" : "inf"; return; }
      char buf[64];
      std::to_chars_result fc;

      // Determine format: use scientific for extreme values or when fixed_format is requested
      // but the value is too small/large for fixed notation
      T abs_val = (value < 0) ? -value : value;
      bool use_scientific = fixed_format
        ? (abs_val != T(0) && (abs_val >= T(1e4) || abs_val < T(1e-2)))
        : (abs_val != T(0) && (abs_val >= T(1e4) || abs_val < T(1e-2)));

      if (use_scientific)
      {
        if (fixed_format)
        {
          fc = std::to_chars(buf, buf + sizeof(buf), value, std::chars_format::scientific, 3);
        }
        else
        {
          // Use shortest round-trip representation (no explicit precision). The C++ standard
          // guarantees a unique shortest decimal that round-trips back to the same double, so
          // the result is platform-independent. A fixed precision instead rounds differently
          // across libc++/libstdc++ (the macOS CI divergence this fixes); shortest round-trip
          // reproduces the existing reference-file values without updates.
          fc = std::to_chars(buf, buf + sizeof(buf), value, std::chars_format::scientific);
        }
      }
      else if (fixed_format)
      {
        fc = std::to_chars(buf, buf + sizeof(buf), value, std::chars_format::fixed, precision);
      }
      else
      {
        fc = std::to_chars(buf, buf + sizeof(buf), value, std::chars_format::fixed, precision);
      }

      if (fc.ec == std::errc{})
      {
        const char* end = fc.ptr;

        // For scientific notation: post-process to match expected format
        // std::to_chars outputs "1.234e+04" but we want "1.234e04"
        const char* e_pos = reinterpret_cast<const char*>(std::memchr(buf, 'e', end - buf));
        if (e_pos)
        {
          const char* e_orig = e_pos; // save before mantissa trimming

          // Trim trailing zeros from mantissa (before 'e')
          const char* dot = reinterpret_cast<const char*>(std::memchr(buf, '.', e_orig - buf));
          if (dot)
          {
            while (e_pos > dot + 1 && *(e_pos - 1) == '0') --e_pos;
            if (e_pos == dot + 1) e_pos = dot + 2; // keep at least one digit after dot
          }
          // Append trimmed mantissa, then 'e'. The shortest round-trip representation can
          // produce an integer mantissa with no decimal point (e.g. "1e+06"); restore the
          // historical "1.0e06" form by ensuring at least one digit after the decimal point.
          target.append(buf, static_cast<size_t>(e_pos - buf));
          if (!dot) target += ".0";
          target += 'e';

          // Fix exponent format: "+04" → "04" (remove '+', keep '-' and zero-padding).
          // std::to_chars uses printf %e style, so the exponent always has >= 2 digits and
          // negative exponents already carry their '-'; only the leading '+' must be removed.
          const char* exp_start = e_orig + 1; // skip 'e'
          if (exp_start < end && *exp_start == '+')
            ++exp_start; // skip '+'
          size_t exp_len = static_cast<size_t>(end - exp_start);
          target.append(exp_start, exp_len);
          return;
        }

        // For fixed notation: trim trailing zeros after decimal point
        const char* dot = reinterpret_cast<const char*>(std::memchr(buf, '.', end - buf));
        if (dot)
        {
          while (end > dot + 1 && *(end - 1) == '0') --end;
          if (end == dot + 1) end = dot + 2; // keep at least one digit after dot
        }
        else if constexpr (std::is_floating_point_v<T>)
        {
          target.append(buf, static_cast<size_t>(end - buf));
          target += ".0";
          return;
        }
        target.append(buf, static_cast<size_t>(end - buf));
      }
      else
      {
        target += std::to_string(static_cast<double>(value));
      }
    }

    // libc++ (e.g. Apple Clang) does not implement the floating-point overloads of
    // std::from_chars — only the integer ones. Detect this and provide a strtod/strtof
    // based fallback with matching semantics (general format, no leading whitespace/'+').
    // All other standard libraries (libstdc++ >= 11, MSVC) use std::from_chars directly.
#if defined(_LIBCPP_VERSION)
  #define OPENMS_NO_FLOAT_FROM_CHARS 1
#endif

#ifdef OPENMS_NO_FLOAT_FROM_CHARS
    template <typename T>
    std::from_chars_result fromCharsFloat(const char* first, const char* last, T& value)
    {
      std::from_chars_result res{first, std::errc::invalid_argument};
      // std::from_chars(general) rejects leading whitespace and a leading '+'; mimic that
      // so behaviour is identical to the std::from_chars path on other platforms.
      if (first == last || *first == '+' || *first == ' ' || *first == '\t' ||
          *first == '\n' || *first == '\r' || *first == '\f' || *first == '\v')
      {
        return res;
      }

      // strtod/strtof need a null-terminated string and may read past 'last', so copy the
      // range into a bounded buffer first.
      const std::size_t n = static_cast<std::size_t>(last - first);
      char stackbuf[64];
      std::unique_ptr<char[]> heapbuf;
      char* buf = stackbuf;
      if (n + 1 > sizeof(stackbuf))
      {
        heapbuf = std::make_unique<char[]>(n + 1);
        buf = heapbuf.get();
      }
      std::memcpy(buf, first, n);
      buf[n] = '\0';

      errno = 0;
      char* parse_end = nullptr;
      if constexpr (std::is_same_v<T, float>)
        value = std::strtof(buf, &parse_end);
      else
        value = std::strtod(buf, &parse_end);

      if (parse_end == buf) // nothing consumed
      {
        return res; // invalid_argument, ptr == first
      }
      res.ptr = first + (parse_end - buf);
      // strtod/strtof set errno==ERANGE on both overflow (returns +/-HUGE_VAL) and
      // underflow (returns a representable subnormal or 0). std::from_chars only reports
      // result_out_of_range on overflow, so match that: accept underflow as success.
      res.ec = (errno == ERANGE && std::isinf(value)) ? std::errc::result_out_of_range : std::errc{};
      return res;
    }
#endif

    /// Thin wrapper dispatching to std::from_chars or the libc++ fallback above.
    template <typename T>
    inline std::from_chars_result parseFloat(const char* first, const char* last, T& value)
    {
#ifdef OPENMS_NO_FLOAT_FROM_CHARS
      return fromCharsFloat(first, last, value);
#else
      return std::from_chars(first, last, value);
#endif
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

    if (f != l && *f == '+') ++f;

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

    if (f != l && *f == '+') ++f;

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
    if (start != l && *start == '+') ++start;
    auto fc = parseFloat(start, l, ret);
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
    if (start != l && *start == '+') ++start;
    auto fc = parseFloat(start, l, ret);
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

    const char* start = begin;
    if (start != end && *start == '+') ++start;

    auto fc = parseFloat(start, end, target);
    if (fc.ec == std::errc{})
    {
      begin = fc.ptr;
      return true;
    }
    return false;
  }

  bool StringUtilsHelper::extractInt(const char*& begin, const char* end, int& target)
  {
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
  // StringUtils — appendToStr implementations (std::to_chars-based)
  //
  // Precision mapping (matches previous karma policies):
  //   float:   6 significant digits (general) or 3 fractional digits (low precision)
  //   double:  15 significant digits (general) or 3 fractional digits (low precision)
  //   ld:      18 significant digits (general) or 3 fractional digits (low precision)
  //
  // NaN output: "NaN" (uppercase) for backward compatibility with karma.
  // =========================================================================

  namespace StringUtils
  {
    // Integer overloads
    void appendToStr(int i, std::string& target)
    {
      char buf[32];
      auto [p, ec] = std::to_chars(buf, buf + sizeof(buf), i);
      if (ec == std::errc{}) target.append(buf, p);
    }
    void appendToStr(unsigned int i, std::string& target)
    {
      char buf[32];
      auto [p, ec] = std::to_chars(buf, buf + sizeof(buf), i);
      if (ec == std::errc{}) target.append(buf, p);
    }
    void appendToStr(short int i, std::string& target)
    {
      char buf[32];
      auto [p, ec] = std::to_chars(buf, buf + sizeof(buf), i);
      if (ec == std::errc{}) target.append(buf, p);
    }
    void appendToStr(short unsigned int i, std::string& target)
    {
      char buf[32];
      auto [p, ec] = std::to_chars(buf, buf + sizeof(buf), i);
      if (ec == std::errc{}) target.append(buf, p);
    }
    void appendToStr(long int i, std::string& target)
    {
      char buf[32];
      auto [p, ec] = std::to_chars(buf, buf + sizeof(buf), i);
      if (ec == std::errc{}) target.append(buf, p);
    }
    void appendToStr(long unsigned int i, std::string& target)
    {
      char buf[32];
      auto [p, ec] = std::to_chars(buf, buf + sizeof(buf), i);
      if (ec == std::errc{}) target.append(buf, p);
    }
    void appendToStr(long long unsigned int i, std::string& target)
    {
      char buf[32];
      auto [p, ec] = std::to_chars(buf, buf + sizeof(buf), i);
      if (ec == std::errc{}) target.append(buf, p);
    }
    void appendToStr(long long signed int i, std::string& target)
    {
      char buf[32];
      auto [p, ec] = std::to_chars(buf, buf + sizeof(buf), i);
      if (ec == std::errc{}) target.append(buf, p);
    }

    // Float overloads — 6 significant digits (general) / 3 fractional digits (fixed)
    void appendToStr(float f, std::string& target)
    {
      appendNumeric(f, target, writtenDigits<float>(), false);
    }
    void appendToStrLowP(float f, std::string& target)
    {
      appendNumeric(f, target, 3, true);
    }

    // Double overloads — 15 significant digits / 3 fractional digits
    void appendToStr(double d, std::string& target)
    {
      appendNumeric(d, target, writtenDigits<double>(), false);
    }
    void appendToStrLowP(double d, std::string& target)
    {
      appendNumeric(d, target, 3, true);
    }

    // Long double overloads
    void appendToStr(long double ld, std::string& target)
    {
      appendNumeric(ld, target, writtenDigits<long double>(), false);
    }
    void appendToStrLowP(long double ld, std::string& target)
    {
      appendNumeric(ld, target, 3, true);
    }

    // DataValue overload
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
