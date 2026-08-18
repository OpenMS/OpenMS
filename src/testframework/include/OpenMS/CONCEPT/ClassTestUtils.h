// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

// Standard-library-only utilities owned by the class-test framework.
//
// The framework deliberately depends on nothing but the C++ standard library, so
// that the tests of every OpenMS library (including ones that do not link
// libOpenMS, like OpenSwathAlgo) can use it. The few utilities it needs are
// therefore implemented here rather than taken from libOpenMS. Two of them
// duplicate libOpenMS behavior on purpose and must not drift:
//
//  * detail::appendNumeric mirrors OpenMS::StringUtils' appendNumeric
//    (src/openms/source/DATASTRUCTURES/StringUtils.cpp) so that
//    TEST_EQUAL(std::string, <number>) stringifies numbers exactly like
//    StringUtils::toStr() does. The ~720-test suite is the oracle for this.
//  * writtenDigits mirrors OpenMS::writtenDigits (OpenMS/CONCEPT/Types.h) for
//    the types TEST_REAL_SIMILAR is used with.

#include <charconv>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <limits>
#include <ostream>
#include <ranges>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

namespace OpenMS
{
  namespace Internal
  {
    namespace ClassTest
    {
      /**
        @brief Number of decimal digits to print for a value in TEST_REAL_SIMILAR.

        The framework's own version of OpenMS::writtenDigits (OpenMS/CONCEPT/Types.h),
        so that tests of libraries that do not link libOpenMS can use TEST_REAL_SIMILAR.
        Class types convertible to double (e.g. DataValue, ParamValue) are printed
        like a double, exactly as the former writtenDigits<DataValue> specialization did.
      */
      template <typename T>
      constexpr int writtenDigits(const T& = T())
      {
        using C = std::remove_cvref_t<T>;
        if constexpr (std::is_floating_point_v<C>)
        {
          return std::numeric_limits<C>::digits10;
        }
        else if constexpr (std::is_integral_v<C>)
        {
          return std::numeric_limits<C>::digits10;
        }
        else if constexpr (std::is_class_v<C> && std::is_convertible_v<const C&, double>)
        {
          return std::numeric_limits<double>::digits10;
        }
        else
        {
          return 6; // default precision of an ostream (27.4.4.1), like OpenMS::writtenDigits
        }
      }

      /**
        @brief Whether TEST_REAL_SIMILAR may convert a value of this type to double.

        True for floating-point types and for class types implicitly convertible to
        double (e.g. DataValue, ParamValue) -- the framework does not need to know
        those types by name. Anything else (integers, strings, ...) is rejected at
        runtime with "does not have a floating point type", as before.
      */
      template <typename T>
      constexpr bool isRealType(const T&)
      {
        using C = std::remove_cvref_t<T>;
        return std::is_floating_point_v<C> || (std::is_class_v<C> && std::is_convertible_v<const C&, double>);
      }

      namespace detail
      {
        /// Remove leading and trailing whitespace (' ', '\\t', '\\n', '\\r'), like StringUtils::trim.
        inline void trim(std::string& s)
        {
          const char* ws = " \t\n\r";
          const size_t b = s.find_first_not_of(ws);
          if (b == std::string::npos)
          {
            s.clear();
            return;
          }
          const size_t e = s.find_last_not_of(ws);
          s = s.substr(b, e - b + 1);
        }

        /// Split on a single character. Mirrors ListUtils::create<std::string>(s, splitter):
        /// an empty input gives an empty result; otherwise empty fields are kept and
        /// fields are NOT trimmed.
        inline std::vector<std::string> split(const std::string& s, char splitter)
        {
          std::vector<std::string> out;
          if (s.empty()) return out;
          size_t begin = 0;
          for (size_t i = 0; i < s.size(); ++i)
          {
            if (s[i] == splitter)
            {
              out.emplace_back(s, begin, i - begin);
              begin = i + 1;
            }
          }
          out.emplace_back(s, begin, s.size() - begin);
          return out;
        }

        /// Join with a separator (display purposes only).
        inline std::string join(const std::vector<std::string>& v, const char* sep)
        {
          std::string out;
          for (size_t i = 0; i < v.size(); ++i)
          {
            if (i) out += sep;
            out += v[i];
          }
          return out;
        }

        /// Append a double/float/long double via std::to_chars.
        /// VERBATIM copy of appendNumeric in src/openms/source/DATASTRUCTURES/StringUtils.cpp --
        /// TEST_EQUAL(std::string, <number>) must format numbers exactly like
        /// StringUtils::toStr(). Keep the two in sync; do not "improve" this one.
        /// NaN is output as "NaN" (uppercase), infinities as "inf"/"-inf".
        /// Trailing zeros are trimmed but at least one digit after '.' is kept.
        template <typename T>
        inline void appendNumeric(T value, std::string& target, int precision, bool fixed_format)
        {
          if (std::isnan(value)) { target += "NaN"; return; }
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
              // Shortest round-trip representation (no explicit precision): the unique shortest
              // decimal that round-trips back to the same double, hence platform-independent.
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

            // For scientific notation: post-process "1.234e+04" into the historical "1.234e04"
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
              target.append(buf, static_cast<size_t>(e_pos - buf));
              if (!dot) target += ".0";
              target += 'e';

              // Fix exponent format: "+04" -> "04" (remove '+', keep '-' and zero-padding).
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

        ///@name The stringification used by TEST_EQUAL(std::string, x).
        /// Mirrors the StringUtils::toStr overload set (which mirrored the converting
        /// constructors of the removed OpenMS::String), so e.g. TEST_EQUAL(some_string, 114)
        /// keeps comparing against "114" and doubles keep their full-precision format.
        //@{
        inline std::string toString(const std::string& s) { return s; }
        inline std::string toString(std::string_view sv) { return std::string(sv); }
        inline std::string toString(const char* s) { return std::string(s); }
        inline std::string toString(char c) { return std::string(1, c); }
        template <typename T>
        inline std::string toString(const T& v)
        {
          static_assert(std::is_arithmetic_v<T>,
                        "TEST_EQUAL(std::string, x): x must be a string, character or number. "
                        "For other types, stringify explicitly in the test.");
          std::string r;
          if constexpr (std::is_same_v<T, bool>)
          {
            r = v ? "1" : "0";
          }
          else if constexpr (std::is_floating_point_v<T>)
          {
            appendNumeric(v, r, std::numeric_limits<T>::digits10, false); // full precision, like toStr(v, true)
          }
          else
          {
            char buf[32];
            auto [p, ec] = std::to_chars(buf, buf + sizeof(buf), v);
            if (ec == std::errc{}) r.append(buf, p);
          }
          return r;
        }
        //@}

        /// Is 'os << t' well-formed? (Checked in this header's context plus ADL, so a
        /// vector<T> counts as streamable only via ADL-found operators, e.g. for T in
        /// namespace OpenMS with <OpenMS/DATASTRUCTURES/ListUtilsIO.h> included.)
        template <typename T>
        concept Streamable = requires(std::ostream& os, const T& t) { os << t; };

        /// Print a value in a TEST_EQUAL/TEST_NOT_EQUAL failure report: directly if
        /// streamable; ranges of printable elements as "[e1, e2, ...]" (the format
        /// OpenMS/DATASTRUCTURES/ListUtilsIO.h prints vectors in); otherwise a placeholder.
        template <typename T>
        void printValue(std::ostream& os, const T& v)
        {
          using C = std::remove_cvref_t<T>;
          if constexpr (Streamable<C>)
          {
            os << v;
          }
          else if constexpr (std::ranges::input_range<C>)
          {
            os << '[';
            bool first = true;
            for (const auto& e : v)
            {
              if (!first) os << ", ";
              first = false;
              using E = std::remove_cvref_t<decltype(e)>;
              if constexpr (std::is_arithmetic_v<E>)
              {
                os << toString(e); // numbers like ListUtilsIO (via toStr)
              }
              else if constexpr (Streamable<E>)
              {
                os << e;
              }
              else
              {
                os << "<?>";
              }
            }
            os << ']';
          }
          else
          {
            os << "<value of unprintable type - provide operator<< or compare via TEST_TRUE>";
          }
        }

        /// Convert a UTF-8 std::string to std::filesystem::path safely on all platforms.
        /// Copy of OpenMS::to_path (OpenMS/SYSTEM/PathUtils.h), kept here so the framework
        /// needs no OpenMS header. On Windows, std::filesystem::path(std::string) uses the
        /// current code page, not UTF-8, so construct from u8string; if the bytes are not
        /// valid UTF-8 (e.g. a filename from argv in the ANSI code page), fall back.
        inline std::filesystem::path to_path(const std::string& s)
        {
#ifdef _WIN32
          try
          {
            return std::filesystem::path(std::u8string(reinterpret_cast<const char8_t*>(s.data()), s.size()));
          }
          catch (const std::system_error&)
          {
            return std::filesystem::path(s);
          }
#else
          return std::filesystem::path(std::u8string(reinterpret_cast<const char8_t*>(s.data()), s.size()));
#endif
        }
      } // namespace detail

      /**
        @brief Stream wrapper printing a floating point value with full precision.

        The framework's own version of OpenMS::precisionWrapper (OpenMS/CONCEPT/PrecisionWrapper.h):
        printing goes through detail::toString, i.e. the same format as StringUtils::toStr(x, true).
        Used by TEST_FILE_SIMILAR's report.
      */
      template <typename FloatingPointType>
      struct PrecisionWrapped
      {
        FloatingPointType value;
      };

      template <typename FloatingPointType>
      inline PrecisionWrapped<FloatingPointType> precisionWrapper(const FloatingPointType value)
      {
        return PrecisionWrapped<FloatingPointType> {value};
      }

      template <typename FloatingPointType>
      inline std::ostream& operator<<(std::ostream& os, const PrecisionWrapped<FloatingPointType>& p)
      {
        os << detail::toString(p.value);
        return os;
      }
    } // namespace ClassTest
  } // namespace Internal
} // namespace OpenMS
