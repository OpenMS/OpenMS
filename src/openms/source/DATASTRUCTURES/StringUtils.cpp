// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg, Chris Bielow $
// $Authors: Marc Sturm, Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/StringUtils.h>

// DataValue can now be included here; it depends on String.h which includes StringUtils.h,
// but since this is a .cpp there is no circularity at compile time.
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/DATASTRUCTURES/ParamValue.h>
#include <OpenMS/SYSTEM/SIMDe.h>

namespace OpenMS
{
  // -------------------------------------------------------------------------
  // Static Boost.Spirit.Qi parser instances
  // -------------------------------------------------------------------------
  boost::spirit::qi::real_parser<double, StringUtilsHelper::real_policies_NANfixed_<double>>
    StringUtilsHelper::parse_double_{};
  boost::spirit::qi::real_parser<float, StringUtilsHelper::real_policies_NANfixed_<float>>
    StringUtilsHelper::parse_float_{};
  boost::spirit::qi::int_parser<>
    StringUtilsHelper::parse_int_{};


  // -------------------------------------------------------------------------
  // StringUtilsHelper — Boost.Spirit.Qi parsing implementations
  // -------------------------------------------------------------------------

  Int32 StringUtilsHelper::toInt32(const std::string& s)
  {
    Int32 ret;
    auto it = s.cbegin();
    if (!boost::spirit::qi::phrase_parse(it, s.cend(),
          boost::spirit::qi::int_, boost::spirit::ascii::space, ret))
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Could not convert string '" + s + "' to an integer value");
    }
    if (it != s.cend())
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Prefix of string '" + s + "' successfully converted to an int32 value. "
        "Additional characters found at position " +
        std::to_string(static_cast<int>(std::distance(s.cbegin(), it) + 1)));
    }
    return ret;
  }

  Int64 StringUtilsHelper::toInt64(const std::string& s)
  {
    Int64 ret;
    auto it = s.cbegin();
    if (!boost::spirit::qi::phrase_parse(it, s.cend(),
          boost::spirit::qi::long_long, boost::spirit::ascii::space, ret))
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Could not convert string '" + s + "' to an int64 value");
    }
    if (it != s.cend())
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Prefix of string '" + s + "' successfully converted to an integer value. "
        "Additional characters found at position " +
        std::to_string(static_cast<int>(std::distance(s.cbegin(), it) + 1)));
    }
    return ret;
  }

  float StringUtilsHelper::toFloat(const std::string& s)
  {
    float ret;
    auto it = s.cbegin();
    if (!boost::spirit::qi::phrase_parse(it, s.cend(),
          parse_float_, boost::spirit::ascii::space, ret))
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Could not convert string '" + s + "' to a float value");
    }
    if (it != s.cend())
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Prefix of string '" + s + "' successfully converted to a float value. "
        "Additional characters found at position " +
        std::to_string(static_cast<int>(std::distance(s.cbegin(), it) + 1)));
    }
    return ret;
  }

  double StringUtilsHelper::toDouble(const std::string& s)
  {
    double ret;
    auto it = s.cbegin();
    if (!boost::spirit::qi::phrase_parse(it, s.cend(),
          parse_double_, boost::spirit::ascii::space, ret))
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Could not convert string '" + s + "' to a double value");
    }
    if (it != s.cend())
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Prefix of string '" + s + "' successfully converted to a double value. "
        "Additional characters found at position " +
        std::to_string(static_cast<int>(std::distance(s.cbegin(), it) + 1)));
    }
    return ret;
  }


  // -------------------------------------------------------------------------
  // StringUtils — DataValue overloads (need DataValue.h, hence not inline)
  // -------------------------------------------------------------------------

  namespace StringUtils
  {
    void appendToStr(const DataValue& d, bool full_precision, std::string& target)
    {
      target += d.toString(full_precision);
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

} // namespace OpenMS
