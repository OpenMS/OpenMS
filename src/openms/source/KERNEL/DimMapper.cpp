// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/DimMapper.h>

#include <OpenMS/DATASTRUCTURES/String.h>

#include <charconv>
#include <cmath>

using namespace std;


namespace OpenMS
{
  namespace
  {
    DimMapper<1> dims({DIM_UNIT::RT});
    DimMapper<1> d({DIM_UNIT::RT});
    DimMapper<1> d2(d);
    bool x = (d == dims);
    Area<2> area(nullptr);

    /// Format a double with fixed precision and comma group separators (matching QLocale::c().toString(value, 'f', precision))
    std::string formatWithGroupSeparators(double value, int precision)
    {
      // Match QLocale::c().toString behavior for special floating-point values
      if (std::isnan(value))
        return "nan";
      if (std::isinf(value))
        return value < 0 ? "-inf" : "inf";

      char buf[512]; // large enough for any fixed-format double (DBL_MAX ~ 309 digits + precision)
      auto [ptr, ec] = std::to_chars(buf, buf + sizeof(buf), value, std::chars_format::fixed, precision);
      if (ec != std::errc{})
        return std::to_string(value); // fallback (should not happen with 512-byte buffer)
      *ptr = '\0';

      std::string result(buf);

      // Locate decimal point
      auto dot_pos = result.find('.');
      std::string integer_part = result.substr(0, dot_pos);
      std::string fractional_part;
      if (dot_pos != std::string::npos)
        fractional_part = result.substr(dot_pos);

      // Handle negative sign
      bool negative = (!integer_part.empty() && integer_part[0] == '-');
      if (negative)
        integer_part = integer_part.substr(1);

      // Insert comma group separators every 3 digits (C locale convention)
      int len = static_cast<int>(integer_part.size());
      if (len > 3)
      {
        std::string formatted;
        int first_group = len % 3;
        if (first_group == 0)
          first_group = 3;
        formatted = integer_part.substr(0, first_group);
        for (int i = first_group; i < len; i += 3)
        {
          formatted += ',';
          formatted += integer_part.substr(i, 3);
        }
        integer_part = formatted;
      }

      return std::string(negative ? "-" : "") + integer_part + fractional_part;
    }
  }

  String DimBase::formattedValue(const ValueType value) const
  {
    return String(this->getDimNameShort()) + ": " + String(formatWithGroupSeparators(value, valuePrecision()));
  }

  String DimBase::formattedValue(const ValueType value, const String& prefix) const
  {
    return prefix + formattedValue(value);
  }

  int DimBase::valuePrecision() const
  {
    // decide on precision depending on unit; add more units if you have some intuition
    constexpr auto precision_for_unit = [](DIM_UNIT u) {
      switch (u)
      {
        case DIM_UNIT::RT:
        case DIM_UNIT::INT:
          return 2;
        case DIM_UNIT::MZ:
          return 8;
        default:
          return 4;
      }
    };
    return precision_for_unit(this->getUnit());
  }
} // namespace OpenMS
