// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/DimMapper.h>

#include <OpenMS/DATASTRUCTURES/String.h>

#include <sstream>
#include <iomanip>

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
  }

  String DimBase::formattedValue(const ValueType value) const
  {
    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss << std::setprecision(valuePrecision()) << value;
    const auto sv = this->getDimNameShort();
    return String(sv.data(), static_cast<String::SizeType>(sv.size())) + ": " + String(oss.str());
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
