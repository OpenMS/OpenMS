// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/MobilityPeak2D.h>

namespace OpenMS
{
  char const * const MobilityPeak2D::dimension_name_short_[] =
  {
    "IM",
    "MZ"
  };

  char const * const MobilityPeak2D::dimension_name_full_[] =
  {
    "ion mobility",
    "mass-to-charge"
  };

  char const * const MobilityPeak2D::dimension_unit_short_[] =
  {
    "?",
    "Th"
  };

  char const * const MobilityPeak2D::dimension_unit_full_[] =
  {
    "?",
    "Thomson"
  };

  char const * MobilityPeak2D::shortDimensionName(UInt const dim)
  {
    return dimension_name_short_[dim];
  }

  char const * MobilityPeak2D::shortDimensionNameIM()
  {
    return dimension_name_short_[IM];
  }

  char const * MobilityPeak2D::shortDimensionNameMZ()
  {
    return dimension_name_short_[MZ];
  }

  char const * MobilityPeak2D::fullDimensionName(UInt const dim)
  {
    return dimension_name_full_[dim];
  }

  char const * MobilityPeak2D::fullDimensionNameIM()
  {
    return dimension_name_full_[IM];
  }

  char const * MobilityPeak2D::fullDimensionNameMZ()
  {
    return dimension_name_full_[MZ];
  }

  char const * MobilityPeak2D::shortDimensionUnit(UInt const dim)
  {
    return dimension_unit_short_[dim];
  }

  char const * MobilityPeak2D::shortDimensionUnitIM()
  {
    return dimension_unit_short_[IM];
  }

  char const * MobilityPeak2D::shortDimensionUnitMZ()
  {
    return dimension_unit_short_[MZ];
  }

  char const * MobilityPeak2D::fullDimensionUnit(UInt const dim)
  {
    return dimension_unit_full_[dim];
  }

  char const * MobilityPeak2D::fullDimensionUnitIM()
  {
    return dimension_unit_full_[IM];
  }

  char const * MobilityPeak2D::fullDimensionUnitMZ()
  {
    return dimension_unit_full_[MZ];
  }

  std::ostream & operator<<(std::ostream & os, const MobilityPeak2D & point)
  {
    os << "IM: " << point.getMobility() <<  " MZ: "  << point.getMZ() << " INT: " << point.getIntensity();
    return os;
  }

} // namespace OpenMS
