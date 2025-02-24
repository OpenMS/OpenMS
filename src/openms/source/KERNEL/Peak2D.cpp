// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/Peak2D.h>

namespace OpenMS
{
  char const * const Peak2D::dimension_name_short_[] =
  {
    "RT",
    "MZ"
  };

  char const * const Peak2D::dimension_name_full_[] =
  {
    "retention time",
    "mass-to-charge"
  };

  char const * const Peak2D::dimension_unit_short_[] =
  {
    "sec",
    "Th"
  };

  char const * const Peak2D::dimension_unit_full_[] =
  {
    "Seconds",
    "Thomson"
  };

  char const * Peak2D::shortDimensionName(UInt const dim)
  {
    return dimension_name_short_[dim];
  }

  char const * Peak2D::shortDimensionNameRT()
  {
    return dimension_name_short_[RT];
  }

  char const * Peak2D::shortDimensionNameMZ()
  {
    return dimension_name_short_[MZ];
  }

  char const * Peak2D::fullDimensionName(UInt const dim)
  {
    return dimension_name_full_[dim];
  }

  char const * Peak2D::fullDimensionNameRT()
  {
    return dimension_name_full_[RT];
  }

  char const * Peak2D::fullDimensionNameMZ()
  {
    return dimension_name_full_[MZ];
  }

  char const * Peak2D::shortDimensionUnit(UInt const dim)
  {
    return dimension_unit_short_[dim];
  }

  char const * Peak2D::shortDimensionUnitRT()
  {
    return dimension_unit_short_[RT];
  }

  char const * Peak2D::shortDimensionUnitMZ()
  {
    return dimension_unit_short_[MZ];
  }

  char const * Peak2D::fullDimensionUnit(UInt const dim)
  {
    return dimension_unit_full_[dim];
  }

  char const * Peak2D::fullDimensionUnitRT()
  {
    return dimension_unit_full_[RT];
  }

  char const * Peak2D::fullDimensionUnitMZ()
  {
    return dimension_unit_full_[MZ];
  }

  std::ostream & operator<<(std::ostream & os, const Peak2D & point)
  {
    os << "RT: " << point.getRT() <<  " MZ: "  << point.getMZ() << " INT: " << point.getIntensity();
    return os;
  }

} // namespace OpenMS
