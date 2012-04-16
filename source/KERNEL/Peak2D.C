// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Stephan Aiche$
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
