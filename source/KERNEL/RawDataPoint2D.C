// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/RawDataPoint2D.h>

namespace OpenMS
{

	char const * const RawDataPoint2D::dimension_name_short [] =
    {
      "RT",
      "MZ"
    };
	
  char const * const RawDataPoint2D::dimension_name_full [] =
    {
      "retention time",
      "mass-to-charge"
    };
  
  char const * const RawDataPoint2D::dimension_unit_short [] =
    {
      "sec",
      "Th"
    };
  
  char const * const RawDataPoint2D::dimension_unit_full [] =
    {
      "Seconds",
      "Thomson"
    };
  
	char const * const RawDataPoint2D::shortDimensionName(DimensionId const dim)
	{
		return dimension_name_short[dim];
	}

	char const * const RawDataPoint2D::shortDimensionNameRT()
	{
		return dimension_name_short[RT];
	}

	char const * const RawDataPoint2D::shortDimensionNameMZ()
	{
		return dimension_name_short[MZ];
	}

	char const * const RawDataPoint2D::fullDimensionName(DimensionId const dim)
	{
		return dimension_name_full[dim];
	}

	char const * const RawDataPoint2D::fullDimensionNameRT()
	{
		return dimension_name_full[RT];
	}

	char const * const RawDataPoint2D::fullDimensionNameMZ()
	{
		return dimension_name_full[MZ];
	}

	char const * const RawDataPoint2D::shortDimensionUnit(DimensionId const dim)
	{
		return dimension_unit_short[dim];
	}

	char const * const RawDataPoint2D::shortDimensionUnitRT()
	{
		return dimension_unit_short[RT];
	}

	char const * const RawDataPoint2D::shortDimensionUnitMZ()
	{
		return dimension_unit_short[MZ];
	}

	char const * const RawDataPoint2D::fullDimensionUnit(DimensionId const dim)
	{
		return dimension_unit_full[dim];
	}

	char const * const RawDataPoint2D::fullDimensionUnitRT()
	{
		return dimension_unit_full[RT];
	}

	char const * const RawDataPoint2D::fullDimensionUnitMZ()
	{
		return dimension_unit_full[MZ];
	}

	///Print the contents to a stream.
	std::ostream& operator << (std::ostream& os, const RawDataPoint2D& point)
	{
		os << "POS: "<< point.getMZ() << " INT: "<<point.getIntensity();
		
		return os;
	}
}
