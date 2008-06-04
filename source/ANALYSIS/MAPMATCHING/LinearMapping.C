// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/LinearMapping.h>

namespace OpenMS
{
	
	LinearMapping::LinearMapping()
		: slope_(1.0),
		  intercept_(0.0)
	{
	}
		
  LinearMapping::~LinearMapping()
  {
  }
    
	LinearMapping::LinearMapping(const LinearMapping& source)
		: slope_(source.slope_),
			intercept_(source.intercept_)
	{
	}

  LinearMapping& LinearMapping::operator = (const LinearMapping& source)
  {
    if (this==&source) return *this;
    
		slope_ = source.slope_;
		intercept_ = source.intercept_;
		
    return *this;
  }   
    
  void LinearMapping::apply(DPosition<1>& pos) const
  {
  	pos[0] = intercept_ + slope_ * pos[0];
  }
      
  void LinearMapping::apply( DoubleReal & pos) const
  {
    pos = intercept_ + slope_ * pos;
  }


  bool LinearMapping::operator==(const LinearMapping& rhs) const
	{
		return slope_==rhs.slope_ && intercept_==rhs.intercept_;
	}
    
	bool LinearMapping::operator!=(const LinearMapping& rhs) const
	{
		return slope_!=rhs.slope_ || intercept_!=rhs.intercept_;
	}

	std::ostream& operator<< (std::ostream& os, LinearMapping const & rhs)
	{
		return os << "-- LinearMapping: slope: " << rhs.getSlope() << " intercept: " << rhs.getIntercept() << " -- ";
	}

} // end of namespace OpenMS

