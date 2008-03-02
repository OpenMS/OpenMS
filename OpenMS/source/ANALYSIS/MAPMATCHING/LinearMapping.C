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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/LinearMapping.h>

namespace OpenMS
{
  // assignment operator
  LinearMapping& LinearMapping::operator = (const LinearMapping& source)
  {
    if (this==&source) return *this;
      
    BaseMapping::operator = (source);   
    param_.setValue("slope",source.slope_);
    param_.setValue("intercept",source.intercept_);
    updateMembers_();
    return *this;
  }   
    
  void LinearMapping::setParam(DoubleReal sl, DoubleReal in)
  {
    param_.setValue("slope",sl);
    param_.setValue("intercept",in);
    updateMembers_();
  }
      
  void LinearMapping::apply(DPosition<1>& pos) const
  {
   pos[0] = intercept_ + slope_ * pos[0];
  }
      
  void LinearMapping::apply( DoubleReal & pos) const
  {
    pos *= slope_;
    pos += intercept_;
  }
		
} // end of namespace OpenMS

