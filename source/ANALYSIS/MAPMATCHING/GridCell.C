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

#include <OpenMS/ANALYSIS/MAPMATCHING/GridCell.h>

using namespace std;

namespace OpenMS
{
  GridCell& GridCell::operator = (const GridCell& rhs)
  {
    if (&rhs==this) return *this;
      
    DRange<2>::operator = (rhs);
    mappings_.clear();
    for (UInt i=0; i<rhs.mappings_.size();++i)
    {
    	mappings_.push_back(new LinearMapping(*(dynamic_cast<LinearMapping*>(rhs.mappings_[i]))));
    }
    return *this;
  }

 std::ostream& operator << (std::ostream& os, const GridCell& grid)
 {
    os << "---------- GridCell BEGIN -----------------" << endl;
    os << "min: " << grid.min() << endl;
    os << "max: " << grid.max() << endl;
    for (UInt i=0; i<grid.getMappings().size();++i)
    {
    	const LinearMapping* m = dynamic_cast<const LinearMapping*>(grid.getMappings()[i]);
    	os << "mapping " << i << ": " << m->getSlope() << " " << m->getIntercept() << endl;
    } 
    os << "---------- GridCell END -------------------" << endl;;
    return os;
  }

}

