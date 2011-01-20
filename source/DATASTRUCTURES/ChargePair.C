// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/ChargePair.h>

namespace OpenMS
{
	

		
  std::ostream& operator << (std::ostream& os, const ChargePair& cp)
  {
    os  << "---------- ChargePair -----------------\n"
		    << "Mass Diff: " << cp.getMassDiff()<< std::endl
		    << "Compomer: " << cp.getCompomer() << std::endl
		    << "Charge: " << cp.getCharge(0) << " : " << cp.getCharge(1) << std::endl
		    << "Element Index: " << cp.getElementIndex(0) << " : " << cp.getElementIndex(1) << std::endl;
    return os;
  }
}

