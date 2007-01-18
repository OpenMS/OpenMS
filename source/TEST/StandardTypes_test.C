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
// $Maintainer: Clemens Groepl  $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/KERNEL/StandardTypes.h>

///////////////////////////

START_TEST(StandardTypes, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

CHECK(StandardTypes)

	STATUS("Note: Here we only check whether the typedefs (as such) are fine.")

RESULT

// super duper macro

#define CHECK_A_TYPEDEF(Type)												\
{																										\
	Type* ptr = 0;																		\
  CHECK(Type())																			\
    ptr = new Type;																	\
    TEST_NOT_EQUAL(ptr, 0)													\
  RESULT																						\
  CHECK(~Type())																		\
    delete ptr;																			\
  RESULT																						\
}


CHECK_A_TYPEDEF(RawDataPoint)
CHECK_A_TYPEDEF(RawDataPoint2D)
CHECK_A_TYPEDEF(RawSpectrum)
CHECK_A_TYPEDEF(RawMap)
CHECK_A_TYPEDEF(Peak)
CHECK_A_TYPEDEF(Peak2D)
CHECK_A_TYPEDEF(PeakSpectrum)
CHECK_A_TYPEDEF(PeakMap)
CHECK_A_TYPEDEF(Feature)
CHECK_A_TYPEDEF(FeatureMap)



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


