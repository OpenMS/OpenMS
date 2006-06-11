// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: DLinearMapping_test.C,v 1.1 2006/04/10 16:11:13 ole_st Exp $
// $Author: ole_st $
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/ANALYSIS/MAPMATCHING/DLinearMapping.h>

#include <vector>

///////////////////////////

START_TEST(DLinearMapping<D>, "$Id: DLinearMapping_test.C,v 1.1 2006/04/10 16:11:13 ole_st Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;


DLinearMapping<10>* d10_ptr = 0;
CHECK(DLinearMapping())
	d10_ptr = new DLinearMapping<10>;
	TEST_NOT_EQUAL(d10_ptr, 0)
RESULT

CHECK(~DLinearMapping())
	delete d10_ptr;
RESULT

CHECK(DLinearMapping(const DLinearMapping& gc))

	DLinearMapping<2> c1(0.0, 10.0);
	DLinearMapping<2> c2(c1);
	
	TEST_EQUAL(c1.getParam()==c2.getParam(),true)	

RESULT

CHECK(operator=)
	DLinearMapping<2> c1(0.0, 10.0);
	DLinearMapping<2> c2 = c1;
	
	TEST_EQUAL(c1.getParam()==c2.getParam(),true)	

RESULT

CHECK(operator==)
	DLinearMapping<2> c1(0.0, 10.0);
	DLinearMapping<2> c2(0.0, 10.0);
	
	TEST_EQUAL(c1==c2,true)	

RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
