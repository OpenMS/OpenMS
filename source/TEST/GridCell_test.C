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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/ANALYSIS/MAPMATCHING/GridCell.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/BaseMapping.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/LinearMapping.h>

#include <vector>

///////////////////////////

START_TEST(GridCell, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

typedef std::vector<BaseMapping*> MappingVector;

GridCell* d10_ptr = 0;
CHECK((GridCell()))
	d10_ptr = new GridCell;
	TEST_NOT_EQUAL(d10_ptr, 0)
RESULT

CHECK((~GridCell()))
	delete d10_ptr;
RESULT

CHECK((GridCell(const GridCell& gc)))

	GridCell c1(0.0, 0.0, 10.0, 10.0);
	BaseMapping* m1 = new LinearMapping(3.0,4.0);
	MappingVector mvec1;
	mvec1.push_back(m1);	
	c1.setMappings(mvec1);
	
	GridCell c2(c1);

	TEST_EQUAL(c1.min(),c2.min());
	TEST_EQUAL(c1.max(),c2.max());
	TEST_EQUAL(c1.getMappings()==c2.getMappings(),true);

RESULT

CHECK((bool operator == (const GridCell& rhs) const))

	GridCell c1(0.0, 0.0, 10.0, 10.0);
	GridCell c2(0.0, 0.0, 10.0, 10.0);

	TEST_EQUAL(c1==c2,true);
	
	c2.setMaxX(4);
	TEST_EQUAL(c1!=c2,true);

RESULT

CHECK((MappingVector& getMappings()))

	GridCell c1;
	BaseMapping* m1 = new LinearMapping(3.0,4.0);
	MappingVector mvec1;
	mvec1.push_back(m1);
	
	c1.setMappings(mvec1);

	MappingVector mvec2 = c1.getMappings();

	TEST_EQUAL(mvec1==mvec2,true);
	
	BaseMapping* m2 = new LinearMapping(5.0,60.0);
	MappingVector mvec3;
	mvec3.push_back(m2);
	c1.setMappings(mvec3);
	
	TEST_EQUAL(c1.getMappings()==mvec2,false);

RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
