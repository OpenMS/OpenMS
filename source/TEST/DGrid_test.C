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
#include <OpenMS/ANALYSIS/MAPMATCHING/DGridCell.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DGrid.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DBaseMapping.h>

#include <string>

///////////////////////////

using namespace std;
using namespace OpenMS;

///////////////////////////

/////////////////////////////////////////////////////////////

START_TEST(DGrid<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


DGrid<2>* pl_ptr = 0;
CHECK((DGrid()))
	pl_ptr = new DGrid<2>();
	TEST_NOT_EQUAL(pl_ptr, 0)
RESULT

CHECK((~DGrid()))
	delete pl_ptr;
RESULT

CHECK((DGrid(const DGrid& grid)))
			
	DGridCell<2> c1(0.0,0.0,2.0,2.0);
	DGridCell<2> c2(3.0,3.0,6.0,6.0);
	
	DGrid<2> grid;
	grid.push_back(c1);
	grid.push_back(c2);
	
	DGrid<2> grid_copy(grid);
	
	TEST_EQUAL(grid_copy.size(),2);
	
	DGrid<2>::const_iterator cit = grid_copy.begin();
	TEST_EQUAL(cit->minX(),0.0);
	TEST_EQUAL(cit->minY(),0.0);		
	
	cit++;
	TEST_EQUAL(cit->maxX(),6.0);
	TEST_EQUAL(cit->maxY(),6.0);
	
RESULT

CHECK((DGrid& operator = (const DGrid& rhs)))

	DGridCell<2> c1(0.0,0.0,2.0,2.0);
	DGridCell<2> c2(3.0,3.0,6.0,6.0);
	
	DGrid<2> grid;
	grid.push_back(c1);
	grid.push_back(c2);
	
	DGrid<2> grid_copy = grid;
	
	TEST_EQUAL(grid_copy.size(),2);
	
	DGrid<2>::const_iterator cit = grid_copy.begin();
	TEST_EQUAL(cit->minX(),0.0);
	TEST_EQUAL(cit->minY(),0.0);		
	
	cit++;
	TEST_EQUAL(cit->maxX(),6.0);
	TEST_EQUAL(cit->maxY(),6.0);
	
RESULT

CHECK((bool operator == (const DGrid& rhs) const))

	DGridCell<2> c1(0.0,0.0,2.0,2.0);
	DGridCell<2> c2(3.0,3.0,6.0,6.0);
	
	DGrid<2> grid;
	grid.push_back(c1);
	grid.push_back(c2);
	
	DGrid<2> grid_copy;
	grid_copy.push_back(c1);
	grid_copy.push_back(c2);
	
	TEST_EQUAL(grid_copy==grid,true);
	
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
