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
#include <OpenMS/ANALYSIS/MAPMATCHING/GridCell.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/Grid.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/BaseMapping.h>

#include <string>

///////////////////////////

using namespace std;
using namespace OpenMS;

///////////////////////////

/////////////////////////////////////////////////////////////

START_TEST(Grid, "$Id: Grid_test.C 1586 2007-03-01 17:59:10Z elange $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


Grid* pl_ptr = 0;
CHECK((Grid()))
	pl_ptr = new Grid();
	TEST_NOT_EQUAL(pl_ptr, 0)
RESULT

CHECK((~Grid()))
	delete pl_ptr;
RESULT

CHECK((Grid(const Grid& grid)))
			
	GridCell c1(0.0,0.0,2.0,2.0);
	GridCell c2(3.0,3.0,6.0,6.0);
	
	Grid grid;
	grid.push_back(c1);
	grid.push_back(c2);
	
	Grid grid_copy(grid);
	
	TEST_EQUAL(grid_copy.size(),2);
	
	Grid::const_iterator cit = grid_copy.begin();
	TEST_EQUAL(cit->minX(),0.0);
	TEST_EQUAL(cit->minY(),0.0);		
	
	cit++;
	TEST_EQUAL(cit->maxX(),6.0);
	TEST_EQUAL(cit->maxY(),6.0);
	
RESULT

CHECK((Grid& operator = (const Grid& rhs)))

	GridCell c1(0.0,0.0,2.0,2.0);
	GridCell c2(3.0,3.0,6.0,6.0);
	
	Grid grid;
	grid.push_back(c1);
	grid.push_back(c2);
	
	Grid grid_copy = grid;
	
	TEST_EQUAL(grid_copy.size(),2);
	
	Grid::const_iterator cit = grid_copy.begin();
	TEST_EQUAL(cit->minX(),0.0);
	TEST_EQUAL(cit->minY(),0.0);		
	
	cit++;
	TEST_EQUAL(cit->maxX(),6.0);
	TEST_EQUAL(cit->maxY(),6.0);
	
RESULT

CHECK((bool operator == (const Grid& rhs) const))

	GridCell c1(0.0,0.0,2.0,2.0);
	GridCell c2(3.0,3.0,6.0,6.0);
	
	Grid grid;
	grid.push_back(c1);
	grid.push_back(c2);
	
	Grid grid_copy;
	grid_copy.push_back(c1);
	grid_copy.push_back(c2);
	
	TEST_EQUAL(grid_copy==grid,true);
	
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
