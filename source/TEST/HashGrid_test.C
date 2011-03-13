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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse, Holger Plattfaut $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/HashGrid.h>
#include <OpenMS/DATASTRUCTURES/GridElement.h>
#include <OpenMS/DATASTRUCTURES/DataPoint.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(HashGrid, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

HashGrid* ptr = 0;
START_SECTION(HashGrid())
{
	ptr = new HashGrid();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~HashGrid())
{
	delete ptr;
}
END_SECTION

START_SECTION((HashGrid(DoubleReal rt_threshold_, DoubleReal mz_threshold_)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void removeElement(GridElement *const element, Int x, Int y)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void removeElement(GridElement *const element)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void removeCell(GridCells::iterator loc)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void insert(GridElement *const element)))
{
  HashGrid tmp;
	DataPoint tmp2;
	tmp.insert(&tmp2);
	TEST_EQUAL(tmp.size(), 1);
}
END_SECTION

START_SECTION((void consoleOut() const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((Size size() const ))
{
	HashGrid tmp;
	DataPoint tmp2;
	tmp.insert(&tmp2);
  TEST_EQUAL(tmp.size(), 1);
}
END_SECTION

START_SECTION((DoubleReal getRTThreshold() const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((DoubleReal getMZThreshold() const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((Int getGridSizeX() const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((Int getGridSizeY() const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((Size getNumberOfElements() const ))
{
  HashGrid tmp;
	DataPoint tmp2;
	tmp.insert(&tmp2);
  TEST_EQUAL(tmp.size(), 1);
}
END_SECTION

START_SECTION((GridCells::iterator begin()))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((GridCells::iterator end()))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((GridCells::iterator find(std::pair< Int, Int > loc)))
{
  NOT_TESTABLE
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

