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
#include <OpenMS/DATASTRUCTURES/DataPoint.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(DataPoint, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DataPoint* ptr = 0;
START_SECTION(DataPoint())
{
	ptr = new DataPoint();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~DataPoint())
{
	delete ptr;
}
END_SECTION

START_SECTION((DataPoint(const DataPoint &copyin)))
{
  DataPoint tmp, tmp2;
	tmp.rt = 50;
	tmp.mz = 0.3;

	tmp2 = tmp;
	TEST_REAL_SIMILAR(tmp.rt, tmp2.rt);
	TEST_REAL_SIMILAR(tmp.mz, tmp2.mz);
}
END_SECTION

START_SECTION((DataPoint& operator=(const DataPoint &rhs)))
{
  DataPoint tmp;
	tmp.rt = 50;
	tmp.mz = 0.3;
 	
	//normal assignment
	DataPoint tmp2;
	tmp2 = tmp;
	TEST_REAL_SIMILAR(tmp2.rt, 50);
	TEST_REAL_SIMILAR(tmp2.mz, 0.3);

	//assignment of empty object
	tmp2 = DataPoint();
	TEST_REAL_SIMILAR(tmp2.rt, 0.0);
	TEST_REAL_SIMILAR(tmp2.mz, 0.0);
}
END_SECTION

START_SECTION((bool operator==(const DataPoint &rhs) const ))
{
  DataPoint tmp, tmp2;
	TEST_EQUAL(tmp == tmp2, true);

	tmp2.rt = 50;
	TEST_EQUAL(tmp == tmp2, false);

	tmp2 = tmp;
	tmp2.mz = 0.1;
	TEST_EQUAL(tmp == tmp2, false);
}
END_SECTION

START_SECTION((bool operator!=(const DataPoint &rhs) const ))
{
  DataPoint tmp, tmp2;
	TEST_EQUAL(tmp != tmp2, false);

	tmp2.rt = 50;
	TEST_EQUAL(tmp != tmp2, true);

	tmp2 = tmp;
	tmp2.mz = 0.1;
	TEST_EQUAL(tmp != tmp2, true);
}
END_SECTION

START_SECTION((bool operator<(const DataPoint &rhs) const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((Int getID() const ))
{
  DataPoint tmp;
  TEST_EQUAL(tmp.getID(), 0);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

