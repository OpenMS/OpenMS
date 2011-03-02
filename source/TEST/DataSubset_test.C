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
#include <OpenMS/DATASTRUCTURES/DataSubset.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(DataSubset, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DataSubset* ptr = 0;
START_SECTION(DataSubset())
{
	ptr = new DataSubset();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~DataSubset())
{
	delete ptr;
}
END_SECTION

START_SECTION((DataSubset(DataPoint &data_point)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((DataSubset(const DataSubset &copy)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((DataSubset(const DataSubset *copy_ptr)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((Int operator<(const DataSubset &el) const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((Int size()))
{
  DataSubset tmp;
	TEST_EQUAL(tmp.size(), 0);
}
END_SECTION

START_SECTION((Int getID() const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((bool operator!=(const DataSubset &el) const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((bool operator==(const DataSubset &el) const ))
{
  NOT_TESTABLE
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

