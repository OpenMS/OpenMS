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
#include <OpenMS/DATASTRUCTURES/QTSILACCluster.h>
#include <OpenMS/DATASTRUCTURES/DataPoint.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(QTSILACCluster, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

QTSILACCluster* ptr = 0;
START_SECTION(QTSILACCluster())
{
	ptr = new QTSILACCluster();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~QTSILACCluster())
{
	delete ptr;
}
END_SECTION

START_SECTION((QTSILACCluster(DataPoint *center_point_)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((DoubleReal getCenterRT()))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((DoubleReal getCenterMZ()))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((Size size() const ))
{
  QTSILACCluster tmp;
	TEST_EQUAL(tmp.size(), 0);
}
END_SECTION

START_SECTION((bool operator<(const QTSILACCluster &cluster) const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void add(DataPoint *element)))
{
  QTSILACCluster tmp;
	DataPoint* tmp2;
	tmp.add(tmp2);
	TEST_EQUAL(tmp.size(), 1);
}
END_SECTION

START_SECTION((std::set<DataPoint*> getClusterMembers()))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((std::pair<DoubleReal,DoubleReal> getDiameters(DataPoint *point)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((bool contains(DataPoint *element)))
{
  NOT_TESTABLE
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

