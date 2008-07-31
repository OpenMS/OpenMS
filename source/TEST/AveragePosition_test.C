// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/MATH/STATISTICS/AveragePosition.h>
///////////////////////////

using namespace OpenMS;
using namespace std;
using OpenMS::Math::AveragePosition;

START_TEST(AveragePosition, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

AveragePosition<3>* ptr = 0;
CHECK(AveragePosition())
{
	ptr = new AveragePosition<3>();
	TEST_NOT_EQUAL(ptr, 0)
}
RESULT

CHECK(~AveragePosition())
{
	delete ptr;
}
RESULT

CHECK((AveragePosition(AveragePosition const &rhs)))
{
	DPosition<4> pos1;
	pos1[0] = 1.0;
	pos1[1] = 2.0;
	pos1[2] = 3.0;
	pos1[3] = 4.0;
	
	DPosition<4> pos2;
	pos2[0] = 13.0;
	pos2[1] = 10.0;
	pos2[2] = 7.0;
	pos2[3] = 4.0;
	
	AveragePosition<4> avg;
	avg.add(pos1,6);
	avg.add(pos2);

	AveragePosition<4> avg_cpy = avg;
	
	TEST_REAL_EQUAL(avg.getPosition()[0],avg_cpy.getPosition()[0]);
	TEST_REAL_EQUAL(avg.getPosition()[1],avg_cpy.getPosition()[1]);
	TEST_REAL_EQUAL(avg.getPosition()[2],avg_cpy.getPosition()[2]);
	TEST_REAL_EQUAL(avg.getPosition()[3],avg_cpy.getPosition()[3]);
	TEST_REAL_EQUAL(avg.getWeight(),avg_cpy.getWeight());
}
RESULT

CHECK((PositionType const& getPosition() const))
{
	DPosition<4> pos1;
	pos1[0] = 1.0;
	pos1[1] = 2.0;
	pos1[2] = 3.0;
	pos1[3] = 4.0;
	
	DPosition<4> pos2;
	pos2[0] = 13.0;
	pos2[1] = 10.0;
	pos2[2] = 7.0;
	pos2[3] = 4.0;
	
	AveragePosition<4> avg;
	avg.add(pos1,-1);
	avg.add(pos2);

	TEST_REAL_EQUAL(avg.getPosition()[0],0);
	TEST_REAL_EQUAL(avg.getPosition()[1],0);
	TEST_REAL_EQUAL(avg.getPosition()[2],0);
	TEST_REAL_EQUAL(avg.getPosition()[3],0);
	TEST_REAL_EQUAL(avg.getWeight(),0);

	avg.add(pos1,4);
	avg.add(pos2);

	TEST_REAL_EQUAL(avg.getPosition()[0],5.8);
	TEST_REAL_EQUAL(avg.getPosition()[1],5.2);
	TEST_REAL_EQUAL(avg.getPosition()[2],4.6);
	TEST_REAL_EQUAL(avg.getPosition()[3],4);
	TEST_REAL_EQUAL(avg.getWeight(),5);
}
RESULT

CHECK((CoordinateType const& getWeight() const))
{
	AveragePosition<1> avg;
	avg.add(DPosition<1>(9),2);
	TEST_REAL_EQUAL(avg.getWeight(),2);
	TEST_REAL_EQUAL(avg.getPosition()[0],9);
	avg.add(DPosition<1>(9),3);
	TEST_REAL_EQUAL(avg.getWeight(),5);
	TEST_REAL_EQUAL(avg.getPosition()[0],9);
	avg.add(DPosition<1>(6),10);
	TEST_REAL_EQUAL(avg.getWeight(),15);
	TEST_REAL_EQUAL(avg.getPosition()[0],7);
}
RESULT

CHECK((void clear()))
{
	DPosition<4> pos1;
	pos1[0] = 1.0;
	pos1[1] = 2.0;
	pos1[2] = 3.0;
	pos1[3] = 4.0;
	AveragePosition<4> avg;
	avg.add(pos1,2);
	TEST_EQUAL(avg.getPosition(),pos1);
	TEST_REAL_EQUAL(avg.getWeight(),2);
	avg.clear();
	TEST_EQUAL(avg.getPosition(),DPosition<4>::zero);
	TEST_EQUAL(avg.getWeight(),0);
}
RESULT

CHECK((void add(PositionType position, CoordinateType const weight=1)))
{
	// already tested above
	NOT_TESTABLE;
}
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



