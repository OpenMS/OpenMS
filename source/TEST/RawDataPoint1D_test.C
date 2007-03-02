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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/KERNEL/RawDataPoint1D.h>

///////////////////////////

START_TEST(RawDataPoint1D<D>, "$Id: DRawDataPoint_test.C 1300 2007-01-18 07:27:04Z marc_sturm $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

RawDataPoint1D* d10_ptr = 0;
CHECK(RawDataPoint1D())
	d10_ptr = new RawDataPoint1D;
	TEST_NOT_EQUAL(d10_ptr, 0)
RESULT

CHECK(~RawDataPoint1D())
	delete d10_ptr;
RESULT

CHECK(const IntensityType& getIntensity() const)
	const RawDataPoint1D p;
	TEST_REAL_EQUAL(p.getIntensity(), 0.0)
RESULT

CHECK(const PositionType& getPos() const)
	const RawDataPoint1D	p;
	TEST_REAL_EQUAL(p.getPos()[0], 0.0)
RESULT

CHECK(IntensityType& getIntensity())
	RawDataPoint1D p;
	TEST_REAL_EQUAL(p.getIntensity(), 0.0)
	p.setIntensity(123.456);
	TEST_REAL_EQUAL(p.getIntensity(), 123.456)
	p.setIntensity(-0.12345);
	TEST_REAL_EQUAL(p.getIntensity(), -0.12345)
	p.setIntensity(0.0);
	TEST_REAL_EQUAL(p.getIntensity(), 0.0)
RESULT

CHECK(void setIntensity(const IntensityType& intensity))
	RawDataPoint1D p;
  p.setIntensity(17.8);
  TEST_REAL_EQUAL(p.getIntensity(), 17.8)
RESULT


CHECK(PositionType& getPos())
	RawDataPoint1D::PositionType pos;
	RawDataPoint1D p;
	pos = p.getPos();
	TEST_REAL_EQUAL(pos[0], 0.0)
	pos[0] = 1.0;
	p.setPos(pos);
	RawDataPoint1D::PositionType pos2(p.getPos());
	TEST_REAL_EQUAL(pos2[0], 1.0)
RESULT

CHECK(RawDataPoint1D(const RawDataPoint1D& p))
	RawDataPoint1D::PositionType pos;
	pos[0] = 21.21;
	RawDataPoint1D p;
	p.setIntensity(123.456);
	p.setPos(pos);
	RawDataPoint1D::PositionType pos2;
	RawDataPoint1D::IntensityType i2;

	RawDataPoint1D copy_of_p(p);
	
	i2 = copy_of_p.getIntensity();
	pos2 = copy_of_p.getPos();
	TEST_REAL_EQUAL(i2, 123.456)

	TEST_REAL_EQUAL(pos2[0], 21.21)
RESULT

CHECK(RawDataPoint1D& operator = (const RawDataPoint1D& rhs))
	RawDataPoint1D::PositionType pos;
	pos[0] = 21.21;
	RawDataPoint1D p;
	p.setIntensity(123.456);
	p.setPos(pos);
	RawDataPoint1D::PositionType pos2;
	RawDataPoint1D::IntensityType i2;

	RawDataPoint1D copy_of_p;
	copy_of_p = p;
	
	i2 = copy_of_p.getIntensity();
	pos2 = copy_of_p.getPos();
	TEST_REAL_EQUAL(i2, 123.456)

	TEST_REAL_EQUAL(pos2[0], 21.21)
RESULT

CHECK(bool operator == (const RawDataPoint1D& rhs) const)
	RawDataPoint1D p1;
	RawDataPoint1D p2(p1);
	TEST_REAL_EQUAL(p1==p2, true)
	
	p1.setIntensity(5);
	TEST_REAL_EQUAL(p1==p2, false)
	p2.setIntensity(5);
	TEST_REAL_EQUAL(p1==p2, true)
	
	p1.getPos()[0]=5;
	TEST_REAL_EQUAL(p1==p2, false)
	p2.getPos()[0]=5;
	TEST_REAL_EQUAL(p1==p2, true)
RESULT

CHECK(bool operator != (const RawDataPoint1D& rhs) const)
	RawDataPoint1D p1;
	RawDataPoint1D p2(p1);
	TEST_REAL_EQUAL(p1!=p2, false)
	
	p1.setIntensity(5);
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.setIntensity(5);
	TEST_REAL_EQUAL(p1!=p2, false)
	
	p1.getPos()[0]=5;
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.getPos()[0]=5;
	TEST_REAL_EQUAL(p1!=p2, false)
RESULT

CHECK([EXTRA] class PositionLess)
	std::vector<RawDataPoint1D > v;
	RawDataPoint1D p;
	
	p.getPos()[0]=3.0;
	v.push_back(p);

	p.getPos()[0]=2.0;
	v.push_back(p);

	p.getPos()[0]=1.0;
	v.push_back(p);
	
	std::sort(v.begin(), v.end(), RawDataPoint1D::PositionLess());
	TEST_REAL_EQUAL(v[0].getPos()[0], 1.0)
	TEST_REAL_EQUAL(v[1].getPos()[0], 2.0)
	TEST_REAL_EQUAL(v[2].getPos()[0], 3.0)
RESULT

CHECK([EXTRA] struct IntensityLess)
	std::vector<RawDataPoint1D > v;
	RawDataPoint1D p;
	
	p.setIntensity(2.5);
	v.push_back(p);

	p.setIntensity(3.5);
	v.push_back(p);

	p.setIntensity(1.5);
	v.push_back(p);
	
	std::sort(v.begin(), v.end(), RawDataPoint1D::IntensityLess());
	TEST_REAL_EQUAL(v[0].getIntensity(), 1.5)
	TEST_REAL_EQUAL(v[1].getIntensity(), 2.5)
	TEST_REAL_EQUAL(v[2].getIntensity(), 3.5)

	v[0]=v[2];
	v[2]=p;
	std::sort(v.begin(), v.end(), RawDataPoint1D::IntensityLess());
	TEST_REAL_EQUAL(v[0].getIntensity(), 1.5)
	TEST_REAL_EQUAL(v[1].getIntensity(), 2.5)
	TEST_REAL_EQUAL(v[2].getIntensity(), 3.5)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
