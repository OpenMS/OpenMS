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

#include <OpenMS/KERNEL/DRawDataPoint.h>

///////////////////////////

START_TEST(DRawDataPoint<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

DRawDataPoint<10>* d10_ptr = 0;
CHECK(DRawDataPoint())
	d10_ptr = new DRawDataPoint<10>;
	TEST_NOT_EQUAL(d10_ptr, 0)
RESULT

CHECK(~DRawDataPoint())
	delete d10_ptr;
RESULT

CHECK(const IntensityType& getIntensity() const)
	const DRawDataPoint<10> p;
	TEST_REAL_EQUAL(p.getIntensity(), 0.0)
RESULT

CHECK(const PositionType& getPosition() const)
	const DRawDataPoint<10>	p;
	TEST_REAL_EQUAL(p.getPosition()[0], 0.0)
	TEST_REAL_EQUAL(p.getPosition()[1], 0.0)
	TEST_REAL_EQUAL(p.getPosition()[2], 0.0)
	TEST_REAL_EQUAL(p.getPosition()[3], 0.0)
	TEST_REAL_EQUAL(p.getPosition()[4], 0.0)
	TEST_REAL_EQUAL(p.getPosition()[5], 0.0)
	TEST_REAL_EQUAL(p.getPosition()[6], 0.0)
	TEST_REAL_EQUAL(p.getPosition()[7], 0.0)
	TEST_REAL_EQUAL(p.getPosition()[8], 0.0)
	TEST_REAL_EQUAL(p.getPosition()[9], 0.0)
RESULT

CHECK(void setPos(CoordinateType const& coordinate))
	DRawDataPoint<1> p0;
  p0.setPos(12345.0);
	const DRawDataPoint<1> p1(p0);
	TEST_REAL_EQUAL(p1.getPos(), 12345.0)
RESULT

CHECK(void setPos(Size const i, const CoordinateType& coordinate))
	DRawDataPoint<4> p0;
  p0.setPos(2,12345.0);
	const DRawDataPoint<4> p1(p0);
	TEST_REAL_EQUAL(p1.getPos(2), 12345.0)
RESULT

CHECK(void setPosition(PositionType const& position))
	DRawDataPoint<2> p0;
  p0.setPos(0,876);
  p0.setPos(1,12345.0);
	DRawDataPoint<2> p1;
  p1.setPosition(p0.getPosition());
	TEST_REAL_EQUAL(p1.getPosition()[0], 876)
	TEST_REAL_EQUAL(p1.getPosition()[1], 12345.0)
RESULT

CHECK(CoordinateType& getPos())
	DRawDataPoint<1> p0;
  p0.setPos(12345.0);
	const DRawDataPoint<1> p1(p0);
	TEST_REAL_EQUAL(p1.getPos(), 12345.0)
RESULT

CHECK(CoordinateType const& getPos() const)
	DRawDataPoint<1> p0;
  p0.setPos(12345.0);
	const DRawDataPoint<1> p1(p0);
	TEST_REAL_EQUAL(p1.getPos(), 12345.0)
RESULT

CHECK(CoordinateType& getPos(Size const i))
	DRawDataPoint<2> p;
	TEST_REAL_EQUAL(p.getPos(1), 0.0)
RESULT

CHECK(CoordinateType const& getPos(Size const i) const)
	const DRawDataPoint<2>	p;
	TEST_REAL_EQUAL(p.getPos(1), 0.0)
RESULT

CHECK(IntensityType& getIntensity())
	DRawDataPoint<3> p;
	TEST_REAL_EQUAL(p.getIntensity(), 0.0)
	p.getIntensity() = 123.456;
	TEST_REAL_EQUAL(p.getIntensity(), 123.456)
	p.getIntensity() = -0.12345;
	TEST_REAL_EQUAL(p.getIntensity(), -0.12345)
	p.getIntensity() = 0.0;
	TEST_REAL_EQUAL(p.getIntensity(), 0.0)
RESULT

CHECK(void setIntensity(const IntensityType& intensity))
	DRawDataPoint<4> p;
  p.setIntensity(17.8);
  TEST_REAL_EQUAL(p.getIntensity(), 17.8)
RESULT


CHECK(PositionType& getPosition())
	DRawDataPoint<3>::PositionType pos;
	DRawDataPoint<3> p;
	pos = p.getPosition();
	TEST_REAL_EQUAL(pos[0], 0.0)
	TEST_REAL_EQUAL(pos[1], 0.0)
	TEST_REAL_EQUAL(pos[2], 0.0)
	pos[0] = 1.0;
	pos[1] = 2.0;
	pos[2] = 3.0;
	p.getPosition() = pos;
	DRawDataPoint<3>::PositionType pos2(p.getPosition());
	TEST_REAL_EQUAL(pos2[0], 1.0)
	TEST_REAL_EQUAL(pos2[1], 2.0)
	TEST_REAL_EQUAL(pos2[2], 3.0)	
RESULT

CHECK(DRawDataPoint(const DRawDataPoint& p))
	DRawDataPoint<3>::PositionType pos;
	pos[0] = 21.21;
	pos[1] = 22.22;
	pos[2] = 23.23;
	DRawDataPoint<3> p;
	p.getIntensity() = 123.456;
	p.getPosition() = pos;
	DRawDataPoint<3>::PositionType pos2;
	DRawDataPoint<3>::IntensityType i2;

	DRawDataPoint<3> copy_of_p(p);
	
	i2 = copy_of_p.getIntensity();
	pos2 = copy_of_p.getPosition();
	TEST_REAL_EQUAL(i2, 123.456)

	TEST_REAL_EQUAL(pos2[0], 21.21)
	TEST_REAL_EQUAL(pos2[1], 22.22)
	TEST_REAL_EQUAL(pos2[2], 23.23)	
RESULT

CHECK(DRawDataPoint& operator = (const DRawDataPoint& rhs))
	DRawDataPoint<3>::PositionType pos;
	pos[0] = 21.21;
	pos[1] = 22.22;
	pos[2] = 23.23;
	DRawDataPoint<3> p;
	p.getIntensity() = 123.456;
	p.getPosition() = pos;
	DRawDataPoint<3>::PositionType pos2;
	DRawDataPoint<3>::IntensityType i2;

	DRawDataPoint<3> copy_of_p;
	copy_of_p = p;
	
	i2 = copy_of_p.getIntensity();
	pos2 = copy_of_p.getPosition();
	TEST_REAL_EQUAL(i2, 123.456)

	TEST_REAL_EQUAL(pos2[0], 21.21)
	TEST_REAL_EQUAL(pos2[1], 22.22)
	TEST_REAL_EQUAL(pos2[2], 23.23)	
RESULT

CHECK(bool operator == (const DRawDataPoint& rhs) const)
	DRawDataPoint<1> p1;
	DRawDataPoint<1> p2(p1);
	TEST_REAL_EQUAL(p1==p2, true)
	
	p1.getIntensity()=5;
	TEST_REAL_EQUAL(p1==p2, false)
	p2.getIntensity()=5;
	TEST_REAL_EQUAL(p1==p2, true)
	
	p1.getPosition()[0]=5;
	TEST_REAL_EQUAL(p1==p2, false)
	p2.getPosition()[0]=5;
	TEST_REAL_EQUAL(p1==p2, true)
RESULT

CHECK(bool operator != (const DRawDataPoint& rhs) const)
	DRawDataPoint<1> p1;
	DRawDataPoint<1> p2(p1);
	TEST_REAL_EQUAL(p1!=p2, false)
	
	p1.getIntensity()=5;
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.getIntensity()=5;
	TEST_REAL_EQUAL(p1!=p2, false)
	
	p1.getPosition()[0]=5;
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.getPosition()[0]=5;
	TEST_REAL_EQUAL(p1!=p2, false)
RESULT

CHECK([EXTRA] class PositionLess)
	std::vector<DRawDataPoint<2> > v;
	DRawDataPoint<2> p;
	
	p.getPosition()[0]=3.0;
	p.getPosition()[1]=2.5;
	v.push_back(p);

	p.getPosition()[0]=2.0;
	p.getPosition()[1]=3.5;
	v.push_back(p);

	p.getPosition()[0]=1.0;
	p.getPosition()[1]=1.5;
	v.push_back(p);
	
	std::sort(v.begin(), v.end(), DRawDataPoint<2>::PositionLess());
	TEST_REAL_EQUAL(v[0].getPosition()[0], 1.0)
	TEST_REAL_EQUAL(v[1].getPosition()[0], 2.0)
	TEST_REAL_EQUAL(v[2].getPosition()[0], 3.0)
	TEST_REAL_EQUAL(v[0].getPosition()[1], 1.5)
	TEST_REAL_EQUAL(v[1].getPosition()[1], 3.5)
	TEST_REAL_EQUAL(v[2].getPosition()[1], 2.5)

	std::sort(v.begin(), v.end(), DRawDataPoint<2>::NthPositionLess<1>());
	TEST_REAL_EQUAL(v[0].getPosition()[1], 1.5)
	TEST_REAL_EQUAL(v[1].getPosition()[1], 2.5)
	TEST_REAL_EQUAL(v[2].getPosition()[1], 3.5)
	TEST_REAL_EQUAL(v[0].getPosition()[0], 1.0)
	TEST_REAL_EQUAL(v[1].getPosition()[0], 3.0)
	TEST_REAL_EQUAL(v[2].getPosition()[0], 2.0)
RESULT

CHECK([EXTRA] struct NthPositionLess)

	std::vector<DRawDataPoint<2> > v;
	DRawDataPoint<2> p;
	
	p.getPosition()[0]=3.0;
	p.getPosition()[1]=2.5;
	v.push_back(p);

	p.getPosition()[0]=2.0;
	p.getPosition()[1]=3.5;
	v.push_back(p);

	p.getPosition()[0]=1.0;
	p.getPosition()[1]=1.5;
	v.push_back(p);

	std::sort(v.begin(), v.end(), DRawDataPoint<2>::NthPositionLess<1>());
	TEST_REAL_EQUAL(v[0].getPosition()[1], 1.5)
	TEST_REAL_EQUAL(v[1].getPosition()[1], 2.5)
	TEST_REAL_EQUAL(v[2].getPosition()[1], 3.5)
	TEST_REAL_EQUAL(v[0].getPosition()[0], 1.0)
	TEST_REAL_EQUAL(v[1].getPosition()[0], 3.0)
	TEST_REAL_EQUAL(v[2].getPosition()[0], 2.0)
RESULT

CHECK([EXTRA] struct NthPositionLess<i>::DIMENSION)
	TEST_EQUAL(DRawDataPoint<7>::NthPositionLess<5>::DIMENSION, 5)
RESULT

CHECK([EXTRA] struct IntensityLess)
	std::vector<DRawDataPoint<2> > v;
	DRawDataPoint<2> p;
	
	p.getIntensity()=2.5;
	v.push_back(p);

	p.getIntensity()=3.5;
	v.push_back(p);

	p.getIntensity()=1.5;
	v.push_back(p);
	
	std::sort(v.begin(), v.end(), DRawDataPoint<2>::IntensityLess());
	TEST_REAL_EQUAL(v[0].getIntensity(), 1.5)
	TEST_REAL_EQUAL(v[1].getIntensity(), 2.5)
	TEST_REAL_EQUAL(v[2].getIntensity(), 3.5)

	v[0]=v[2];
	v[2]=p;
	std::sort(v.begin(), v.end(), DRawDataPoint<2>::IntensityLess());
	TEST_REAL_EQUAL(v[0].getIntensity(), 1.5)
	TEST_REAL_EQUAL(v[1].getIntensity(), 2.5)
	TEST_REAL_EQUAL(v[2].getIntensity(), 3.5)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
