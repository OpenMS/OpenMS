// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/KERNEL/RawDataPoint1D.h>

///////////////////////////

START_TEST(RawDataPoint1D<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

RawDataPoint1D* d10_ptr = 0;
CHECK((RawDataPoint1D()))
	d10_ptr = new RawDataPoint1D;
	TEST_NOT_EQUAL(d10_ptr, 0)
RESULT

CHECK((~RawDataPoint1D()))
	delete d10_ptr;
RESULT

CHECK((IntensityType getIntensity() const))
	TEST_REAL_EQUAL(RawDataPoint1D().getIntensity(), 0.0)
RESULT

CHECK((PositionType const& getPosition() const))
	TEST_REAL_EQUAL(RawDataPoint1D().getPosition()[0], 0.0)
RESULT

CHECK((CoordinateType getMZ() const))
	TEST_REAL_EQUAL(RawDataPoint1D().getMZ(), 0.0)
RESULT

CHECK((CoordinateType getPos() const))
	TEST_REAL_EQUAL(RawDataPoint1D().getPos(), 0.0)
RESULT

CHECK((void setIntensity(IntensityType intensity)))
	RawDataPoint1D p;
 	p.setIntensity(17.8);
 	TEST_REAL_EQUAL(p.getIntensity(), 17.8)
RESULT

CHECK((void setPosition(PositionType const &position)))
	RawDataPoint1D::PositionType pos;
	pos[0] = 1.0;
	RawDataPoint1D p;
	p.setPosition(pos);
	TEST_REAL_EQUAL(p.getPosition()[0], 1.0)
RESULT

CHECK((PositionType& getPosition()))
	RawDataPoint1D::PositionType pos;
	pos[0] = 1.0;
	RawDataPoint1D p;
	p.getPosition() = pos;
	TEST_REAL_EQUAL(p.getPosition()[0], 1.0)
RESULT

CHECK((void setMZ(CoordinateTypemz)))
	RawDataPoint1D p;
	p.setMZ(5.0);
	TEST_REAL_EQUAL(p.getMZ(), 5.0)
RESULT

CHECK((void setPos(CoordinateTypepos)))
	RawDataPoint1D p;
	p.setPos(5.0);
	TEST_REAL_EQUAL(p.getPos(), 5.0)
RESULT

CHECK((RawDataPoint1D(const RawDataPoint1D& p)))
	RawDataPoint1D::PositionType pos;
	pos[0] = 21.21;
	RawDataPoint1D p;
	p.setIntensity(123.456);
	p.setPosition(pos);
	RawDataPoint1D::PositionType pos2;
	RawDataPoint1D::IntensityType i2;

	RawDataPoint1D copy_of_p(p);
	
	i2 = copy_of_p.getIntensity();
	pos2 = copy_of_p.getPosition();
	TEST_REAL_EQUAL(i2, 123.456)

	TEST_REAL_EQUAL(pos2[0], 21.21)
RESULT

CHECK((RawDataPoint1D& operator = (const RawDataPoint1D& rhs)))
	RawDataPoint1D::PositionType pos;
	pos[0] = 21.21;
	RawDataPoint1D p;
	p.setIntensity(123.456);
	p.setPosition(pos);
	RawDataPoint1D::PositionType pos2;
	RawDataPoint1D::IntensityType i2;

	RawDataPoint1D copy_of_p;
	copy_of_p = p;
	
	i2 = copy_of_p.getIntensity();
	pos2 = copy_of_p.getPosition();
	TEST_REAL_EQUAL(i2, 123.456)

	TEST_REAL_EQUAL(pos2[0], 21.21)
RESULT

CHECK((bool operator == (const RawDataPoint1D& rhs) const))
	RawDataPoint1D p1;
	RawDataPoint1D p2(p1);
	TEST_REAL_EQUAL(p1==p2, true)
	
	p1.setIntensity(5);
	TEST_REAL_EQUAL(p1==p2, false)
	p2.setIntensity(5);
	TEST_REAL_EQUAL(p1==p2, true)
	
	p1.getPosition()[0]=5;
	TEST_REAL_EQUAL(p1==p2, false)
	p2.getPosition()[0]=5;
	TEST_REAL_EQUAL(p1==p2, true)
RESULT

CHECK((bool operator != (const RawDataPoint1D& rhs) const))
	RawDataPoint1D p1;
	RawDataPoint1D p2(p1);
	TEST_REAL_EQUAL(p1!=p2, false)
	
	p1.setIntensity(5);
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.setIntensity(5);
	TEST_REAL_EQUAL(p1!=p2, false)
	
	p1.getPosition()[0]=5;
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.getPosition()[0]=5;
	TEST_REAL_EQUAL(p1!=p2, false)
RESULT

CHECK([EXTRA] class PositionLess)
	std::vector<RawDataPoint1D > v;
	RawDataPoint1D p;
	
	p.getPosition()[0]=3.0;
	v.push_back(p);

	p.getPosition()[0]=2.0;
	v.push_back(p);

	p.getPosition()[0]=1.0;
	v.push_back(p);
	
	std::sort(v.begin(), v.end(), RawDataPoint1D::PositionLess());
	TEST_REAL_EQUAL(v[0].getPosition()[0], 1.0)
	TEST_REAL_EQUAL(v[1].getPosition()[0], 2.0)
	TEST_REAL_EQUAL(v[2].getPosition()[0], 3.0)
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
