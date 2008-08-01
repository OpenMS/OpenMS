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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/KERNEL/Peak1D.h>

///////////////////////////

START_TEST(Peak1D<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

Peak1D* d10_ptr = 0;
CHECK((Peak1D()))
	d10_ptr = new Peak1D;
	TEST_NOT_EQUAL(d10_ptr, 0)
RESULT

CHECK((~Peak1D()))
	delete d10_ptr;
RESULT

CHECK((IntensityType getIntensity() const))
	TEST_REAL_EQUAL(Peak1D().getIntensity(), 0.0)
RESULT

CHECK((PositionType const& getPosition() const))
	TEST_REAL_EQUAL(Peak1D().getPosition()[0], 0.0)
RESULT

CHECK((CoordinateType getMZ() const))
	TEST_REAL_EQUAL(Peak1D().getMZ(), 0.0)
RESULT

CHECK((CoordinateType getPos() const))
	TEST_REAL_EQUAL(Peak1D().getPos(), 0.0)
RESULT

CHECK((void setIntensity(IntensityType intensity)))
	Peak1D p;
 	p.setIntensity(17.8);
 	TEST_REAL_EQUAL(p.getIntensity(), 17.8)
RESULT

CHECK((void setPosition(PositionType const &position)))
	Peak1D::PositionType pos;
	pos[0] = 1.0;
	Peak1D p;
	p.setPosition(pos);
	TEST_REAL_EQUAL(p.getPosition()[0], 1.0)
RESULT

CHECK((PositionType& getPosition()))
	Peak1D::PositionType pos;
	pos[0] = 1.0;
	Peak1D p;
	p.getPosition() = pos;
	TEST_REAL_EQUAL(p.getPosition()[0], 1.0)
RESULT

CHECK((void setMZ(CoordinateTypemz)))
	Peak1D p;
	p.setMZ(5.0);
	TEST_REAL_EQUAL(p.getMZ(), 5.0)
RESULT

CHECK((void setPos(CoordinateTypepos)))
	Peak1D p;
	p.setPos(5.0);
	TEST_REAL_EQUAL(p.getPos(), 5.0)
RESULT

CHECK((Peak1D(const Peak1D& p)))
	Peak1D::PositionType pos;
	pos[0] = 21.21;
	Peak1D p;
	p.setIntensity(123.456);
	p.setPosition(pos);
	Peak1D::PositionType pos2;
	Peak1D::IntensityType i2;

	Peak1D copy_of_p(p);
	
	i2 = copy_of_p.getIntensity();
	pos2 = copy_of_p.getPosition();
	TEST_REAL_EQUAL(i2, 123.456)

	TEST_REAL_EQUAL(pos2[0], 21.21)
RESULT

CHECK((Peak1D& operator = (const Peak1D& rhs)))
	Peak1D::PositionType pos;
	pos[0] = 21.21;
	Peak1D p;
	p.setIntensity(123.456);
	p.setPosition(pos);
	Peak1D::PositionType pos2;
	Peak1D::IntensityType i2;

	Peak1D copy_of_p;
	copy_of_p = p;
	
	i2 = copy_of_p.getIntensity();
	pos2 = copy_of_p.getPosition();
	TEST_REAL_EQUAL(i2, 123.456)

	TEST_REAL_EQUAL(pos2[0], 21.21)
RESULT

CHECK((bool operator == (const Peak1D& rhs) const))
	Peak1D p1;
	Peak1D p2(p1);
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

CHECK((bool operator != (const Peak1D& rhs) const))
	Peak1D p1;
	Peak1D p2(p1);
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
	std::vector<Peak1D > v;
	Peak1D p;
	
	p.getPosition()[0]=3.0;
	v.push_back(p);

	p.getPosition()[0]=2.0;
	v.push_back(p);

	p.getPosition()[0]=1.0;
	v.push_back(p);
	
	std::sort(v.begin(), v.end(), Peak1D::PositionLess());
	TEST_REAL_EQUAL(v[0].getPosition()[0], 1.0)
	TEST_REAL_EQUAL(v[1].getPosition()[0], 2.0)
	TEST_REAL_EQUAL(v[2].getPosition()[0], 3.0)
RESULT

CHECK([EXTRA] struct IntensityLess)
	std::vector<Peak1D > v;
	Peak1D p;
	
	p.setIntensity(2.5);
	v.push_back(p);

	p.setIntensity(3.5);
	v.push_back(p);

	p.setIntensity(1.5);
	v.push_back(p);
	
	std::sort(v.begin(), v.end(), Peak1D::IntensityLess());
	TEST_REAL_EQUAL(v[0].getIntensity(), 1.5)
	TEST_REAL_EQUAL(v[1].getIntensity(), 2.5)
	TEST_REAL_EQUAL(v[2].getIntensity(), 3.5)

	v[0]=v[2];
	v[2]=p;
	std::sort(v.begin(), v.end(), Peak1D::IntensityLess());
	TEST_REAL_EQUAL(v[0].getIntensity(), 1.5)
	TEST_REAL_EQUAL(v[1].getIntensity(), 2.5)
	TEST_REAL_EQUAL(v[2].getIntensity(), 3.5)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
