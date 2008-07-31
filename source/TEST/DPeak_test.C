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

#include <OpenMS/KERNEL/DPeak.h>

///////////////////////////

START_TEST(DPeak<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

DPeak<10>* d10_ptr = 0;
CHECK((DPeak()))
	d10_ptr = new DPeak<10>;
	TEST_NOT_EQUAL(d10_ptr, 0)
RESULT

CHECK((~DPeak()))
	delete d10_ptr;
RESULT

CHECK((IntensityType getIntensity() const))
	const DPeak<10> p;
	TEST_REAL_EQUAL(p.getIntensity(), 0.0)
RESULT

CHECK((void setIntensity(IntensityTypeintensity)))
	DPeak<3> p;
	TEST_REAL_EQUAL(p.getIntensity(), 0.0)
	p.setIntensity(123.456);
	TEST_REAL_EQUAL(p.getIntensity(), 123.456)
	p.setIntensity(-0.12345);
	TEST_REAL_EQUAL(p.getIntensity(), -0.12345)
	p.setIntensity(0.0);
	TEST_REAL_EQUAL(p.getIntensity(), 0.0)
RESULT

CHECK((const PositionType& getPosition() const))
	const DPeak<10>	p;
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

CHECK((void setPosition(PositionType const& position)))
	DPosition<2> p;
  	p[0] = 876;
  	p[1] = 12345.0;
	DPeak<2> p1;
  	p1.setPosition(p);
	TEST_REAL_EQUAL(p1.getPosition()[0], 876)
	TEST_REAL_EQUAL(p1.getPosition()[1], 12345.0)
RESULT

CHECK((PositionType& getPosition()))
	DPeak<2> p0;
  	p0.getPosition()[0] = 876;
  	p0.getPosition()[1] = 12345.0;
	DPeak<2> p1;
  	p1.setPosition(p0.getPosition());
	TEST_REAL_EQUAL(p1.getPosition()[0], 876)
	TEST_REAL_EQUAL(p1.getPosition()[1], 12345.0)
RESULT

CHECK((DPeak(const DPeak& p)))
	DPeak<3>::PositionType pos;
	pos[0] = 21.21;
	pos[1] = 22.22;
	pos[2] = 23.23;
	DPeak<3> p;
	p.setIntensity(123.456);
	p.setPosition(pos);
	DPeak<3>::PositionType pos2;
	DPeak<3>::IntensityType i2;

	DPeak<3> copy_of_p(p);
	
	i2 = copy_of_p.getIntensity();
	pos2 = copy_of_p.getPosition();
	TEST_REAL_EQUAL(i2, 123.456)

	TEST_REAL_EQUAL(pos2[0], 21.21)
	TEST_REAL_EQUAL(pos2[1], 22.22)
	TEST_REAL_EQUAL(pos2[2], 23.23)	
RESULT

CHECK((DPeak& operator = (const DPeak& rhs)))
	DPeak<3>::PositionType pos;
	pos[0] = 21.21;
	pos[1] = 22.22;
	pos[2] = 23.23;
	DPeak<3> p;
	p.setIntensity(123.456);
	p.setPosition(pos);
	DPeak<3>::PositionType pos2;
	DPeak<3>::IntensityType i2;

	DPeak<3> copy_of_p;
	copy_of_p = p;
	
	i2 = copy_of_p.getIntensity();
	pos2 = copy_of_p.getPosition();
	TEST_REAL_EQUAL(i2, 123.456)

	TEST_REAL_EQUAL(pos2[0], 21.21)
	TEST_REAL_EQUAL(pos2[1], 22.22)
	TEST_REAL_EQUAL(pos2[2], 23.23)	
RESULT

CHECK((bool operator == (const DPeak& rhs) const))
	DPeak<1> p1;
	DPeak<1> p2(p1);
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

CHECK((bool operator != (const DPeak& rhs) const))
	DPeak<1> p1;
	DPeak<1> p2(p1);
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
	std::vector<DPeak<2> > v;
	DPeak<2> p;
	
	p.getPosition()[0]=3.0;
	p.getPosition()[1]=2.5;
	v.push_back(p);

	p.getPosition()[0]=2.0;
	p.getPosition()[1]=3.5;
	v.push_back(p);

	p.getPosition()[0]=1.0;
	p.getPosition()[1]=1.5;
	v.push_back(p);
	
	std::sort(v.begin(), v.end(), DPeak<2>::PositionLess());
	TEST_REAL_EQUAL(v[0].getPosition()[0], 1.0)
	TEST_REAL_EQUAL(v[1].getPosition()[0], 2.0)
	TEST_REAL_EQUAL(v[2].getPosition()[0], 3.0)
	TEST_REAL_EQUAL(v[0].getPosition()[1], 1.5)
	TEST_REAL_EQUAL(v[1].getPosition()[1], 3.5)
	TEST_REAL_EQUAL(v[2].getPosition()[1], 2.5)

	std::sort(v.begin(), v.end(), DPeak<2>::NthPositionLess<1>());
	TEST_REAL_EQUAL(v[0].getPosition()[1], 1.5)
	TEST_REAL_EQUAL(v[1].getPosition()[1], 2.5)
	TEST_REAL_EQUAL(v[2].getPosition()[1], 3.5)
	TEST_REAL_EQUAL(v[0].getPosition()[0], 1.0)
	TEST_REAL_EQUAL(v[1].getPosition()[0], 3.0)
	TEST_REAL_EQUAL(v[2].getPosition()[0], 2.0)
RESULT

CHECK([EXTRA] struct NthPositionLess)

	std::vector<DPeak<2> > v;
	DPeak<2> p;
	
	p.getPosition()[0]=3.0;
	p.getPosition()[1]=2.5;
	v.push_back(p);

	p.getPosition()[0]=2.0;
	p.getPosition()[1]=3.5;
	v.push_back(p);

	p.getPosition()[0]=1.0;
	p.getPosition()[1]=1.5;
	v.push_back(p);

	std::sort(v.begin(), v.end(), DPeak<2>::NthPositionLess<1>());
	TEST_REAL_EQUAL(v[0].getPosition()[1], 1.5)
	TEST_REAL_EQUAL(v[1].getPosition()[1], 2.5)
	TEST_REAL_EQUAL(v[2].getPosition()[1], 3.5)
	TEST_REAL_EQUAL(v[0].getPosition()[0], 1.0)
	TEST_REAL_EQUAL(v[1].getPosition()[0], 3.0)
	TEST_REAL_EQUAL(v[2].getPosition()[0], 2.0)
RESULT

CHECK([EXTRA] struct NthPositionLess<i>::DIMENSION)
	TEST_EQUAL(DPeak<7>::NthPositionLess<5>::DIMENSION, 5)
RESULT

CHECK([EXTRA] struct IntensityLess)
	std::vector<DPeak<2> > v;
	DPeak<2> p;
	
	p.setIntensity(2.5);
	v.push_back(p);

	p.setIntensity(3.5);
	v.push_back(p);

	p.setIntensity(1.5);
	v.push_back(p);
	
	std::sort(v.begin(), v.end(), DPeak<2>::IntensityLess());
	TEST_REAL_EQUAL(v[0].getIntensity(), 1.5)
	TEST_REAL_EQUAL(v[1].getIntensity(), 2.5)
	TEST_REAL_EQUAL(v[2].getIntensity(), 3.5)

	v[0]=v[2];
	v[2]=p;
	std::sort(v.begin(), v.end(), DPeak<2>::IntensityLess());
	TEST_REAL_EQUAL(v[0].getIntensity(), 1.5)
	TEST_REAL_EQUAL(v[1].getIntensity(), 2.5)
	TEST_REAL_EQUAL(v[2].getIntensity(), 3.5)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
