// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: DPeak_test.C,v 1.4 2006/03/28 08:03:34 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/KERNEL/DPeak.h>

///////////////////////////

START_TEST(DPeak<D>, "$Id: DPeak_test.C,v 1.4 2006/03/28 08:03:34 marc_sturm Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

DPeak<10>* d10_ptr = 0;
CHECK(DPeak())
	d10_ptr = new DPeak<10>;
	TEST_NOT_EQUAL(d10_ptr, 0)
RESULT

CHECK(~DPeak())
	delete d10_ptr;
RESULT

CHECK(const IntensityType& getIntensity() const)
	const DPeak<10> p;
	TEST_REAL_EQUAL(p.getIntensity(), 0.0)
RESULT

CHECK(const PositionType& getPosition() const)
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

CHECK(IntensityType& getIntensity())
	DPeak<3> p;
	TEST_REAL_EQUAL(p.getIntensity(), 0.0)
	p.getIntensity() = 123.456;
	TEST_REAL_EQUAL(p.getIntensity(), 123.456)
	p.getIntensity() = -0.12345;
	TEST_REAL_EQUAL(p.getIntensity(), -0.12345)
	p.getIntensity() = 0.0;
	TEST_REAL_EQUAL(p.getIntensity(), 0.0)
RESULT

CHECK(PositionType& getPosition())
	DPeak<3>::PositionType pos;
	DPeak<3> p;
	pos = p.getPosition();
	TEST_REAL_EQUAL(pos[0], 0.0)
	TEST_REAL_EQUAL(pos[1], 0.0)
	TEST_REAL_EQUAL(pos[2], 0.0)
	pos[0] = 1.0;
	pos[1] = 2.0;
	pos[2] = 3.0;
	p.getPosition() = pos;
	DPeak<3>::PositionType pos2(p.getPosition());
	TEST_REAL_EQUAL(pos2[0], 1.0)
	TEST_REAL_EQUAL(pos2[1], 2.0)
	TEST_REAL_EQUAL(pos2[2], 3.0)	
RESULT

CHECK(DPeak(const DPeak<D>& p))
	DPeak<3>::PositionType pos;
	pos[0] = 21.21;
	pos[1] = 22.22;
	pos[2] = 23.23;
	DPeak<3> p;
	p.getIntensity() = 123.456;
	p.getPosition() = pos;
	p.setMetaValue("cluster_id",4711);
	
	DPeak<3>::PositionType pos2;
	DPeak<3>::IntensityType i2;

	DPeak<3> copy_of_p(p);
	i2 = copy_of_p.getIntensity();
	pos2 = copy_of_p.getPosition();

	TEST_REAL_EQUAL(i2, 123.456)

	TEST_REAL_EQUAL(pos2[0], 21.21)
	TEST_REAL_EQUAL(pos2[1], 22.22)
	TEST_REAL_EQUAL(pos2[2], 23.23)	
	
	TEST_EQUAL(p.getMetaValue("cluster_id"),DataValue(4711));
RESULT

CHECK(DPeak& operator = (const DPeak& rhs))
DPeak<3>::PositionType pos;
	pos[0] = 21.21;
	pos[1] = 22.22;
	pos[2] = 23.23;
	DPeak<3> p;
	p.getIntensity() = 123.456;
	p.getPosition() = pos;
	p.setMetaValue("cluster_id",4712);
	
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

CHECK(DPeak& operator == (const DPeak& rhs))
	DPeak<1> p1;
	DPeak<1> p2(p1);
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

CHECK(DPeak& operator != (const DPeak& rhs))
	DPeak<1> p1;
	DPeak<1> p2(p1);
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

CHECK(meta info with copy constructor)
	DPeak<1> p;
	p.setMetaValue(2,std::string("bla"));
 	DPeak<1> p2(p);
	TEST_EQUAL(p.getMetaValue(2), "bla")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
 	p.setMetaValue(2,std::string("bluff"));
	TEST_EQUAL(p.getMetaValue(2), "bluff")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
RESULT

CHECK(meta info with assignment)
	DPeak<1> p;
	p.setMetaValue(2,std::string("bla"));
 	DPeak<1> p2 = p;
	TEST_EQUAL(p.getMetaValue(2), "bla")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
 	p.setMetaValue(2,std::string("bluff"));
	TEST_EQUAL(p.getMetaValue(2), "bluff")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
