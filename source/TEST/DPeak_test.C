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

CHECK((DPeak(DPeak const& p)))
	DPeak<3> p;
	p.getIntensity() = 123.456;
	p.setMetaValue("cluster_id",4711);
	
	DPeak<3> copy_of_p(p);

	TEST_REAL_EQUAL(copy_of_p.getIntensity(), 123.456)
	TEST_EQUAL(copy_of_p.getMetaValue("cluster_id"),DataValue(4711));
RESULT

CHECK((DPeak& operator = (const DPeak& rhs)))
	DPeak<3> p;
	p.getIntensity() = 123.456;
	p.setMetaValue("cluster_id",4711);
	
	DPeak<3> copy_of_p;
	copy_of_p = p;

	TEST_REAL_EQUAL(copy_of_p.getIntensity(), 123.456)
	TEST_EQUAL(copy_of_p.getMetaValue("cluster_id"),DataValue(4711));
RESULT

CHECK((bool operator == (const DPeak& rhs) const))
	DPeak<1> p1, p2;
	TEST_REAL_EQUAL(p1==p2, true)
	
	p1.getIntensity()=5;
	TEST_REAL_EQUAL(p1==p2, false)
	p2.getIntensity()=5;
	TEST_REAL_EQUAL(p1==p2, true)

	p1.setMetaValue("cluster_id",4711);
	TEST_REAL_EQUAL(p1==p2, false)
	p1.removeMetaValue("cluster_id");
	TEST_REAL_EQUAL(p1==p2, true)		
RESULT

CHECK((bool operator != (const DPeak& rhs) const))
	DPeak<1> p1, p2;
	TEST_REAL_EQUAL(p1!=p2, false)
	
	p1.getIntensity()=5;
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.getIntensity()=5;
	TEST_REAL_EQUAL(p1!=p2, false)

	p1.setMetaValue("cluster_id",4711);
	TEST_REAL_EQUAL(p1!=p2, true)
	p1.removeMetaValue("cluster_id");
	TEST_REAL_EQUAL(p1!=p2, false)	
RESULT

CHECK(([EXTRA] meta info with copy constructor))
	DPeak<1> p;
	p.setMetaValue(2,std::string("bla"));
 	DPeak<1> p2(p);
	TEST_EQUAL(p.getMetaValue(2), "bla")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
 	p.setMetaValue(2,std::string("bluff"));
	TEST_EQUAL(p.getMetaValue(2), "bluff")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
RESULT

CHECK(([EXTRA] meta info with assignment))
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
