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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/KERNEL/RichPeak2D.h>

///////////////////////////

START_TEST(RichPeak2D<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

RichPeak2D* d10_ptr = 0;
RichPeak2D* d10_nullPointer = 0;
START_SECTION((RichPeak2D()))
	d10_ptr = new RichPeak2D;
  TEST_NOT_EQUAL(d10_ptr, d10_nullPointer)
END_SECTION

START_SECTION((~RichPeak2D()))
	delete d10_ptr;
END_SECTION

START_SECTION((RichPeak2D(const RichPeak2D &p)))
	RichPeak2D p;
	p.setIntensity(123.456f);
	p.setMetaValue("cluster_id",4711);
	
	RichPeak2D copy_of_p(p);

	TEST_REAL_SIMILAR(copy_of_p.getIntensity(), 123.456f)
	TEST_EQUAL(copy_of_p.getMetaValue("cluster_id"),DataValue(4711));
END_SECTION

START_SECTION((RichPeak2D(const Peak2D &p)))
	Peak2D p;
	p.setIntensity(123.456f);
	
	RichPeak2D copy_of_p(p);

	TEST_REAL_SIMILAR(copy_of_p.getIntensity(), 123.456f)
END_SECTION		
		
START_SECTION((RichPeak2D& operator=(const RichPeak2D &rhs)))
	RichPeak2D p;
	p.setIntensity(123.456f);
	p.setMetaValue("cluster_id",4711);
	
	RichPeak2D copy_of_p;
	copy_of_p = p;

	TEST_REAL_SIMILAR(copy_of_p.getIntensity(), 123.456f)
	TEST_EQUAL(copy_of_p.getMetaValue("cluster_id"),DataValue(4711));
END_SECTION
		
START_SECTION((RichPeak2D& operator=(const Peak2D &rhs)))
	Peak2D p;
	p.setIntensity(123.456f);
	
	RichPeak2D copy_of_p;
	copy_of_p.setMetaValue("cluster_id",4711);
	copy_of_p = p;

	TEST_REAL_SIMILAR(copy_of_p.getIntensity(), 123.456f)
	TEST_EQUAL(copy_of_p.isMetaEmpty(), true);
END_SECTION
		
START_SECTION((bool operator == (const RichPeak2D& rhs) const))
	RichPeak2D p1, p2;
	TEST_EQUAL(p1==p2, true)
	
	p1.setIntensity(5.0f);
	TEST_EQUAL(p1==p2, false)
	p2.setIntensity(5.0f);
	TEST_EQUAL(p1==p2, true)

	p1.setMetaValue("cluster_id",4711);
	TEST_EQUAL(p1==p2, false)
	p1.removeMetaValue("cluster_id");
	TEST_EQUAL(p1==p2, true)		
END_SECTION

START_SECTION((bool operator != (const RichPeak2D& rhs) const))
	RichPeak2D p1, p2;
	TEST_EQUAL(p1!=p2, false)
	
	p1.setIntensity(5.0f);
	TEST_EQUAL(p1!=p2, true)
	p2.setIntensity(5.0f);
	TEST_EQUAL(p1!=p2, false)

	p1.setMetaValue("cluster_id",4711);
	TEST_EQUAL(p1!=p2, true)
	p1.removeMetaValue("cluster_id");
	TEST_EQUAL(p1!=p2, false)	
END_SECTION

START_SECTION(([EXTRA] meta info with copy constructor))
	RichPeak2D p;
	p.setMetaValue(2,String("bla"));
 	RichPeak2D p2(p);
	TEST_EQUAL(p.getMetaValue(2), "bla")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
 	p.setMetaValue(2,String("bluff"));
	TEST_EQUAL(p.getMetaValue(2), "bluff")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
END_SECTION

START_SECTION(([EXTRA] meta info with assignment))
	RichPeak2D p;
	p.setMetaValue(2,String("bla"));
 	RichPeak2D p2 = p;
	TEST_EQUAL(p.getMetaValue(2), "bla")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
 	p.setMetaValue(2,String("bluff"));
	TEST_EQUAL(p.getMetaValue(2), "bluff")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
