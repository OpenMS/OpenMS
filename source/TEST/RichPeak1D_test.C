// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/KERNEL/RichPeak1D.h>

///////////////////////////

START_TEST(RichPeak1D<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

RichPeak1D* d10_ptr = 0;
START_SECTION((RichPeak1D()))
	d10_ptr = new RichPeak1D;
	TEST_NOT_EQUAL(d10_ptr, 0)
END_SECTION

START_SECTION((~RichPeak1D()))
	delete d10_ptr;
END_SECTION

START_SECTION((RichPeak1D(const RichPeak1D &p)))
	RichPeak1D p;
	p.setIntensity(123.456f);
	p.setMetaValue("cluster_id",4711);
	
	RichPeak1D copy_of_p(p);

	TEST_REAL_SIMILAR(copy_of_p.getIntensity(), 123.456)
	TEST_EQUAL(copy_of_p.getMetaValue("cluster_id"),DataValue(4711));
END_SECTION

START_SECTION((RichPeak1D& operator=(const RichPeak1D &rhs)))
	RichPeak1D p;
	p.setIntensity(123.456f);
	p.setMetaValue("cluster_id",4711);
	
	RichPeak1D copy_of_p;
	copy_of_p = p;

	TEST_REAL_SIMILAR(copy_of_p.getIntensity(), 123.456)
	TEST_EQUAL(copy_of_p.getMetaValue("cluster_id"),DataValue(4711));
END_SECTION

START_SECTION((bool operator == (const RichPeak1D& rhs) const))
	RichPeak1D p1, p2;
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

START_SECTION((bool operator != (const RichPeak1D& rhs) const))
	RichPeak1D p1, p2;
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
	RichPeak1D p;
	p.setMetaValue(2,String("bla"));
 	RichPeak1D p2(p);
	TEST_EQUAL(p.getMetaValue(2), "bla")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
 	p.setMetaValue(2,String("bluff"));
	TEST_EQUAL(p.getMetaValue(2), "bluff")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
END_SECTION

START_SECTION(([EXTRA] meta info with assignment))
	RichPeak1D p;
	p.setMetaValue(2,String("bla"));
 	RichPeak1D p2 = p;
	TEST_EQUAL(p.getMetaValue(2), "bla")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
 	p.setMetaValue(2,String("bluff"));
	TEST_EQUAL(p.getMetaValue(2), "bluff")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
