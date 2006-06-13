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
// $Id: Histogram_test.C,v 1.3 2006/03/28 08:03:34 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/MATH/STATISTICS/Histogram.h>
#include <iostream>
#include <vector>

using namespace OpenMS;
using namespace OpenMS::Math;
using namespace std;

///////////////////////////

START_TEST(Distribution, "$Id: Histogram_test.C,v 1.3 2006/03/28 08:03:34 marc_sturm Exp $")

/////////////////////////////////////////////////////////////

Histogram<float,float>* dis_ptr = 0;

CHECK( Histogram())
	dis_ptr = new Histogram<float,float>();
	TEST_NOT_EQUAL(dis_ptr, 0)
RESULT

CHECK( ~Histogram() )
	delete dis_ptr;
RESULT

Histogram<float,float> d(0, 10, 1);

CHECK( Histogram(const Histogram&))
	Histogram<float, float> d2(d);
	TEST_EQUAL(d == d2, true)
RESULT

CHECK( min() )
	TEST_REAL_EQUAL(d.min(), 0.0)
RESULT

CHECK( max() ) 
	TEST_REAL_EQUAL(d.max(), 10.0)
RESULT

CHECK( binSize() )
	TEST_REAL_EQUAL(d.binSize(), 1)
RESULT

CHECK( size() )
	TEST_EQUAL(d.size(), 10)
RESULT

CHECK( inc(double val, double inrement) )
	d.inc(5, 250.3);
	TEST_REAL_EQUAL(d.maxValue(), 250.3)
	TEST_EXCEPTION(Exception::OutOfRange, d.inc(10.1, 250.3))
RESULT

CHECK( minValue() )
	TEST_REAL_EQUAL(d.minValue(), 0.0)
RESULT

CHECK( maxValue() )
	TEST_REAL_EQUAL(d.maxValue(), 250.3)
RESULT

CHECK( binValue(UnsignedInt index) )
	TEST_REAL_EQUAL(d.binValue(UnsignedInt(5)), 250.3)
	TEST_EXCEPTION(Exception::OutOfRange, d.binValue(UnsignedInt(11)))
RESULT
	
CHECK( binValue(double val) )
	TEST_REAL_EQUAL(d.binValue(5.5f), 250.3)
	TEST_EXCEPTION(Exception::OutOfRange, d.binValue(10.1f))
RESULT
	
CHECK( set(double min, double max, double bin_size) )
	d.set(1, 11, 2);
	TEST_REAL_EQUAL(d.min(), 1)
	TEST_REAL_EQUAL(d.max(), 11)
	TEST_REAL_EQUAL(d.size(), 5)
	TEST_REAL_EQUAL(d.binSize(), 2)
RESULT

CHECK(bool operator == )
	Histogram<float, float> dist(1, 11, 2);
	TEST_EQUAL(d == dist, true)
RESULT

CHECK(bool operator != )
	Histogram<float, float> dist(1, 12, 2);
	TEST_EQUAL(d != dist, true)
RESULT

CHECK(Histogram& operator = )
	Histogram<float, float> dist;
	dist = d;
	TEST_EQUAL(d == dist, true)
RESULT
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
