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

#include <OpenMS/MATH/STATISTICS/Histogram.h>
#include <iostream>
#include <vector>

using namespace OpenMS;
using namespace OpenMS::Math;
using namespace std;

///////////////////////////

START_TEST(Histogram, "$Id$")

/////////////////////////////////////////////////////////////

Histogram<float,float>* dis_ptr = 0;

CHECK(Histogram())
	dis_ptr = new Histogram<float,float>();
	TEST_NOT_EQUAL(dis_ptr, 0)
RESULT

CHECK(~Histogram())
	delete dis_ptr;
RESULT

Histogram<float,float> d(0, 10, 1);

CHECK(Histogram(const Histogram& histogram))
	Histogram<float, float> d2(d);
	TEST_EQUAL(d == d2, true)
RESULT

CHECK(BinSizeType min() const)
	TEST_REAL_EQUAL(d.min(), 0.0)
RESULT

CHECK(BinSizeType max() const)
	TEST_REAL_EQUAL(d.max(), 10.0)
RESULT

CHECK(BinSizeType binSize() const)
	TEST_REAL_EQUAL(d.binSize(), 1)
RESULT

CHECK(Size size() const)
	TEST_EQUAL(d.size(), 10)
RESULT

CHECK(Histogram(BinSizeType min, BinSizeType max, BinSizeType bin_size) throw(Exception::OutOfRange))
	Histogram<float, float> d3(5.5, 7.7, 0.2);
	TEST_REAL_EQUAL(d3.min(), 5.5)
	TEST_REAL_EQUAL(d3.max(), 7.7)
	TEST_REAL_EQUAL(d3.binSize(), 0.2)
RESULT

CHECK(ValueType minValue() const)
	TEST_REAL_EQUAL(d.minValue(), 0.0)
RESULT

CHECK(ValueType maxValue() const)
	TEST_REAL_EQUAL(d.maxValue(), 0.0)
RESULT

CHECK(ValueType operator [] (UnsignedInt index) const throw(Exception::IndexOverflow))
	d.set(4, 14, 2);
	TEST_EQUAL(d.size(),5);
	TEST_REAL_EQUAL(d[0],0.0);
	TEST_REAL_EQUAL(d[1],0.0);
	TEST_REAL_EQUAL(d[2],0.0);
	TEST_REAL_EQUAL(d[3],0.0);
	TEST_REAL_EQUAL(d[4],0.0);
	TEST_EXCEPTION(Exception::IndexOverflow, d[5])
RESULT

CHECK(void inc(BinSizeType val, ValueType increment=1) throw(Exception::OutOfRange))
	TEST_EXCEPTION(Exception::OutOfRange, d.inc(3.9, 250.3))
	TEST_EXCEPTION(Exception::OutOfRange, d.inc(14.1, 250.3))
	d.inc(4, 1.0);
	d.inc(5.9, 1.0);
	TEST_REAL_EQUAL(d[0],2.0);
	TEST_REAL_EQUAL(d[1],0.0);
	TEST_REAL_EQUAL(d[2],0.0);
	TEST_REAL_EQUAL(d[3],0.0);
	TEST_REAL_EQUAL(d[4],0.0);
	
	d.inc(8.0, 45.0);
	d.inc(8.1, 1.0);
	d.inc(9.9, 4.0);

	TEST_REAL_EQUAL(d[0],2.0);
	TEST_REAL_EQUAL(d[1],0.0);
	TEST_REAL_EQUAL(d[2],50.0);
	TEST_REAL_EQUAL(d[3],0.0);
	TEST_REAL_EQUAL(d[4],0.0);

	d.inc(12.0, 1.0);
	d.inc(13.1, 2.0);
	d.inc(14.0, 3.0);	

	TEST_REAL_EQUAL(d[0],2.0);
	TEST_REAL_EQUAL(d[1],0.0);
	TEST_REAL_EQUAL(d[2],50.0);
	TEST_REAL_EQUAL(d[3],0.0);
	TEST_REAL_EQUAL(d[4],6.0);
RESULT

CHECK(ConstIterator begin() const)
	Histogram<float,float>::ConstIterator it = d.begin();
	TEST_REAL_EQUAL(*it, 2.0)
RESULT

CHECK(ConstIterator end() const)
	Histogram<float,float>::ConstIterator it = d.begin();
	TEST_REAL_EQUAL(*it,2.0);
	++it;
	TEST_REAL_EQUAL(*it,0.0);
	++it;
	TEST_REAL_EQUAL(*it,50.0);
	++it;
	TEST_REAL_EQUAL(*it,0.0);
	++it;
	TEST_REAL_EQUAL(*it,6.0);
	++it;
	TEST_EQUAL(it==d.end(),true);
RESULT

CHECK(ValueType binValue(BinSizeType val) const throw(Exception::OutOfRange))
	TEST_EXCEPTION(Exception::OutOfRange, d.binValue(3.9))
	TEST_REAL_EQUAL(d.binValue(4.0),2.0);
	TEST_REAL_EQUAL(d.binValue(5.9),2.0);
	TEST_REAL_EQUAL(d.binValue(6.0),0.0);
	TEST_REAL_EQUAL(d.binValue(7.9),0.0);
	TEST_REAL_EQUAL(d.binValue(8.0),50.0);
	TEST_REAL_EQUAL(d.binValue(9.9),50.0);
	TEST_REAL_EQUAL(d.binValue(10.0),0.0);
	TEST_REAL_EQUAL(d.binValue(11.9),0.0);
	TEST_REAL_EQUAL(d.binValue(12.0),6.0);
	TEST_REAL_EQUAL(d.binValue(14.0),6.0);
	TEST_EXCEPTION(Exception::OutOfRange, d.binValue(14.1))
RESULT
	
CHECK(void set(BinSizeType min, BinSizeType max, BinSizeType bin_size) throw(Exception::OutOfRange))
	d.set(1, 11, 2);
	TEST_REAL_EQUAL(d.min(), 1)
	TEST_REAL_EQUAL(d.max(), 11)
	TEST_REAL_EQUAL(d.size(), 5)
	TEST_REAL_EQUAL(d.binSize(), 2)
RESULT

CHECK(bool operator == (const Histogram& histogram) const)
	Histogram<float, float> dist(1, 11, 2);
	TEST_EQUAL(d == dist, true)
RESULT

CHECK(bool operator != (const Histogram& histogram) const)
	Histogram<float, float> dist(1, 12, 2);
	TEST_EQUAL(d != dist, true)
RESULT

CHECK(Histogram& operator = (const Histogram& histogram))
	Histogram<float, float> dist;
	dist = d;
	TEST_EQUAL(d == dist, true)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
